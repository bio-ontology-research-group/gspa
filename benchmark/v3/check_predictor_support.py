"""Real EC-claim CLI replay plus deliberate identity/range faults, not accuracy truth."""
import argparse
from collections import Counter
import copy
import json
from pathlib import Path
import shutil
import subprocess
import sys
import numpy as np

from check_model_input_support import predictions
from migrate_cohort_feature_contract import hashes, check_hashes

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
import train_context_model as T


def save(path,value):
    path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')


def run(before,after,original_models,models,claims,layout,out):
    out.mkdir(parents=True,exist_ok=False)
    java=Path(shutil.which('java')).resolve()
    paths=[before,after,claims,layout,java,Path(__file__),Path(T.__file__),
           Path(__file__).with_name('check_model_input_support.py')]
    paths+=list(original_models.glob('model_*.json'))+list(models.glob('model_*.json'))
    identity=hashes(paths);save(out/'inputs.json',identity)
    original=[r for r in map(json.loads,claims.read_text().splitlines()) if r['function_type']=='EC']
    if not original:raise ValueError('no EC inputs')
    cases={};summaries={};runs=[];baseline_features={};old_scores={}
    faults=['verified','tool_version','reference_sha256','score_model','missing','weaker_missing','zero_scores','legacy_support']
    for case in faults:
        rows=copy.deepcopy(original)
        if case in ('tool_version','reference_sha256','score_model'):
            for row in rows:row['metadata']['predictor_contract'][case]='e'*64 if case.endswith('_sha256') else 'deliberate-test-fault'
        if case=='missing':
            for row in rows:row['metadata'].pop('predictor_contract')
        if case=='weaker_missing':
            weaker=copy.deepcopy(rows)
            for row in weaker:row['metadata'].pop('predictor_contract');row['raw_score']=0.
            rows+=weaker
        if case=='zero_scores':
            for row in rows:row['raw_score']=0.
        p=out/(case+'.claims.jsonl');p.write_text(''.join(json.dumps(r,sort_keys=True)+'\n' for r in rows))
        cases[case]=p;identity.update(hashes([p]))
    save(out/'inputs.json',identity)
    for case,p in cases.items():
        baseline=None
        for engine,jar in [('baseline',after),('before',before),('after',after)]:
            for kind in (['baseline'] if engine=='baseline' else ['logistic','gbdt']):
                directory=out/case/(engine+'-'+kind);directory.mkdir(parents=True)
                output=directory/'predictions.tsv';features=directory/'features'
                cmd=[str(java),'-Xmx2g','-jar',str(jar),'integrate','--claims',str(p),'--out',str(output),
                     '--lite','--enable-priors','none','--genome-id','GCA_000007005.1','--genome-layout',str(layout),
                     '--feature-export-dir',str(features)]
                model_path=(original_models if engine=='before' or case=='legacy_support' else models)/f'model_{kind}.json'
                if engine!='baseline':cmd+=['--learned-model',str(model_path)]
                with (directory/'run.log').open('w') as log:
                    result=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,timeout=180)
                runs.append(dict(case=case,engine=engine,kind=kind,command=cmd,exit_code=result.returncode))
                save(out/'runs.json',runs)
                if result.returncode:raise ValueError('CLI failed: '+str(directory))
                ledger=json.loads(Path(str(output)+'.stages.json').read_text())
                stages={s['name']:s for s in ledger['stages']}
                if any(s['status'] in ('failed','not_run') for s in stages.values()) or any(
                        stages[k]['status']!='executed' for k in ('reference_preparation','evidence_integration')):
                    raise ValueError('incomplete stage ledger')
                fs=[json.loads(s) for s in (features/'candidate_features.jsonl').read_text().splitlines()]
                projection={r['candidate_key']:r['features'] for r in fs}
                if case=='verified':baseline_features[engine+'-'+kind]=projection
                elif case!='zero_scores' and projection!=baseline_features[engine+'-'+kind]:
                    raise ValueError('identity-only fault changed numerical features')
                current=predictions(output)
                if engine=='baseline':baseline=current;continue
                report=json.loads(Path(str(output)+'.learned.json').read_text())
                if report['model']['sha256']!=T.sha256_file(model_path):raise ValueError('wrong model identity')
                records={r['candidate_key']:r for r in report['candidates']}
                if set(records)!=set(projection) or len(records)!=len(fs):raise ValueError('candidate coverage changed')
                if engine=='before':
                    scores={k:r['probability'] for k,r in records.items()}
                    if not any(v is not None for v in scores.values()):raise ValueError('before control has no supported candidates')
                    if case=='verified':old_scores[kind]=scores
                    elif case!='zero_scores' and scores!=old_scores[kind]:raise ValueError('before identity fault changed scores')
                    summaries[case+'-'+engine+'-'+kind]=dict(scored=sum(v is not None for v in scores.values()),
                        reasons=dict(Counter(r['unsupported_reason'] or 'scored' for r in records.values())));continue
                art=json.loads(model_path.read_text());reasons=Counter();error=0.
                for f in fs:
                    r=records[f['candidate_key']];x=np.array([f['features'][n] for n in art['contract']['features']])
                    reason=T.unsupported_input_reason(art,x,f['sequence_sources'],f.get('quality_evidence'),
                        f['cohort_context_mode'],f['sequence_contracts'])
                    if reason!=r['unsupported_reason'] or r['sequence_contracts']!=f['sequence_contracts']:
                        raise ValueError('Python/JVM support or provenance differs')
                    reasons[reason or 'scored']+=1
                    if reason:
                        if r['score_status']!='unsupported' or any(r[k] is not None for k in (
                                'probability','margin','above_display_threshold','locus_probability','locus_null_probability')) or f['model']['probability'] is not None:
                            raise ValueError('unsupported candidate received probability')
                        key='|'.join([f['protein_id'],f['function_type'],f['function_id']])
                        if key not in baseline or key not in current or current[key]!=baseline[key]:
                            raise ValueError('unsupported legacy prediction changed')
                    else:
                        prob=T.predict_proba(art,x,f['sequence_sources'],f.get('quality_evidence'),
                            f['cohort_context_mode'],f['sequence_contracts'])
                        error=max(error,abs(prob-r['probability']),abs(prob-f['model']['probability']))
                        if abs(prob-old_scores[kind][f['candidate_key']])>1e-12:
                            raise ValueError('preserved scoring parameters changed a supported probability')
                if error>1e-12:raise ValueError('Python/JVM probabilities differ')
                expected=('unseen_predictor_contract' if case in ('tool_version','reference_sha256','score_model') else
                    'candidate_predictor_contract_unavailable' if case in ('missing','weaker_missing') else
                    'sequence_score_outside_training_calibration_range' if case=='zero_scores' else
                    'model_predictor_contract_support_unavailable' if case=='legacy_support' else None)
                if expected:
                    previously_supported=sum(v is not None for v in old_scores[kind].values())
                    if reasons[expected]<previously_supported or any(k not in (expected,
                            'unseen_source_function_role_or_availability_profile') for k in reasons) or reasons.get(
                            'unseen_source_function_role_or_availability_profile',0)>len(fs)-previously_supported:
                        raise ValueError('fault did not abstain as expected: '+str(reasons))
                if case=='verified' and not reasons['scored']:raise ValueError('no positive scoring control')
                summaries[case+'-'+engine+'-'+kind]=dict(reasons=dict(reasons),max_probability_error=error)
    check_hashes(identity)
    result=dict(status='verified',scope=__doc__,input_claims=len(original),cli_runs=len(runs),inputs=identity,
                summaries=summaries,outputs=hashes([p for p in out.rglob('*') if p.is_file()]))
    save(out/'verification.json',result)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('before','after','original-models','models','claims','layout','out'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();r=run(*[getattr(a,k).resolve() for k in ('before','after','original_models','models','claims','layout','out')])
    print(json.dumps({k:r[k] for k in ('status','input_claims','cli_runs')}))
