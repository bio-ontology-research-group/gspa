"""Descriptive matched comparison of verified production outputs on frozen supervision.

This previously examined cohort measures reference-annotation support under its
closed-world policy. Unconfirmed rows are not experimentally proven negatives.
All test rows, including learned abstentions with rule fallback, remain included.
No model promotion, biological calibration, unseen-family or null claim is made.
"""
import argparse
from collections import Counter
import csv
import json
import math
from pathlib import Path
import sys

import fit_observed_gap_models as B
import v3_metrics as M
F=B.F
ARMS=('rule_no_gaps','rule_observed_gaps','logistic','gbdt')
POLICY=dict(schema='gspa.observed-gap-comparison/1',scope=__doc__,arms=list(ARMS),
    thresholds='Separate max-F1 thresholds from complete calibration rows; ties choose higher threshold. Also report common fixed 0.5.',
    population='All frozen test supervision rows within each transfer regime; no supported-only selection or cross-regime pooling.',
    bootstrap=dict(unit='paired genome',replicates=10000,seed=20261027),
    uncertainty='Descriptive percentile resampling of four previously examined test genomes; not a calibrated population confidence claim.')


def save(path,value):
    path.write_text(json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n')


def probability(value):
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or not 0<=value<=1:
        raise ValueError('invalid probability')
    return value


def function_key(row):
    term=row['function_id'];kind=row['function_type']
    if kind=='EC' and term.startswith('EC:'):term=term[3:]
    return row['protein_id'],kind,term


def read_scores(path):
    result={}
    with path.open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            key=function_key(dict(row,function_type=row['type']))
            if key in result:raise ValueError('duplicate prediction identity')
            result[key]=probability(float(row['posterior_prob']))
    return result


def attach_scores(rows, by_arm):
    """Require all registered supervision identities, never silently intersect."""
    result=[];seen=set()
    if set(by_arm)!=set(ARMS):raise ValueError('missing comparison arm')
    for row in rows:
        if row['split'] not in ('train','calibration','test'):raise ValueError('unknown split')
        if type(row['label']) is not int or row['label'] not in (0,1):raise ValueError('invalid label')
        gid=row['genome_id'];key=function_key(row);identity=(gid,key)
        if identity in seen:raise ValueError('duplicate supervised function identity')
        seen.add(identity)
        scores={}
        for arm in ARMS:
            if gid not in by_arm[arm] or key not in by_arm[arm][gid]:raise ValueError('missing supervised prediction')
            scores[arm]=probability(by_arm[arm][gid][key])
        result.append(dict(example_id=row['example_id'],genome_id=gid,split=row['split'],
            label=row['label'],label_basis=row['label_basis'],scores=scores))
    return result


def select_thresholds(calibration):
    if not calibration or any(r['split']!='calibration' for r in calibration):
        raise ValueError('threshold selection requires calibration rows only')
    if {r['label'] for r in calibration}!={0,1}:raise ValueError('both calibration classes required')
    return {arm:F.T.choose_threshold(F.T.np.asarray([r['scores'][arm] for r in calibration]),
        F.T.np.asarray([r['label'] for r in calibration]),'max_f1') for arm in ARMS}


def metrics(rows,arm,threshold):
    if not rows:raise ValueError('empty evaluation population')
    bins=[[] for _ in range(10)];tp=fp=fn=tn=0;loss=0.
    for r in rows:
        p=probability(r['scores'][arm]);y=r['label'];call=p>=threshold
        tp+=call and y==1;fp+=call and y==0;fn+=not call and y==1;tn+=not call and y==0
        loss+=(p-y)**2;bins[min(int(p*10),9)].append((p,y))
    return dict(n=len(rows),positives=tp+fn,tp=tp,fp=fp,fn=fn,tn=tn,
        precision=tp/(tp+fp) if tp+fp else None,recall=tp/(tp+fn) if tp+fn else None,
        f1=2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else None,
        brier=loss/len(rows),ece_10bin=sum(abs(sum(p-y for p,y in b)) for b in bins)/len(rows),
        threshold=threshold)


def paired(rows,arm,base,thresholds,n_boot=10000):
    vectors={};counts=Counter()
    for gid in sorted({r['genome_id'] for r in rows}):
        subset=[r for r in rows if r['genome_id']==gid];corrected=harmed=changed=0
        for r in subset:
            a=r['scores'][arm]>=thresholds[arm];b=r['scores'][base]>=thresholds[base]
            corrected+=a==r['label'] and b!=r['label'];harmed+=a!=r['label'] and b==r['label']
            changed+=a!=b
        delta=metrics(subset,arm,thresholds[arm])['brier']-metrics(subset,base,thresholds[base])['brier']
        vectors[gid]=[len(subset),delta,1,corrected,harmed,changed]
        counts.update(corrected=corrected,harmed=harmed,changed=changed)
    boot=M.genome_bootstrap(vectors,lambda v:v[1]/v[2],n_boot,POLICY['bootstrap']['seed'])
    # Undefined intervals are explicit nulls, never nonstandard JSON NaN.
    boot={k:(None if isinstance(v,float) and not math.isfinite(v) else v) for k,v in boot.items()}
    return dict(reference_corrected=counts['corrected'],reference_harmed=counts['harmed'],
        net_reference_corrections=counts['corrected']-counts['harmed'],changed_calls=counts['changed'],
        mean_genome_brier_delta=boot,per_genome={g:dict(rows=v[0],brier_delta=v[1],
        reference_corrected=v[3],reference_harmed=v[4],changed_calls=v[5]) for g,v in vectors.items()})


def evaluate_rows(test,thresholds,n_boot=10000):
    if not test or any(r['split']!='test' for r in test):raise ValueError('test rows only')
    result=dict(rows=len(test),genomes=sorted({r['genome_id'] for r in test}),
        label_basis=dict(Counter(r['label_basis'] for r in test)),policies={})
    for name,cuts in [('calibration_selected',{a:thresholds[a]['value'] for a in ARMS}),('fixed_0.5',{a:.5 for a in ARMS})]:
        result['policies'][name]=dict(arms={a:metrics(test,a,cuts[a]) for a in ARMS},comparisons={
            f'{arm}_vs_{base}':paired(test,arm,base,cuts,n_boot)
            for base in ('rule_no_gaps','rule_observed_gaps') for arm in ARMS if arm!=base})
    return result


def register(training,replay,out):
    fitting=B.verify(training);comparison=Path(fitting['comparison']);protocol=B.G.verify(comparison)
    production=json.loads((replay/'verification.json').read_text())
    if production.get('status')!='verified':raise ValueError('complete verified production replay required')
    F.check_hashes(production['inputs']);F.check_hashes(production['outputs'])
    export=json.loads((comparison/'verification.json').read_text());F.check_hashes(export['outputs'])
    sources=[Path(m.__file__).resolve() for m in tuple(sys.modules.values()) if getattr(m,'__file__',None)
        and Path(m.__file__).suffix=='.py' and Path(m.__file__).resolve().is_relative_to(F.ROOT)]
    inputs=F.hashes(sources+[replay/'verification.json',comparison/'verification.json',training/'training.json'])
    inputs.update(production['inputs']);inputs.update(production['outputs']);inputs.update(export['outputs'])
    plans={}
    for mode in B.G.MODES:
        table=training/mode/'default/table/training_table.jsonl';inputs.update(F.hashes([table]));paths={}
        for arm in ARMS:
            if arm.startswith('rule_'):
                manifest=Path(protocol['plans'][mode][arm[5:]]['manifest'])
                with manifest.open() as f:paths[arm]={r['genome_id']:r['out_path'] for r in csv.DictReader(f,delimiter='\t')}
            else:paths[arm]={gid:str(replay/mode/arm/'inference'/(gid+'.integrated.tsv')) for gid in protocol['genomes']}
            if set(paths[arm])!=set(protocol['genomes']):raise ValueError('incomplete genome coverage')
        plans[mode]=dict(table=str(table),paths=paths)
    out.mkdir(parents=True,exist_ok=False);save(out/'protocol.json',dict(POLICY,plans=plans))
    save(out/'registration.json',dict(status='registered',registered_at=F.datetime.now(F.timezone.utc).isoformat(),
        protocol_sha256=F.T.sha256_file(out/'protocol.json'),inputs=inputs))


def run(out):
    reg=json.loads((out/'registration.json').read_text());F.check_hashes(reg['inputs'])
    if F.T.sha256_file(out/'protocol.json')!=reg['protocol_sha256']:raise ValueError('comparison protocol changed')
    protocol=json.loads((out/'protocol.json').read_text())
    if any(protocol[k]!=v for k,v in POLICY.items()):raise ValueError('analysis policy changed')
    if (out/'thresholds.json').exists():raise FileExistsError('refusing analysis overwrite')
    prepared={};thresholds={}
    for mode,plan in protocol['plans'].items():
        with Path(plan['table']).open() as f:rows=[json.loads(s) for s in f]
        scores={a:{g:read_scores(Path(p)) for g,p in paths.items()} for a,paths in plan['paths'].items()}
        prepared[mode]=attach_scores(rows,scores)
        thresholds[mode]=select_thresholds([r for r in prepared[mode] if r['split']=='calibration'])
    # All policies are frozen before the first test metric is calculated.
    save(out/'thresholds.json',thresholds);threshold_hash=F.T.sha256_file(out/'thresholds.json')
    results={mode:evaluate_rows([r for r in rows if r['split']=='test'],thresholds[mode]) for mode,rows in prepared.items()}
    save(out/'results.json',dict(scope=__doc__,results=results))
    with (out/'scored-supervision.jsonl').open('x') as f:
        for mode,rows in prepared.items():
            for row in rows:f.write(json.dumps(dict(mode=mode,**row),allow_nan=False)+'\n')
    F.check_hashes(reg['inputs'])
    if F.T.sha256_file(out/'thresholds.json')!=threshold_hash:raise ValueError('thresholds changed during evaluation')
    if F.T.sha256_file(out/'protocol.json')!=reg['protocol_sha256']:raise ValueError('comparison protocol changed during evaluation')
    save(out/'verification.json',dict(status='verified',inputs=reg['inputs'],outputs=F.hashes([p for p in out.iterdir() if p.is_file()])))
    return results


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['register','run'])
    for k in ('training','replay','out'):p.add_argument('--'+k,type=Path,required=k=='out')
    a=p.parse_args()
    if a.action=='register':
        if not a.training or not a.replay:p.error('training and replay required')
        register(a.training.resolve(),a.replay.resolve(),a.out.resolve())
    else:run(a.out.resolve())
