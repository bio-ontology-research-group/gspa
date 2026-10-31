"""Audit every candidate from twelve frozen observed-gap context-mode models.

At most two independent arms run concurrently; each uses the registered JVM
limits and four Python audit workers. All commands are registered before any
execution. Unsupported fallback, taxonomy and supplied gaps remain explicit.
"""
import argparse,csv,json
from concurrent.futures import ThreadPoolExecutor,ProcessPoolExecutor,as_completed
from multiprocessing import get_context
from pathlib import Path
import sys
import run_observed_gap_context_modes as R
import check_family_excluded_models as C
F=R.F


def make_plans(protocol,training,out,source_files):
    plans={};inputs={}
    for transfer in R.G.MODES:
        plans[transfer]={}
        for mode in R.MODES:
            baseline=protocol['plans'][transfer][mode]
            with Path(baseline['manifest']).open() as f:
                r=csv.DictReader(f,delimiter='\t');fields=r.fieldnames;original=list(r)
            if len(original)!=len(protocol['genomes']) or {r['genome_id'] for r in original}!=set(protocol['genomes']):raise ValueError('baseline genome coverage differs')
            plans[transfer][mode]={}
            for kind in F.KINDS:
                arm=out/transfer/mode/kind;arm.mkdir(parents=True,exist_ok=False)
                inference=arm/'inference';inference.mkdir();manifest=arm/'manifest.tsv'
                with manifest.open('w') as f:
                    w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
                    w.writerows(dict(r,out_path=str(inference/(r['genome_id']+'.integrated.tsv'))) for r in original)
                model=training/transfer/mode/'models'/f'model_{kind}.json';command=list(baseline['command'])
                for flag,value in [('--batch-manifest',manifest),('--feature-export-dir',arm/'features'),
                    ('--batch-stats-out',inference/'batch_stats.tsv'),('--features-out',inference/'rlc-features.tsv'),
                    ('--suggestions-out',inference/'rlc-suggestions.tsv'),('--learned-model',model)]:F.replace(command,flag,value)
                F.save(arm/'command.json',command)
                inputs.update(F.V.describe_run(command,manifest,sources=source_files)['inputs']);inputs.update(F.hashes([arm/'command.json']))
                plans[transfer][mode][kind]=dict(arm=str(arm),manifest=str(manifest),command=command,inference=str(inference),
                    features=str(arm/'features'),model=str(model),baseline=baseline,original_rows=original)
    return plans,inputs


def audit_arm(transfer,mode,kind,plan,protocol,source_files):
    F.V.run_batch(plan['command'],Path(plan['manifest']),Path(plan['inference']),sources=source_files)
    summaries={}
    with ProcessPoolExecutor(max_workers=4,mp_context=get_context('spawn')) as pool:
        jobs=[]
        for row in plan['original_rows']:
            gid=row['genome_id'];current=Path(plan['inference'])/(gid+'.integrated.tsv')
            if R.G.R.canonical(R.G.R.jsonl(Path(str(current)+'.gaps.jsonl')))!=R.G.R.canonical(R.G.R.jsonl(Path(protocol['gaps'][gid]['path']))):raise ValueError('learned context run changed supplied gaps')
            jobs.append((gid,pool.submit(C.audit_genome,Path(row['out_path']),current,
                Path(plan['baseline']['features'])/gid/'candidate_features.jsonl',
                Path(plan['features'])/gid/'candidate_features.jsonl',Path(plan['model']))))
        for gid,future in jobs:summaries[gid]=future.result()
    F.save(Path(plan['arm'])/'audit.json',summaries)
    print(json.dumps(dict(transfer=transfer,mode=mode,model=kind,genomes=len(summaries),candidates=sum(s['candidates'] for s in summaries.values()))),flush=True)
    return summaries


def run(training,out):
    protocol=R.verify(training);frozen=json.loads((training/'training.json').read_text());export=json.loads((training/'export.json').read_text())
    if frozen.get('status')!='verified' or export.get('status')!='verified':raise ValueError('all frozen models and exports must be verified')
    F.check_hashes(frozen['outputs']);F.check_hashes(export['outputs'])
    for transfer in R.G.MODES:
        for mode in R.MODES:
            for kind in F.KINDS:
                if str((training/transfer/mode/'models'/f'model_{kind}.json').resolve()) not in frozen['outputs']:raise ValueError('missing frozen model')
    out.mkdir(parents=True,exist_ok=False);source_files=R.sources()
    plans,inputs=make_plans(protocol,training,out,source_files)
    inputs.update(F.hashes([training/'protocol.json',training/'registration.json',training/'export.json',training/'training.json']))
    inputs.update(json.loads((training/'registration.json').read_text())['inputs']);inputs.update(export['outputs']);inputs.update(frozen['outputs'])
    F.save(out/'registration.json',dict(status='registered',registered_at=F.datetime.now(F.timezone.utc).isoformat(),
        scope=__doc__,plans=plans,inputs=inputs,max_concurrent_arms=2,audit_workers_per_arm=4))
    summaries={t:{m:{} for m in R.MODES} for t in R.G.MODES}
    with ThreadPoolExecutor(max_workers=2) as pool:
        jobs={pool.submit(audit_arm,t,m,k,p,protocol,source_files):(t,m,k) for t,modes in plans.items() for m,kinds in modes.items() for k,p in kinds.items()}
        try:
            for future in as_completed(jobs):
                t,m,k=jobs[future];summaries[t][m][k]=future.result()
        except BaseException:
            for future in jobs:future.cancel()
            raise
    F.check_hashes(inputs);R.verify(training)
    result=dict(status='verified',schema='gspa.observed-gap-context-model-verification/1',scope=__doc__,
        genomes=len(protocol['genomes']),summaries=summaries,inputs=inputs,outputs=F.hashes([p for p in out.rglob('*') if p.is_file()]))
    F.save(out/'verification.json',result)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('training','out'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();r=run(a.training.resolve(),a.out.resolve());print(json.dumps({k:r[k] for k in ('schema','status','genomes')}))
