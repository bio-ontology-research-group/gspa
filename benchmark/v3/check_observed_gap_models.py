"""Run four frozen gap-enabled models through full joint-cohort production.

Every candidate is compared with the original Python reader. Unsupported and
taxon-forbidden outputs must preserve the matched gap-enabled rule baseline.
No truth labels enter the score comparisons; source files are hashed for integrity.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context
import csv
import json
from pathlib import Path
import sys

import fit_observed_gap_models as B
import check_family_excluded_models as C
F=B.F
G=B.G


def prepare_plans(protocol, training, out, sources):
    """Register every mode/model command before the first production invocation."""
    plans={};inputs={}
    for mode in G.MODES:
        baseline=protocol['plans'][mode]['observed_gaps']
        with Path(baseline['manifest']).open() as f:
            reader=csv.DictReader(f,delimiter='\t');fields=reader.fieldnames;original=list(reader)
        if len(original)!=len(protocol['genomes']) or {r['genome_id'] for r in original}!=set(protocol['genomes']):
            raise ValueError('baseline genome coverage differs')
        plans[mode]={}
        for kind in F.KINDS:
            arm=out/mode/kind;arm.mkdir(parents=True,exist_ok=False)
            inference=arm/'inference';inference.mkdir()
            manifest=arm/'manifest.tsv'
            rows=[dict(r,out_path=str(inference/(r['genome_id']+'.integrated.tsv'))) for r in original]
            with manifest.open('w') as f:
                w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
            model=training/mode/'default/models'/f'model_{kind}.json'
            command=list(baseline['command'])
            for flag,value in [('--batch-manifest',manifest),('--feature-export-dir',arm/'features'),
                ('--batch-stats-out',inference/'batch_stats.tsv'),('--learned-model',model),
                ('--features-out',inference/'rlc-features.tsv'),('--suggestions-out',inference/'rlc-suggestions.tsv')]:
                F.replace(command,flag,value)
            F.save(arm/'command.json',command)
            spec=F.V.describe_run(command,manifest,sources=sources)
            inputs.update(spec['inputs']);inputs.update(F.hashes([arm/'command.json']))
            plans[mode][kind]=dict(arm=str(arm),inference=str(inference),features=str(arm/'features'),
                manifest=str(manifest),model=str(model),command=command,baseline=baseline,original_rows=original)
    return plans,inputs


def run(training,out):
    fitting=B.verify(training)
    frozen=json.loads((training/'training.json').read_text())
    if frozen.get('status')!='verified':raise ValueError('verified frozen models required')
    F.check_hashes(frozen['outputs'])
    comparison=Path(fitting['comparison']);protocol=G.verify(comparison)
    export=json.loads((comparison/'verification.json').read_text())
    if export.get('status')!='verified':raise ValueError('complete verified gap export required')
    F.check_hashes(export['outputs'])
    out.mkdir(parents=True,exist_ok=False)
    sources=[Path(m.__file__).resolve() for m in tuple(sys.modules.values())
             if getattr(m,'__file__',None) and Path(m.__file__).suffix=='.py'
             and Path(m.__file__).resolve().is_relative_to(F.ROOT)]
    plans,inputs=prepare_plans(protocol,training,out,sources)
    inputs.update(F.hashes([training/'training.json',training/'registration.json',training/'protocol.json',
                           comparison/'verification.json',comparison/'registration.json',comparison/'protocol.json']))
    inputs.update(frozen['outputs']);inputs.update(export['outputs'])
    inputs.update(json.loads((training/'registration.json').read_text())['inputs'])
    F.save(out/'registration.json',dict(status='registered',registered_at=F.datetime.now(F.timezone.utc).isoformat(),
        scope=__doc__,plans=plans,inputs=inputs))
    summaries={}
    for mode in G.MODES:
        summaries[mode]={}
        for kind in F.KINDS:
            plan=plans[mode][kind]
            F.V.run_batch(plan['command'],Path(plan['manifest']),Path(plan['inference']),sources=sources)
            summaries[mode][kind]={}
            with ProcessPoolExecutor(max_workers=4,mp_context=get_context('spawn')) as pool:
                jobs=[]
                for row in plan['original_rows']:
                    gid=row['genome_id'];current=Path(plan['inference'])/(gid+'.integrated.tsv')
                    if G.R.canonical(G.R.jsonl(Path(str(current)+'.gaps.jsonl')))!=G.R.canonical(G.R.jsonl(Path(protocol['gaps'][gid]['path']))):
                        raise ValueError('learned run changed supplied gap records')
                    jobs.append((gid,pool.submit(C.audit_genome,Path(row['out_path']),current,
                        Path(plan['baseline']['features'])/gid/'candidate_features.jsonl',
                        Path(plan['features'])/gid/'candidate_features.jsonl',Path(plan['model']))))
                for gid,future in jobs:
                    result=future.result();summaries[mode][kind][gid]=result
                    print(json.dumps(dict(mode=mode,model=kind,genome=gid,audit=result)),flush=True)
            F.save(Path(plan['arm'])/'audit.json',summaries[mode][kind])
    F.check_hashes(inputs);B.verify(training);G.verify(comparison)
    result=dict(schema='gspa.observed-gap-model-verification/1',status='verified',scope=__doc__,
        genomes=len(protocol['genomes']),summaries=summaries,inputs=inputs,
        outputs=F.hashes([p for p in out.rglob('*') if p.is_file()]))
    F.save(out/'verification.json',result)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--training',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();result=run(a.training.resolve(),a.out.resolve())
    print(json.dumps({k:result[k] for k in ('schema','status','genomes')}))
