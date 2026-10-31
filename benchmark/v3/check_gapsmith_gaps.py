"""Audit before/after gap preparation against every row of a verified tool observation."""
import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import subprocess
import sys
from check_hmmer_route import digest


def check(observation, ec2go, before, after, out):
    receipt=json.loads((observation/'observation.json').read_text())
    if receipt['status']!='observed' or receipt['exit_code']:
        raise ValueError('external execution is not verified')
    for name,sha in receipt['output_hashes'].items():
        if digest(observation/name)!=sha: raise ValueError('external output changed')
    reactions=observation/'outputs/GCA_000012525.1-all-Reactions.tbl'
    pathways=observation/'outputs/GCA_000012525.1-all-Pathways.tbl'
    files=[reactions,pathways,ec2go,before,after,observation/'observation.json',Path(__file__).resolve(),
           Path(__file__).with_name('check_hmmer_route.py').resolve(),Path(sys.executable).resolve()]
    initial={str(p):digest(p) for p in files}
    out.mkdir(parents=True,exist_ok=False)
    runs=[]; records={}
    for arm,script in [('before',before),('after',after)]:
        output=out/(arm+'.jsonl')
        command=[sys.executable,str(script),'--pathways-tbl',str(pathways),'--reactions-tbl',str(reactions),
                 '--ec2go',str(ec2go),'--out',str(output),'--tag','GCA_000012525.1']
        with (out/(arm+'.log')).open('w') as log:
            result=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,timeout=30)
        runs.append(dict(arm=arm,command=command,exit_code=result.returncode))
        (out/'runs.json').write_text(json.dumps(runs,indent=2)+'\n')
        if result.returncode: raise ValueError('converter failed')
        records[arm]=[json.loads(line) for line in output.read_text().splitlines()]
    groups=defaultdict(list)
    for row in csv.DictReader(reactions.open(),delimiter='\t'):
        groups[(row['pathway'].strip('|'),row['rxn'])].append(row)
    active={r['ID'].strip('|') for r in csv.DictReader(pathways.open(),delimiter='\t') if r['Prediction']=='TRUE'}
    def identity(r): return (r['pathway_id'].removeprefix('MetaCyc:'),r['reaction_id'])
    def reasons(group):
        statuses={r['status'] for r in group}
        return [reason for reason,condition in [
            ('strong_hit_present','good_blast' in statuses),('reference_unknown','no_seq_data' in statuses),
            ('spontaneous',any(r['spont']=='TRUE' or r['status']=='spontaneous' for r in group)),
            ('complex_requires_subunit_inference',any(r['is_complex']=='TRUE' for r in group))] if condition]
    invalid=Counter()
    for record in records['before']:
        invalid.update(reasons(groups[identity(record)]))
    after_keys=set()
    for record in records['after']:
        key=identity(record);group=groups[key]
        if key[0] not in active or reasons(group): raise ValueError('ineligible candidate gap')
        if record['gapseq_guessed']!=any(r['status']=='bad_blast' for r in group): raise ValueError('order-dependent guess flag')
        after_keys.add(key)
    eligible={key for key,group in groups.items() if key[0] in active and not reasons(group)}
    if after_keys!=eligible: raise ValueError('eligible reaction missing or extra')
    if len({json.dumps(r,sort_keys=True) for r in records['after']})!=len(records['after']): raise ValueError('duplicate gap')
    if initial!={p:digest(Path(p)) for p in initial}:raise ValueError('inputs changed')
    report=dict(schema='gspa.gapsmith-gap-audit/1',status='verified',runs=runs,inputs=initial,
        before_records=len(records['before']),after_records=len(records['after']),eligible_reactions=len(eligible),
        before_ineligible_reason_counts=dict(invalid),reason_counts_overlap=True,
        scope='Search-supported candidate gap preparation on one training genome; no biological gap accuracy claim.',
        outputs={str(p.relative_to(out)):digest(p) for p in out.iterdir() if p.is_file()})
    (out/'verification.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:report[k] for k in ('status','before_records','after_records','eligible_reactions','before_ineligible_reason_counts')}))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('observation','ec2go','before','after','out'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();check(*(getattr(a,k).resolve() for k in ('observation','ec2go','before','after','out')))
