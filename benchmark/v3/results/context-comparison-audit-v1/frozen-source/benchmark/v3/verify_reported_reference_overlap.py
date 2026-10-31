"""Verify full sequence identity for reported perfect-identity/coverage hits.

This is a lower-bound overlap audit of named hits, not a full reference-family
exclusion audit. It reads no functional truth labels.
"""
import argparse
from collections import defaultdict
import hashlib
import json
import os
from pathlib import Path


def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def sequences(path):
    result={};pid=None;sequence=[]
    def finish():
        if pid is None: return
        value=''.join(sequence).upper()
        if not value: raise ValueError('empty sequence: '+pid)
        if pid in result and result[pid]!=value: raise ValueError('conflicting sequence ID: '+pid)
        result[pid]=value
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            finish();pid=line[1:].split()[0];sequence=[]
        elif line.strip():
            if pid is None: raise ValueError('sequence before FASTA header')
            sequence.append(''.join(line.split()))
    finish();return result


def verify(plan):
    query={}
    hashes={}
    for name,expected in plan.get('source_inputs',{}).items():
        hashes[name]=digest(Path(name))
        if hashes[name]!=expected: raise ValueError('reported-hit source changed')
    for gid,entry in plan['proteomes'].items():
        path=Path(entry['path']);hashes[str(path)]=digest(path)
        if hashes[str(path)]!=entry['sha256']: raise ValueError('proteome file changed')
        query[gid]=sequences(path)
    byfile=defaultdict(list)
    for record in plan['pairs']: byfile[record['reference_file']].append(record)
    results=[]
    for filename,records in sorted(byfile.items()):
        path=Path(filename);hashes[filename]=digest(path)
        if hashes[filename]!=plan['reference_hashes'][filename]: raise ValueError('reference file changed')
        reference=sequences(path)
        for record in records:
            target=query[record['genome_id']][record['protein_id']]
            observed=reference[record['reference_id']]
            sequence_hash=lambda text:hashlib.sha256(text.encode('ascii')).hexdigest()
            if sequence_hash(target)!=record['query_sequence_sha256']:
                raise ValueError('query identity differs from frozen plan')
            results.append(dict(record,query_length=len(target),reference_length=len(observed),
                                reference_sequence_sha256=sequence_hash(observed),exact_sequence=(target==observed)))
    if hashes!={p:digest(Path(p)) for p in hashes}: raise ValueError('input changed while auditing')
    counts={}
    for gid in sorted(plan['proteomes']):
        rows=[r for r in results if r['genome_id']==gid]
        counts[gid]=dict(pairs=len(rows),proteins=len({r['protein_id'] for r in rows}),
                         exact_pairs=sum(r['exact_sequence'] for r in rows),
                         exact_proteins=len({r['protein_id'] for r in rows if r['exact_sequence']}))
    return dict(schema='gspa.reported-reference-overlap/1',status='verified',scope=__doc__,
                normalization='Uppercase sequence with FASTA whitespace removed; terminal stops retained.',
                counts=counts,results=results,inputs=hashes)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--plan-sha256',required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if not os.environ.get('SLURM_JOB_ID'): raise RuntimeError('run the reference scan through Slurm')
    if digest(a.plan)!=a.plan_sha256: raise ValueError('plan changed')
    result=verify(json.loads(a.plan.read_text()))
    if digest(a.plan)!=a.plan_sha256: raise ValueError('plan changed while scanning')
    result.update(plan_sha256=a.plan_sha256,source_sha256=digest(Path(__file__)),job_id=os.environ['SLURM_JOB_ID'])
    with a.out.open('x') as stream:json.dump(result,stream,indent=2,sort_keys=True);stream.write('\n')
    print(json.dumps({'status':'verified','counts':result['counts']},sort_keys=True))
