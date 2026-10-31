"""Independent audit of published real-probe bytes and all sequence evidence."""
from collections import defaultdict
import hashlib,json,os,re
from pathlib import Path
ROOT=Path('/data/hohndor/gspa-v3/multisource-comparison-probe-v4/jobs/10403')
OUT=Path('/data/hohndor/gspa-v3/source-probe-audit-v2/jobs')/os.environ['SLURM_JOB_ID']

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
 return h.hexdigest()

assert (ROOT/'exit-code.txt').read_text().strip()=='0'
path=ROOT/'results/verification.json';probe_sha=sha(path);r=json.loads(path.read_text())
assert r['status']=='verified' and r['schema']=='gspa.source-context-probe/2'
assert r['ec_identity']==dict(combined_candidates=29512,raw_claim_files_unchanged=True,shared_source_candidates=511,verified_both_source_candidates=511)
assert len(r['outputs'])==145
outputs={ROOT/'results'/Path(p).relative_to(r['results_root']):h for p,h in r['outputs'].items()}
for p,h in outputs.items():assert sha(p)==h,str(p)
print('All 145 published output hashes match.',flush=True)
gid='GCA_000007005.1'

def key(row):
 kind=row['function_type'];fid=row['function_id']
 if kind=='EC':
  fid=fid.removeprefix('EC:')
  assert re.fullmatch(r'[1-7]\.(?:\d+|-)\.(?:\d+|-)\.(?:\d+|n\d+|-)',fid)
  fid='EC:'+fid
 return row['protein_id'],kind,fid

profiles={};raw={}
for arm,source in [('diamond_homology','swissprot_diamond_baseline'),('clean','clean')]:
 p=ROOT/'results/claims'/arm/(gid+'.claims.jsonl')
 raw[arm]=p.read_bytes();rows={}
 for line in raw[arm].splitlines():
  row=json.loads(line);k=key(row)
  assert row['source']==source and k not in rows
  assert isinstance(row['metadata']['predictor_contract'],dict)
  rows[k]={source:dict(score=row['raw_score'],contract=row['metadata']['predictor_contract'])}
 profiles[arm]=rows
combined=ROOT/'results/claims/clean_diamond_homology'/(gid+'.claims.jsonl')
assert combined.read_bytes() in (raw['diamond_homology']+raw['clean'],raw['clean']+raw['diamond_homology'])
merged=defaultdict(dict)
for rows in profiles.values():
 for k,sources in rows.items():merged[k].update(sources)
profiles['clean_diamond_homology']=dict(merged)
shared=set(profiles['clean'])&set(profiles['diamond_homology'])
assert len(shared)==511 and all(k[1]=='EC' for k in shared)
assert {a:len(p) for a,p in profiles.items()}=={'diamond_homology':26328,'clean':3695,'clean_diamond_homology':29512}
checks={}
for arm,expected in profiles.items():
 for mode in ['SEQONLY','UNIVERSAL','WITHOUT_RLC']:
  p=(ROOT/'results/without-rlc'/arm/'UNIVERSAL' if mode=='WITHOUT_RLC' else ROOT/'results/runs'/arm/mode)/'features'/gid/'candidate_features.jsonl'
  assert p in outputs
  seen=set();zero=0;both=0;digest=hashlib.sha256()
  with p.open('rb') as f:
   for line in f:
    digest.update(line);row=json.loads(line);k=key(row)
    assert k not in seen and k in expected
    assert row['function_id']==k[2] and row['genome_id']==gid
    if k[1]=='EC':assert row['role']=='CATALYTIC'
    seen.add(k);sources=expected[k];v=row['features']
    assert row['sequence_sources']==sorted(sources)
    contracts=[dict(source=s,contract=sources[s]['contract']) for s in sorted(sources)]
    assert row['sequence_contracts']==contracts
    assert v['seq.n_sources_reporting']==v['seq.n_groups_reporting']==len(sources)
    assert v['seq.max_raw_any']==max(x['score'] for x in sources.values())
    for source,column,group in [('clean','clean','ml_protein_seq'),('swissprot_diamond_baseline','other','homology')]:
     score=sources.get(source,{}).get('score',0)
     for prefix in ['seq.src.'+column,'seq.group.'+group]:
      assert v[prefix+'.max_raw']==score
      assert v[prefix+'.reported']==int(source in sources)
    if mode=='SEQONLY':
     for name,value in v.items():
      if name.startswith(('ctx.','cohort.','asm.','fba.')):assert value==0;zero+=1
    both+=int(len(sources)==2)
  assert digest.hexdigest()==outputs[p]
  assert seen==set(expected)
  assert both==(511 if arm=='clean_diamond_homology' else 0)
  checks[arm+'/'+mode]=dict(candidates=len(seen),both_source_candidates=both,zero_context_values=zero)
  print(arm,mode,checks[arm+'/'+mode],flush=True)
assert sha(path)==probe_sha
for arm,data in raw.items():assert hashlib.sha256(data).hexdigest()==outputs[ROOT/'results/claims'/arm/(gid+'.claims.jsonl')]
assert sha(combined)==outputs[combined]
report=dict(status='verified',probe_receipt_sha256=probe_sha,published_outputs_checked=len(outputs),checks=checks,
 input_sequence_files_unchanged=True,method='Independent stdlib reader; no GSPA ingestion or benchmark audit modules imported.',
 scope='All nine one-genome feature exports and published bytes. Neither function truth, learned accuracy, nor recurrence across genomes.')
with (OUT/'verification.json').open('x') as f:json.dump(report,f,indent=2);f.write('\n')
