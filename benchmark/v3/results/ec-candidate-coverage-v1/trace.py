import csv,json,hashlib,collections
from pathlib import Path
cache=Path('/home/leechuck/.cache/gspa-v3-data'); out=Path('/home/leechuck/tmp/gspa-candidate-coverage-v1')
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
protocol=read(cache/'observed-gap-context-modes-v1/protocol.json'); bound=read(cache/'predictor-contract-v1/inputs.json')
ledger={r['path']:r['expected'] for r in map(json.loads,Path('benchmark/v3/results/context-comparison-audit-v1/job-10223/audit/file-checks.jsonl').read_text().splitlines())}
paths={'reference':cache/'training-cohort/baseline/reference_labels.tsv','hits':cache/'training-cohort/baseline/hits.tsv','proteins':cache/'training-cohort/prepared/proteins.tsv','labels':cache/'reviewed-thrs-locus-v1/labels.tsv','splits':Path('benchmark/v3/results/training-cohort/splits.tsv').resolve()}
inputs={}
for k,p in paths.items():
 h=sha(p);assert h==(bound if str(p) in bound else ledger)[str(p)];inputs[str(p)]=h
splits={r['genome_id']:r['split'] for r in csv.DictReader(paths['splits'].open(),delimiter='\t')}
labels=[r for r in csv.DictReader(paths['labels'].open(),delimiter='\t') if r['evidence']=='ECO:0000269' and r['label']=='1' and r['function_id'].startswith('EC:')]
refs=collections.defaultdict(set)
for r in csv.DictReader(paths['reference'].open(),delimiter='\t'):refs[r['function_id']].add(r['accession'])
seqs={(r['genome_id'],r['genome_protein_id']):r['protein_id'] for r in csv.DictReader(paths['proteins'].open(),delimiter='\t')}
hits=collections.defaultdict(list)
for line in paths['hits'].read_text().splitlines():
 q,s,p,c,e,b=line.split('\t');hits[q].append(dict(subject=s,identity=float(p),coverage=float(c),evalue=float(e),bitscore=float(b)))
claims={}; claimrecords=0
for mode,files in protocol['claims'].items():
 claims[mode]={}
 for gid,info in files.items():
  p=Path(info['path']);h=sha(p);assert h==info['sha256']==ledger[str(p)];inputs[str(p)]=h
  wanted={r['protein_id'] for r in labels if r['genome_id']==gid};values={}
  for line in p.open():
   r=json.loads(line);claimrecords+=1
   if r['function_type']=='EC' and r['protein_id'] in wanted:
    key=(r['protein_id'],'EC:'+r['function_id'].removeprefix('EC:'));assert key not in values;values[key]=r['raw_score']
  claims[mode][gid]=values
rows=[]
for label in labels:
 gid,pid,ec=label['genome_id'],label['protein_id'],label['function_id'];sid=seqs[gid,pid];support=[h for h in hits[sid] if h['subject'] in refs[ec]]
 for mode in claims:
  eligible=[h for h in support if mode=='homology' or not(h['identity']>=30 and h['coverage']>=80)]
  expected=max([round(h['identity']/100*min(h['coverage'],100)/100,6) for h in eligible]+[0.0]);present=(pid,ec) in claims[mode][gid]
  assert present==(expected>0)
  if present:assert claims[mode][gid][pid,ec]==expected
  stage=('emitted' if present else 'no_reference_annotation' if not refs[ec] else 'no_matching_search_hit' if not support else 'family_excluded' if not eligible else 'zero_score')
  rows.append(dict(genome_id=gid,protein_id=pid,sequence_id=sid,function_id=ec,split=splits[gid],mode=mode,reference_accessions=sorted(refs[ec]),query_hit_count=len(hits[sid]),matching_hits=support,eligible_hits=eligible,claim_score=expected if present else None,status=stage))
summary={m:{s:dict(collections.Counter(r['status'] for r in rows if r['mode']==m and r['split']==s)) for s in ['train','calibration','test']} for m in claims}
for p,h in inputs.items():assert sha(Path(p))==h
report=dict(status='verified',inputs=inputs,reference_pairs=len(labels),rows=len(rows),claims_scanned=claimrecords,summary=summary,interpretation='Exact-reference EC diagnostic only. Missing search hit does not distinguish sensitivity/evalue/max-target effects. Reference absence does not prove no enzyme or absent annotation in full Swiss-Prot. No target labels enter generation; no prediction changed.')
(out/'trace-rows.jsonl').write_text(''.join(json.dumps(r,sort_keys=True)+'\n' for r in rows));(out/'trace-report.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='inputs'},indent=2))
