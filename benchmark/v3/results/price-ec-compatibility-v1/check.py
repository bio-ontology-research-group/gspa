"""Whole small Price cell: latest JAR versus retained pre-correction outputs."""
import argparse, hashlib, json, sys
from pathlib import Path
parser=argparse.ArgumentParser()
parser.add_argument('--out',type=Path,required=True)
a=parser.parse_args()
base=Path('/data/hohndor/gspa-v3/price-production-v3')
sys.path.insert(0,str(base/'code'))
import run_price_matrix_capacity as C
jar=Path('/data/hohndor/gspa-v3/ec-prefix-identity-v3/jobs/10402/gspa-3.0.0.jar')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(jar)=='37aeee2b565951dd3a4e049a23aacd2e8155f29bd8b6216fed7006e4273965e6'
plan=json.loads((base/'collection/plan.json').read_text())
reg=json.loads((base/'collection/registration.json').read_text())
assert sha(base/'collection/plan.json')==reg['plan_sha256']
cell=plan['cells'][0]
assert cell['claims']==3481 and cell['assembly']=='GCF_000236665.1' and cell['predictor']=='topec'
old=base/'collection/attempts/0/job-10292/production'
inputs={str(jar):sha(jar),str(Path(__file__)):sha(Path(__file__)),str(base/'collection/plan.json'):sha(base/'collection/plan.json')}
for p in (base/'code').glob('*.py'):
 assert reg['inputs'][str(p)]==sha(p)
 inputs[str(p)]=sha(p)
for arm in ('stage0','stage1'):
 receipt=json.loads((old/arm/'verified_run.json').read_text())
 assert C.V.output_hashes(receipt['spec'])==receipt['artifacts']
 inputs.update({str(p):sha(p) for p in (old/arm).iterdir() if p.is_file() and p.name!='.run.lock'})
C.run(Path(cell['part']),Path(plan['layout']),jar,Path(plan['reference']),a.out/'production',scratch=a.out.parent,java_heap='6g',timeout=1800)

def canonical_ledger(value):
 if isinstance(value,list): return [canonical_ledger(v) for v in value]
 if not isinstance(value,dict): return value
 out={}
 for k,v in value.items():
  if k=='elapsed_seconds': continue
  if k=='function_key':
   parts=v.rsplit('|',2)
   if len(parts)==3 and parts[1]=='EC': parts[2]=parts[2].removeprefix('EC:')
   out[k]='|'.join(parts)
  elif k=='function_id' and value.get('function_type')=='EC': out[k]=v.removeprefix('EC:')
  else: out[k]=canonical_ledger(v)
 return out

results={}
for arm in ('stage0','stage1'):
 name=cell['assembly']+'.integrated.tsv'
 before,after=old/arm/name,a.out/'production'/arm/name
 tsv_equal=before.read_bytes()==after.read_bytes()
 before_ledger=json.loads(Path(str(before)+'.stages.json').read_text())
 after_ledger=json.loads(Path(str(after)+'.stages.json').read_text())
 semantic_equal=canonical_ledger(before_ledger)==canonical_ledger(after_ledger)
 results[arm]=dict(prediction_bytes_equal=tsv_equal,ledger_equal_after_explicit_prefix_and_timing_normalization=semantic_equal,candidates=3481)
 assert tsv_equal and semantic_equal,results
assert all(sha(Path(p))==digest for p,digest in inputs.items())
outputs={str(p):sha(p) for p in a.out.rglob('*') if p.is_file()}
(a.out/'verification.json').write_text(json.dumps(dict(status='verified',scope='One complete TopEC cell, both legacy prior settings; not all Price cells or biological accuracy.',results=results,inputs=inputs,outputs=outputs,results_root=str(a.out)),indent=2)+'\n')
print(json.dumps(results),flush=True)
