"""Verify every retained prediction/ledger; external inputs are receipt-bound only."""
import csv,hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parent
OLD=Path('/data/hohndor/gspa-v3/price-production-v3/collection/attempts/0/job-10292/production')
NEW=Path('/scratch/10405/price-ec/results/production')
JARS={'before':('/data/hohndor/gspa-v3/price-production-v3/gspa-3.0.0.jar','0606bec0bf009bc6b9c1301a28ae8d21998730c36138c24cb5720a44eb31a50a'),
      'after':('/data/hohndor/gspa-v3/ec-prefix-identity-v3/jobs/10402/gspa-3.0.0.jar','37aeee2b565951dd3a4e049a23aacd2e8155f29bd8b6216fed7006e4273965e6')}
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()

def normalize(value):
 if isinstance(value,list):return [normalize(v) for v in value]
 if not isinstance(value,dict):return value
 result={}
 for k,v in value.items():
  if k=='elapsed_seconds':continue
  if k=='function_key':
   parts=v.rsplit('|',2)
   if len(parts)==3 and parts[1]=='EC':parts[2]=parts[2].removeprefix('EC:')
   result[k]='|'.join(parts)
  elif k=='function_id' and value.get('function_type')=='EC':result[k]=v.removeprefix('EC:')
  else:result[k]=normalize(v)
 return result

def verify():
 summary={}
 for arm in ['stage0','stage1']:
  receipts={};rows={};ledgers={}
  name='GCF_000236665.1.integrated.tsv'
  for side,remote in [('before',OLD),('after',NEW)]:
   folder=ROOT/side/arm
   receipt=json.loads((folder/'verified_run.json').read_text());receipts[side]=receipt
   assert receipt['returncode']==0
   assert receipt['spec']['inputs'][JARS[side][0]]==JARS[side][1]
   expected={str(remote/arm/name),str(remote/arm/(name+'.stages.json'))}
   assert set(receipt['artifacts'])==expected
   for path,digest in receipt['artifacts'].items():
    assert sha(ROOT/side/Path(path).relative_to(remote))==digest
   manifest=str(remote/arm/'manifest.tsv')
   assert sha(folder/'manifest.tsv')==receipt['spec']['inputs'][manifest]
   command=receipt['spec']['command']
   assert ('--enable-priors' in command)==(arm=='stage0')
   if arm=='stage0':assert command[command.index('--enable-priors')+1]=='none'
   with (folder/name).open() as f:
    rows[side]=list(csv.DictReader(f,delimiter='\t'))
   assert len(rows[side])==3481
   assert len({(r['protein_id'],r['type'],r['function_id']) for r in rows[side]})==3481
   ledgers[side]=json.loads((folder/(name+'.stages.json')).read_text())
  a,b=(receipts[s]['spec']['inputs'] for s in ['before','after'])
  assert all(a[k]==b[k] for k in a.keys()&b.keys())
  assert a.keys()-b.keys()=={JARS['before'][0],str(OLD/arm/'manifest.tsv')}
  assert b.keys()-a.keys()=={JARS['after'][0],str(NEW/arm/'manifest.tsv')}
  for x,y in zip(rows['before'],rows['after']):
   assert x['type']==y['type']=='EC'
   assert not x['function_id'].startswith('EC:')
   assert y['function_id']=='EC:'+x['function_id']
   assert {k:v for k,v in x.items() if k!='function_id'}=={k:v for k,v in y.items() if k!='function_id'}
  assert normalize(ledgers['before'])==normalize(ledgers['after'])
  summary[arm]=dict(candidates=3481,identifier_prefix_changes=3481,other_tsv_field_changes=0,
    semantic_ledger_equal=True,shared_receipt_input_hashes_equal=len(a.keys()&b.keys()))
 return dict(status='verified',arms=summary,scope='One complete TopEC cell; both prior settings. No all-cell or accuracy conclusion.')

if __name__=='__main__':
 manifest=ROOT/'archive-hashes.json'
 if manifest.exists():
  for name,h in json.loads(manifest.read_text()).items():assert sha(ROOT/name)==h,name
 print(json.dumps(verify(),indent=2,sort_keys=True))
