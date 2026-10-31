from pathlib import Path
import gzip, hashlib, json, shutil

repo=Path('/home/leechuck/Public/software/gspa')
data=Path.home()/'.cache/gspa-v3-data/nutrient-environments-v1'
scratch=Path(__file__).parent
out=repo/'benchmark/v3/results/nutrient-environments-v1'
out.mkdir(exist_ok=False)

def digest(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()

index={}
def add(path,key):
 h=digest(path);size=path.stat().st_size;stored=key+'.gz' if size>512*1024 else key
 target=out/stored;target.parent.mkdir(parents=True,exist_ok=True)
 if stored!=key:
  with path.open('rb') as src,target.open('wb') as dst,gzip.GzipFile(filename='',mode='wb',mtime=0,compresslevel=1,fileobj=dst) as z:shutil.copyfileobj(src,z,1<<20)
 else:shutil.copyfile(path,target)
 index[key]=dict(original=str(path),bytes=size,sha256=h,stored=stored)

final=json.loads((data/'derived-v3/report.json').read_text())
assert digest(repo/'benchmark/v3/prepare_nutrient_environments.py')==final['implementation_sha256']
for name,h in final['source_hashes'].items():assert digest(data/'raw'/name)==h,name
for name,h in final['outputs'].items():assert digest(data/'derived-v3'/name)==h,name
preflight=json.loads((data/'models/registered-preflight.json').read_text())
assert digest(data/'derived-v3/inputs.jsonl')==preflight['inputs_sha256']
assert digest(repo/'benchmark/v3/register_nutrient_models.py')==preflight['implementation_sha256']
for r in preflight['models'].values():assert digest(Path(r['path']))==r['sha256']
for p in sorted(data.rglob('*')):
 if p.is_file():add(p,'data/'+str(p.relative_to(data)))
add(Path(preflight['models']['Ec']['path']),'data/models/iML1515.xml')
for p in sorted(scratch.rglob('*')):
 if p.is_file() and p.name!='archive.log':add(p,'verification/'+str(p.relative_to(scratch)))
for name in ['benchmark/v3/prepare_nutrient_environments.py','benchmark/v3/register_nutrient_models.py',
             'benchmark/v3/tests/test_nutrient_environments.py','benchmark/v3/prepare_metabolic_benchmarks.py',
             'scripts/metabolic_experiments.py','docs/benchmarks/v3-nutrient-environments.md',
             'tasks/todo.md','tasks/v3-acceptance.md']:
 add(repo/name,'source/'+name)
(out/'storage-index.json').write_text(json.dumps(index,indent=2,sort_keys=True)+'\n')
for key,r in index.items():
 opener=gzip.open if key!=r['stored'] else open;h=hashlib.sha256()
 with opener(out/r['stored'],'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 assert h.hexdigest()==r['sha256'],key
report=dict(status='verified',payloads=len(index),original_bytes=sum(r['bytes'] for r in index.values()),stored_payload_bytes=sum((out/r['stored']).stat().st_size for r in index.values()))
(out/'archive-verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report))
