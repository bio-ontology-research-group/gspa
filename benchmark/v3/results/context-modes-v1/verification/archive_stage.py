from pathlib import Path
import gzip,hashlib,json,shutil,sys
repo=Path('/home/leechuck/Public/software/gspa');cache=Path.home()/'.cache/gspa-v3-data';scratch=Path(__file__).parent
out=repo/'benchmark/v3/results/context-modes-v1';out.mkdir(exist_ok=True)
index_path=out/'storage-index.json'
index=json.loads(index_path.read_text()) if index_path.exists() else {}
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def add(p,key):
 size=p.stat().st_size;sha=digest(p);stored=key+('.gz' if size>512*1024 else '')
 if key in index:
  assert index[key]['sha256']==sha,key
  return
 dest=out/stored;dest.parent.mkdir(parents=True,exist_ok=True)
 if stored!=key:
  with p.open('rb') as src,dest.open('wb') as dst,gzip.GzipFile(filename='',mode='wb',compresslevel=1,mtime=0,fileobj=dst) as z:shutil.copyfileobj(src,z,1<<20)
 else:shutil.copyfile(p,dest)
 index[key]=dict(original=str(p),sha256=sha,bytes=size,stored=stored)
def tree(root,key):
 for p in sorted(root.rglob('*')):
  if p.is_file() and '__pycache__' not in p.parts and not (p.name.startswith('archive-') and p.suffix=='.log'):add(p,key+'/'+str(p.relative_to(root)))
stage=sys.argv[1]
if stage=='acquisition':
 tree(cache/'cohort-pathway-search-archaea-v2','acquisition/archaea-searches')
 tree(cache/'archaea-reference-v1','acquisition/reference-receipts')
 add(scratch/'complete-acquisition-audit.json','acquisition/full-cohort-audit.json')
 add(scratch/'archaeal-acquisition-completed.txt','acquisition/reference-job-completed.txt')
 add(scratch/'search-10131-state.txt','acquisition/search-array-completed.txt')
elif stage=='exports':
 for name,key in [('cohort-context-modes-v1-inputs','inputs/initial'),('cohort-context-modes-v2-inputs','inputs/final'),
                  ('cohort-context-modes-v1','exports/initial'),('cohort-context-modes-v2','exports/final')]:
  tree(cache/name,key)
elif stage=='verification':
 tree(scratch,'verification')
 for name in json.loads((scratch/'source-files.json').read_text()):add(repo/name,'source/'+name)
else:raise ValueError(stage)
index_path.write_text(json.dumps(index,indent=2,sort_keys=True)+'\n')
for key,r in index.items():
 h=hashlib.sha256();p=out/r['stored'];opener=gzip.open if key!=r['stored'] else open
 with opener(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 assert h.hexdigest()==r['sha256'],key
print(json.dumps(dict(status='verified',stage=stage,payloads=len(index),original_bytes=sum(v['bytes'] for v in index.values()),stored_bytes=sum((out/v['stored']).stat().st_size for v in index.values()))),flush=True)
