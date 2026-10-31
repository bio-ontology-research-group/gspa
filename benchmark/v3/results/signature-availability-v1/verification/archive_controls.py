from pathlib import Path
import gzip,hashlib,json,shutil
repo=Path('/home/leechuck/Public/software/gspa');cache=Path.home()/'.cache/gspa-v3-data';scratch=Path(__file__).parent
out=repo/'benchmark/v3/results/signature-availability-v1';index=json.loads((out/'storage-index.json').read_text())
root=cache/'cohort-signature-availability-controls-v1'
assert (root/'all-genome-duplicates.json').is_file()
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def add(p,key):
 assert key not in index,key
 size=p.stat().st_size;stored=key+('.gz' if size>512*1024 else '');h=digest(p)
 dest=out/stored;dest.parent.mkdir(parents=True,exist_ok=True)
 if stored!=key:
  with p.open('rb') as src,dest.open('wb') as dst,gzip.GzipFile(filename='',mode='wb',compresslevel=1,mtime=0,fileobj=dst) as z:shutil.copyfileobj(src,z,1<<20)
 else:shutil.copyfile(p,dest)
 index[key]=dict(original=str(p),sha256=h,bytes=size,stored=stored)
for p in sorted(root.rglob('*')):
 if p.is_file():add(p,'controls/'+str(p.relative_to(root)))
for name in ['control-evaluation.log','controls-logistic.log','controls-gbdt.log','finish_controls.py','archive_controls.py']:
 add(scratch/name,'verification/'+name)
(out/'storage-index.json').write_text(json.dumps(index,indent=2,sort_keys=True)+'\n')
for key,r in index.items():
 h=hashlib.sha256();p=out/r['stored'];opener=gzip.open if key!=r['stored'] else open
 with opener(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 assert h.hexdigest()==r['sha256'],key
print(json.dumps({'status':'verified','payloads':len(index),'original_bytes':sum(r['bytes'] for r in index.values()),'stored_bytes':sum((out/r['stored']).stat().st_size for r in index.values())}),flush=True)
