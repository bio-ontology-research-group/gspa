from pathlib import Path
import gzip,hashlib,json,shutil,csv
repo=Path('/home/leechuck/Public/software/gspa');cache=Path.home()/'.cache/gspa-v3-data'
scratch=Path(__file__).parent;source=cache/'cohort-pathway-replay-v1'
out=repo/'benchmark/v3/results/cohort-pathway-replay-v1';out.mkdir(exist_ok=False)
def digest(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
 return h.hexdigest()
v=json.loads((source/'verification.json').read_text())
for p,h in v['inputs'].items():
 if digest(Path(p))!=h:raise ValueError('changed input '+p)
for p,h in v['outputs'].items():
 if digest(source/p)!=h:raise ValueError('changed output '+p)
index={}
def archive(p,key):
 h=digest(p);size=p.stat().st_size;stored=key+('.gz' if size>512*1024 else '')
 target=out/stored;target.parent.mkdir(parents=True,exist_ok=True)
 if stored!=key:
  with p.open('rb') as src,target.open('wb') as dst,gzip.GzipFile(filename='',mode='wb',compresslevel=1,mtime=0,fileobj=dst) as gz:shutil.copyfileobj(src,gz,1<<20)
 else:shutil.copyfile(p,target)
 index[key]=dict(original=str(p),bytes=size,sha256=h,stored=stored)
for p in sorted(source.rglob('*')):
 if p.is_file():archive(p,'cli/'+str(p.relative_to(source)))
for p in sorted((cache/'cohort-reported-overlap-v1').rglob('*')):
 if p.is_file():archive(p,'reference-overlap/'+str(p.relative_to(cache/'cohort-reported-overlap-v1')))
for name in ['replay_cohort_pathways.py','audit_cohort_pathways.py','observe_gapsmith.py','verify_reported_reference_overlap.py','tests/test_replay_cohort_pathways.py','tests/test_reported_reference_overlap.py']:
 archive(repo/'benchmark/v3'/name,'source/'+name)
for name in ['cli.log','python-final.log','overlap-submission.txt','overlap-first-job.txt','overlap-retry-submission.txt','overlap-retry-job.txt','compare_repeated_search.py','repeated-search-difference.json','archive_replay.py']:
 archive(scratch/name,'verification/'+name)
archive(repo/'benchmark/v3/results/training-cohort/splits.tsv','reference-overlap/cohort-splits.tsv')
# Source inspection is a hypothesis only; this checkout is not proven to build the installed binary.
for name in ['gapsmith-find/src/runner.rs','gapsmith-align/src/diamond.rs']:
 archive(repo.parent/'gapsmith/crates'/name,'upstream-inspection/'+name)
(out/'storage-index.json').write_text(json.dumps(index,indent=2,sort_keys=True)+'\n')
for key,record in index.items():
 p=out/record['stored'];h=hashlib.sha256()
 opener=gzip.open if record['stored']!=key else open
 with opener(p,'rb') as f:
  for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
 if h.hexdigest()!=record['sha256']:raise ValueError('archive payload mismatch '+key)
print(json.dumps({'status':'verified','payloads':len(index),'original_bytes':sum(r['bytes'] for r in index.values()),'stored_bytes':sum((out/r['stored']).stat().st_size for r in index.values())}),flush=True)
