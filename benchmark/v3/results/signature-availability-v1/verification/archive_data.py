from pathlib import Path
import hashlib,gzip,json,shutil
repo=Path('/home/leechuck/Public/software/gspa');cache=Path.home()/'.cache/gspa-v3-data';scratch=Path(__file__).parent
out=repo/'benchmark/v3/results/signature-availability-v1';out.mkdir(exist_ok=False)
index={}
def digest(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
 return h.hexdigest()
def add(p,key):
 size=p.stat().st_size;h=digest(p);stored=key+('.gz' if size>512*1024 else '')
 dest=out/stored;dest.parent.mkdir(parents=True,exist_ok=True)
 if stored!=key:
  with p.open('rb') as src,dest.open('wb') as dst,gzip.GzipFile(filename='',mode='wb',compresslevel=1,mtime=0,fileobj=dst) as z:shutil.copyfileobj(src,z,1<<20)
 else:shutil.copyfile(p,dest)
 index[key]=dict(original=str(p),sha256=h,bytes=size,stored=stored)
for name,key in [('cohort-signature-before-v1','migration/before'),('cohort-signature-availability-v1','migration/after'),
                 ('signature-availability-cli-v1','cli/initial'),('signature-availability-cli-v2','cli/final')]:
 root=cache/name
 for p in sorted(root.rglob('*')):
  if p.is_file():add(p,key+'/'+str(p.relative_to(root)))
for kind in ['logistic','gbdt']:add(cache/f'cohort-signature-availability-v1/models/model_{kind}.json',f'models/model_{kind}.json')
for name in ['gspa-core/src/main/groovy/gspa/integration/cohort/CohortInference.groovy',
 'gspa-core/src/main/groovy/gspa/integration/contextbridge/CohortContextBuilder.groovy',
 'gspa-core/src/test/groovy/gspa/integration/cohort/CohortSignatureAvailabilitySpec.groovy',
 'gspa-core/src/test/groovy/gspa/integration/workflow/LocalReactionContextSpec.groovy',
 'gspa-core/src/test/groovy/gspa/integration/workflow/CohortWorkflowSpec.groovy',
 'benchmark/v3/check_signature_availability.py','docs/concepts/cohort-context.md']:
 add(repo/name,'source/'+name)
# Freeze every repository script bound by the model export, so original evidence
# does not depend on later edits to the live checkout.
identity=json.loads((cache/'cohort-signature-availability-v1/inputs.json').read_text())
for name,h in identity.items():
 p=Path(name)
 if p.suffix=='.py' and p.is_relative_to(repo):
  assert digest(p)==h
  add(p,'source/bound/'+str(p.relative_to(repo)))
for p in sorted(scratch.rglob('*')):
 if p.is_file() and not p.name.startswith('controls-') and p.name!='archive.log':add(p,'verification/'+str(p.relative_to(scratch)))
(out/'storage-index.json').write_text(json.dumps(index,indent=2,sort_keys=True)+'\n')
for key,r in index.items():
 h=hashlib.sha256();p=out/r['stored'];opener=gzip.open if key!=r['stored'] else open
 with opener(p,'rb') as f:
  for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
 assert h.hexdigest()==r['sha256'],key
print(json.dumps({'status':'verified','payloads':len(index),'original_bytes':sum(v['bytes'] for v in index.values()),'stored_bytes':sum((out/v['stored']).stat().st_size for v in index.values())}),flush=True)
