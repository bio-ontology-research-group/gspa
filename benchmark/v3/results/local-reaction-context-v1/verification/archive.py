from pathlib import Path
import gzip,hashlib,json,shutil
root=Path.cwd(); cache=Path.home()/'.cache/gspa-v3-data'; scratch=Path.home()/'tmp/gspa-local-context-v1'
out=root/'benchmark/v3/results/local-reaction-context-v1';out.mkdir(exist_ok=False)
index={}
def save(src,relative):
 dst=out/relative;dst.parent.mkdir(parents=True,exist_ok=True)
 data=src.read_bytes(); sha=hashlib.sha256(data).hexdigest()
 if len(data)>524288:
  dst=Path(str(dst)+'.gz')
  with dst.open('wb') as raw:
   with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0,compresslevel=6) as gz: gz.write(data)
 else: dst.write_bytes(data)
 index[str(relative)]=dict(stored=str(dst.relative_to(out)),sha256=sha,bytes=len(data),original=str(src))
for src in scratch.rglob('*'):
 if src.is_file(): save(src,Path('verification')/src.relative_to(scratch))
for name,label in [('local-reaction-context-v1','failed-cli'),('local-reaction-context-v2','cli')]:
 for src in (cache/name).rglob('*'):
  if src.is_file(): save(src,Path(label)/src.relative_to(cache/name))
 print('Archived',label,flush=True)
receipt=json.loads((cache/'local-reaction-context-v2/verification.json').read_text())
for i,(p,h) in enumerate(receipt['inputs'].items()):
 src=Path(p)
 if src.suffix=='.jar' or src.name=='java' or src.is_relative_to(cache/'local-reaction-context-v2'): continue
 assert hashlib.sha256(src.read_bytes()).hexdigest()==h
 save(src,Path('inputs')/f'{i:03d}-{src.name}')
source=['gspa-core/src/main/groovy/gspa/integration/'+n for n in ['contextbridge/CohortContextBuilder.groovy','contextbridge/CohortStage.groovy','learned/CohortFeatureColumns.groovy','learned/ExternalContext.groovy','learned/ContextFeatureExtractor.groovy','learned/FeatureExporter.groovy','workflow/IntegrationWorkflow.groovy']]
source+=['gspa-core/src/test/groovy/gspa/integration/workflow/'+n for n in ['CohortWorkflowSpec.groovy','LocalReactionContextSpec.groovy','ContextFixtures.groovy']]
source+=['benchmark/v3/'+n for n in ['check_local_reaction_context.py','check_hmmer_route.py','prepare_training_data.py','tests/test_training_review_fixes.py','tests/test_nextflow_cohort.py']]
source+=['docs/benchmarks/v3-local-reaction-context.md']
source_hashes={}
for p in source:
 save(root/p,Path('source')/p);source_hashes[p]=hashlib.sha256((root/p).read_bytes()).hexdigest()
(out/'source-files.json').write_text(json.dumps(source_hashes,indent=2,sort_keys=True)+'\n')
for module in ['gspa-core','gspa-cli','gspa-predictors']:
 for src in (root/module/'build/test-results/test').glob('TEST-*.xml'): save(src,Path('test-results')/module/src.name)
(out/'storage-index.json').write_text(json.dumps(index,indent=2,sort_keys=True)+'\n')
print('Archived',len(index),'files',flush=True)
