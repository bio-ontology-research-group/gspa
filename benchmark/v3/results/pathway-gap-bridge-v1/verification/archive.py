from pathlib import Path
import gzip,hashlib,json
root=Path.cwd();cache=Path.home()/'.cache/gspa-v3-data';scratch=Path.home()/'tmp/gspa-pathway-gap-bridge-v1'
out=root/'benchmark/v3/results/pathway-gap-bridge-v1';out.mkdir(exist_ok=False)
index={}
def save(src,relative):
 dst=out/relative;dst.parent.mkdir(parents=True,exist_ok=True)
 data=src.read_bytes();sha=hashlib.sha256(data).hexdigest()
 if len(data)>524288:
  dst=Path(str(dst)+'.gz')
  with dst.open('wb') as raw:
   with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0,compresslevel=6) as gz:gz.write(data)
 else:dst.write_bytes(data)
 index[str(relative)]=dict(stored=str(dst.relative_to(out)),sha256=sha,bytes=len(data),original=str(src))
for p in scratch.rglob('*'):
 if p.is_file() and p.name!='archive.log':save(p,Path('verification')/p.relative_to(scratch))
for name,label in [('pathway-gap-before-v1','before-cli'),('pathway-gap-after-v1','after-cli')]:
 for p in (cache/name).rglob('*'):
  if p.is_file():save(p,Path(label)/p.relative_to(cache/name))
 print('Archived',label,flush=True)
r=json.loads((cache/'pathway-gap-after-v1/verification.json').read_text())
for i,(p,h) in enumerate(r['inputs'].items()):
 src=Path(p)
 if src.suffix=='.jar' or src.name=='java' or src.is_relative_to(cache/'pathway-gap-after-v1'):continue
 assert hashlib.sha256(src.read_bytes()).hexdigest()==h
 save(src,Path('inputs')/f'{i:03d}-{src.name}')
source=['gspa-core/src/main/groovy/gspa/integration/'+n for n in ['contextbridge/PathwayGapBuilder.groovy','IntegrationState.groovy','MetabolicGap.groovy','workflow/IntegrationWorkflow.groovy','workflow/IntegrationWorkflowConfig.groovy','workflow/IntegrationReferenceLoader.groovy','learned/ContextFeatureExtractor.groovy']]
source+=['gspa-core/src/main/resources/workflow/integration-workflow-config.keys.json','gspa-cli/src/main/groovy/gspa/cli/IntegrateCommand.groovy']
source+=['gspa-predictors/src/main/groovy/gspa/predictor/'+n for n in ['AnnotationPipeline.groovy','pathway/GapseqPredictor.groovy']]
source+=['gspa-core/src/test/groovy/gspa/integration/workflow/'+n for n in ['PathwayGapBuilderSpec.groovy','PathwayGapWorkflowSpec.groovy']]
source+=['gspa-predictors/src/test/groovy/gspa/predictor/pathway/PathwayObservationSpec.groovy']
source+=['benchmark/v3/'+n for n in ['check_pathway_gap_route.py','check_hmmer_route.py','v3_runner.py','tests/test_verified_batch.py','verified_batch.py','tests/test_nextflow_cohort.py']]
source+=['docs/benchmarks/v3-pathway-gap-bridge.md','docs/concepts/context-workflow.md']
source_hashes={}
for p in source:
 save(root/p,Path('source')/p);source_hashes[p]=hashlib.sha256((root/p).read_bytes()).hexdigest()
(out/'source-files.json').write_text(json.dumps(source_hashes,indent=2,sort_keys=True)+'\n')
for module in ['gspa-core','gspa-cli','gspa-predictors']:
 for p in (root/module/'build/test-results/test').glob('TEST-*.xml'):save(p,Path('test-results')/module/p.name)
(out/'storage-index.json').write_text(json.dumps(index,indent=2,sort_keys=True)+'\n')
print('Archived',len(index),'files',flush=True)
