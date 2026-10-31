from pathlib import Path
import gzip
import hashlib
import json

root = Path.cwd()
cache = Path.home()/'.cache/gspa-v3-data'
scratch = Path.home()/'tmp/gspa-rlc-mapping-v1'
out = root/'benchmark/v3/results/rlc-mapping-v1'
out.mkdir(exist_ok=False)
index = {}

def save(src, relative):
    dst = out/relative
    dst.parent.mkdir(parents=True, exist_ok=True)
    data = src.read_bytes()
    if len(data) > 524288:
        dst = Path(str(dst)+'.gz')
        with dst.open('wb') as raw:
            with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0, compresslevel=6) as gz:
                gz.write(data)
    else:
        dst.write_bytes(data)
    index[str(relative)] = dict(stored=str(dst.relative_to(out)), bytes=len(data),
                               sha256=hashlib.sha256(data).hexdigest(), original=str(src))

for p in scratch.rglob('*'):
    if p.is_file() and p.name != 'archive.log':
        save(p, Path('verification')/p.relative_to(scratch))
names = ['rlc-mapping-v1','rlc-mapping-missing-evidence-v1','rlc-mapping-gap-bridge-v1']
inputs = {}
for name in names:
    directory = cache/name
    for p in directory.rglob('*'):
        if p.is_file(): save(p, Path('cli')/name/p.relative_to(directory))
    receipt = json.loads((directory/'verification.json').read_text())
    for p,h in receipt['inputs'].items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest() == h, p
        inputs[p] = h
    for p,h in receipt['outputs'].items():
        assert hashlib.sha256((directory/p).read_bytes()).hexdigest() == h, p
for i,(p,h) in enumerate(sorted(inputs.items())):
    src = Path(p)
    if src.suffix == '.jar' or src.name == 'java' or any(src.is_relative_to(cache/n) for n in names): continue
    save(src, Path('inputs')/f'{i:03d}-{src.name}')

source = ['gspa-core/src/main/groovy/gspa/'+n for n in [
    'ontology/PathwayLoader.groovy','integration/IntegrationState.groovy',
    'integration/contextbridge/PathwayGapBuilder.groovy',
    'integration/workflow/IntegrationReferenceLoader.groovy',
    'integration/workflow/IntegrationWorkflow.groovy',
    'integration/suggester/ReactionLocalContextSuggester.groovy']]
source += ['gspa-core/src/test/groovy/gspa/'+n for n in [
    'ontology/EcToGoAlternativesSpec.groovy','integration/workflow/RlcMappingWorkflowSpec.groovy',
    'integration/suggester/ReactionLocalContextSuggesterSpec.groovy',
    'integration/workflow/PathwayGapBuilderSpec.groovy','integration/workflow/PathwayGapWorkflowSpec.groovy',
    'integration/workflow/ContextFixtures.groovy','integration/workflow/ReactionContextEligibilitySpec.groovy']]
source += ['benchmark/v3/'+n for n in ['check_rlc_mapping.py','check_rlc_missing_evidence.py',
    'check_pathway_gap_route.py','check_hmmer_route.py','tests/test_nextflow_cohort.py']]
source += ['docs/benchmarks/v3-rlc-mapping.md','docs/concepts/context-workflow.md',
           'tasks/todo.md','tasks/v3-acceptance.md']
source_hashes = {}
for p in source:
    save(root/p, Path('source')/p)
    source_hashes[p] = hashlib.sha256((root/p).read_bytes()).hexdigest()
for module in ['gspa-core','gspa-cli','gspa-predictors']:
    for p in (root/module/'build/test-results/test').glob('TEST-*.xml'):
        save(p, Path('test-results')/module/p.name)
(out/'source-files.json').write_text(json.dumps(source_hashes, indent=2, sort_keys=True)+'\n')
(out/'storage-index.json').write_text(json.dumps(index, indent=2, sort_keys=True)+'\n')
print('Archived', len(index), 'uncompressed payloads', flush=True)
