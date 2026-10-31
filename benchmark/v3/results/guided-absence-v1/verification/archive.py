from pathlib import Path
import gzip, hashlib, json, shutil

repo = Path('/home/leechuck/Public/software/gspa')
scratch = Path(__file__).parent
out = repo / 'benchmark/v3/results/guided-absence-v1'
out.mkdir(exist_ok=True)
index = {}

def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''): h.update(block)
    return h.hexdigest()

def add(p, key):
    stored = key + '.gz' if p.stat().st_size > 524288 else key
    dest = out / stored
    dest.parent.mkdir(parents=True, exist_ok=True)
    if stored == key: shutil.copy2(p, dest)
    else:
        with p.open('rb') as src, dest.open('wb') as dst, gzip.GzipFile(filename='', mode='wb', fileobj=dst, mtime=0) as z:
            shutil.copyfileobj(src, z)
    index[key] = dict(source=str(p), stored=stored, bytes=p.stat().st_size, sha256=sha(p))

for p in sorted(scratch.rglob('*')):
    if p.is_file() and '__pycache__' not in p.parts and p.name != 'archive.log':
        add(p, 'verification/' + str(p.relative_to(scratch)))
sources = [
    'gspa-core/src/main/groovy/gspa/integration/contextbridge/' + n + '.groovy'
    for n in ['GuidedAbsenceBuilder', 'GuidedAbsenceStage', 'ExperimentStage', 'AbsenceComparisonBridge', 'SbmlModelInfo',
              'SidecarRun', 'PathwayGapBuilder']]
sources += ['gspa-core/src/main/groovy/gspa/integration/workflow/IntegrationWorkflow.groovy',
    'gspa-core/src/main/groovy/gspa/integration/workflow/IntegrationWorkflowConfig.groovy',
    'gspa-core/src/main/groovy/gspa/integration/learned/ContextFeatureExtractor.groovy',
    'gspa-core/src/test/groovy/gspa/integration/workflow/ExperimentWorkflowSpec.groovy',
    'gspa-core/src/test/groovy/gspa/integration/workflow/GuidedAbsenceBuilderSpec.groovy',
    'gspa-core/src/test/resources/experiments/make_fixtures.py',
    'scripts/metabolic_experiments.py', 'scripts/absence_experiments.py',
    'benchmark/v3/check_absence_cobra.py', 'benchmark/v3/run_guided_absence_mechanism.py',
    'docs/benchmarks/v3-guided-absence.md', 'docs/benchmarks/v3-absence-comparison.md',
    'docs/scenarios/metabolic-experiments.md', 'tasks/todo.md', 'tasks/v3-acceptance.md']
for name in sources: add(repo / name, 'source/' + name)
for key, r in index.items():
    p = out / r['stored']
    data = gzip.decompress(p.read_bytes()) if r['stored'] != key else p.read_bytes()
    assert len(data) == r['bytes'] and hashlib.sha256(data).hexdigest() == r['sha256'], key
    assert sha(Path(r['source'])) == r['sha256'], key
(out / 'storage-index.json').write_text(json.dumps(index, indent=2, sort_keys=True) + '\n')
report = dict(status='verified', payloads=len(index), original_bytes=sum(r['bytes'] for r in index.values()),
    stored_bytes=sum((out / r['stored']).stat().st_size for r in index.values()), focused_tests=50,
    jvm_tests=1265, jvm_optional_skips=6, python_mechanism_tests=89, cli_runs=6,
    independent_phenotype_comparisons=100, whole_item_accepted=False)
(out / 'archive-verification.json').write_text(json.dumps(report, indent=2) + '\n')
(out / 'README.md').write_text('''# Guided absence mechanism: verified scope

Automatic construction connects run gaps and retained raw candidates to separate
absence comparisons with shared budgets. These are synthetic mechanism checks,
not biological absence classification or gene-function validation.

Final CLI evidence is in `verification/cli-replay/`: six packaged CLI runs,
including a baseline, and 100 independent COBRApy/GLPK phenotype comparisons.
All agree, while annotation and learned-feature exports remain byte-identical.
The full JVM gate has 1,265 tests, no failures and six optional skips; 50 focused
tests and 89 Python mechanism tests pass. Exact commands and package versions
are retained. Four optional JVM skips require the installed-HMMER environment
variable; two are existing optional core tests.

`focused-initial` retains the mapping callback's set/list failure.
`focused-second` retains a test assertion using nonexistent annotation fields.
`focused-third.log` and `focused-final-xml/` are the passing focused results.
`full-jvm-xml/` and `full-jvm.log` are the passing full-suite results.

Source snapshots, inputs, outputs, child logs and construction provenance are
retained. The immutable engine location and SHA-256 are in
`verification/engine.json`; the JAR is cached locally, not duplicated here.
The replay runner and receipts retain original workspace paths: remap them when
restoring on another host. The source report documents the input contract and
reproduction command. Large files use deterministic gzip; `storage-index.json`
records their uncompressed hashes and sizes. Every stored payload was expanded
and checked against both its recorded hash and the original source.

Run `sha256sum -c SHA256SUMS` to check stored files. No whole audit item is accepted.
''')
files = sorted(p for p in out.rglob('*') if p.is_file() and p.name != 'SHA256SUMS')
(out / 'SHA256SUMS').write_text(''.join(sha(p) + '  ' + str(p.relative_to(out)) + '\n' for p in files))
report['checksummed_files'] = len(files)
print(json.dumps(report))
