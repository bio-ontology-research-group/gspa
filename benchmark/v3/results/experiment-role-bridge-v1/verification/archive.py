from pathlib import Path
import gzip, hashlib, json, shutil, xml.etree.ElementTree as ET

repo = Path('/home/leechuck/Public/software/gspa')
scratch = Path(__file__).parent
out = repo / 'benchmark/v3/results/experiment-role-bridge-v1'

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def totals(folder):
    result = dict(tests=0, failures=0, errors=0, skipped=0)
    for p in folder.glob('TEST-*.xml'):
        r = ET.parse(p).getroot()
        for k in result: result[k] += int(r.get(k, 0))
    return result

before = totals(scratch / 'before-xml')
focused = totals(scratch / 'after-xml')
full = totals(scratch / 'full-jvm-xml')
assert before['failures'] == 12
assert focused == dict(tests=105, failures=0, errors=0, skipped=0)
assert full == dict(tests=1278, failures=0, errors=0, skipped=6)
replay = json.loads((scratch / 'cli-replay-final/verification.json').read_text())
assert replay['status'] == 'verified' and len(replay['cases']) == 9
assert sum(o['oracle_comparisons'] for c in replay['cases'] for o in c['observations'].values()) == 90
out.mkdir(exist_ok=False)
index = {}

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
sources = ['gspa-core/src/main/groovy/gspa/integration/contextbridge/' + n + '.groovy'
    for n in ['MaskedAssignmentBuilder', 'SupportedAssignmentHypotheses', 'AbsenceComparisonBridge',
              'GuidedAbsenceBuilder', 'ExperimentStage', 'GuidedAbsenceStage', 'SbmlModelInfo', 'SidecarRun']]
sources += ['gspa-core/src/main/groovy/gspa/integration/learned/ContextFeatureExtractor.groovy',
    'gspa-core/src/main/groovy/gspa/integration/learned/CandidateRole.groovy']
sources += ['gspa-core/src/test/groovy/gspa/integration/workflow/' + n + '.groovy'
    for n in ['MaskedAssignmentBuilderSpec', 'GuidedAbsenceBuilderSpec', 'ExperimentWorkflowSpec']]
sources += ['gspa-core/src/test/groovy/gspa/integration/learned/ClaimRoleBridgeSpec.groovy',
    'gspa-core/src/test/resources/experiments/make_fixtures.py',
    'scripts/metabolic_experiments.py', 'scripts/masked_assignment_experiments.py', 'scripts/absence_experiments.py',
    'benchmark/v3/check_absence_cobra.py', 'benchmark/v3/check_masked_assignment_cobra.py',
    'benchmark/v3/run_experiment_role_bridge.py', 'docs/benchmarks/v3-experiment-role-bridge.md',
    'docs/benchmarks/v3-masked-assignment-experiments.md', 'docs/benchmarks/v3-guided-absence.md',
    'tasks/todo.md', 'tasks/v3-acceptance.md']
for name in sources: add(repo / name, 'source/' + name)
for key, r in index.items():
    p = out / r['stored']
    data = gzip.decompress(p.read_bytes()) if r['stored'] != key else p.read_bytes()
    assert len(data) == r['bytes'] and hashlib.sha256(data).hexdigest() == r['sha256'], key
    assert sha(Path(r['source'])) == r['sha256'], key
(out / 'storage-index.json').write_text(json.dumps(index, indent=2, sort_keys=True) + '\n')
summary = dict(status='verified', payloads=len(index), original_bytes=sum(r['bytes'] for r in index.values()),
    stored_bytes=sum((out / r['stored']).stat().st_size for r in index.values()), failing_before_cases=12,
    focused_tests=focused, full_jvm_tests=full, final_cli_runs=18, independent_phenotype_comparisons=90,
    whole_item_accepted=False)
(out / 'archive-verification.json').write_text(json.dumps(summary, indent=2) + '\n')
(out / 'README.md').write_text('''# Experiment role handoff: verification evidence

Raw candidate roles now survive masked assignment construction. Complex-subunit
observations require supported partner clauses. Noncatalytic-only candidates
cannot become enzyme assignments, and malformed role metadata fails validation.
Explicit and guided absence comparisons also reject a standalone activity from
mixed catalytic/subunit observations. No scientific accuracy claim is made.

Twelve failing-before regression cases are retained in `verification/before-xml/`.
All 105 focused cases pass in `after-xml/`. The full JVM suite passes 1,278 tests
with six optional skips, recorded in `full-jvm-xml/` and `full-jvm.log`.

Final CLI evidence is `verification/cli-replay-final/`: nine matched old/new
cases, 18 invocations and 90 independent GLPK comparisons. Five invalid masked
assignments now abstain; one mixed-role absence request fails validation without
publishing a result. Eight successfully completed pairs preserve annotation
bytes. Valid default/supported cases retain identical phenotype cells; masked
public problems are unchanged. Proposed hypotheses retain role provenance.

`verification/cli-replay/` is an incomplete initial harness run: its unknown
assembly fixture omitted the required nullable completeness_fraction field.
The initial harness and failure are retained; use the final replay for results.

Sources, exact commands, package versions, inputs, logs, test XML and old/new CLI
outputs are included. Engine hashes and immutable local cache paths are in the
receipts; JARs are not duplicated. Replays preserve original workspace paths;
remap paths on another host. The source report documents the scope and runner.
Large files use deterministic gzip. The storage index records original sizes
and hashes; every payload was checked after expansion against its original file.
Run `sha256sum -c SHA256SUMS` to verify stored files.

All role, partner and model fixtures are synthetic. Biological role accuracy,
independent partner evidence and calibrated recovery remain open. No whole
audit item is accepted.
''')
files = sorted(p for p in out.rglob('*') if p.is_file() and p.name != 'SHA256SUMS')
(out / 'SHA256SUMS').write_text(''.join(sha(p) + '  ' + str(p.relative_to(out)) + '\n' for p in files))
summary['checksummed_files'] = len(files)
print(json.dumps(summary))
