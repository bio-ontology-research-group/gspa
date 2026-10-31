from pathlib import Path
import gzip, hashlib, json, shutil

repo = Path('/home/leechuck/Public/software/gspa')
scratch = Path(__file__).parent
cache = Path.home() / '.cache/gspa-v3-data/nutrient-community-v1'
out = repo / 'benchmark/v3/results/nutrient-community-checkpoint-v1'
out.mkdir(exist_ok=False)
index = {}

def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda: f.read(1048576), b''): h.update(chunk)
    return h.hexdigest()

def add(p, key):
    stored = key + '.gz' if p.stat().st_size > 524288 else key
    dest = out / stored
    dest.parent.mkdir(parents=True, exist_ok=True)
    if stored == key: shutil.copy2(p, dest)
    else:
        with p.open('rb') as src, dest.open('wb') as dst, gzip.GzipFile(filename='', mode='wb', mtime=0, fileobj=dst) as z:
            shutil.copyfileobj(src, z)
    index[key] = dict(source=str(p), stored=stored, bytes=p.stat().st_size, sha256=sha(p))

assert 'Ran 7 tests' in (scratch / 'focused-final.log').read_text()
assert 'Ran 645 tests' in (scratch / 'benchmark-gate.log').read_text()
assert (scratch / 'benchmark-gate.log').read_text().rstrip().endswith('OK (skipped=13)')
smoke = json.loads((scratch / 'real-model-smoke.json').read_text())
assert len(smoke) == 12 and sum(len(c['checks']) for c in smoke) == 60
assert all(r['agrees'] for c in smoke for r in c['checks'])
for p in sorted(scratch.rglob('*')):
    if p.is_file() and p.name != 'archive.log' and '__pycache__' not in p.parts:
        add(p, 'verification/' + str(p.relative_to(scratch)))
package = cache / 'package-v1'
for p in sorted(package.rglob('*')):
    if p.is_file(): add(p, 'package/' + str(p.relative_to(package)))
for name in ['benchmark/v3/tests/test_nutrient_community.py', 'docs/benchmarks/v3-nutrient-community.md',
             'docs/benchmarks/v3-nutrient-phenotypes.md', 'docs/benchmarks/v3-nutrient-environments.md',
             'tasks/todo.md', 'tasks/v3-acceptance.md']:
    add(repo / name, 'source/' + name)
raw = cache.parent / 'nutrient-environments-v1'
for name in ['derived-v3/inputs.jsonl', 'models/registered-preflight.json']:
    add(raw / name, 'inputs/' + name)
for key, r in index.items():
    p = out / r['stored']
    data = gzip.decompress(p.read_bytes()) if r['stored'] != key else p.read_bytes()
    assert len(data) == r['bytes'] and hashlib.sha256(data).hexdigest() == r['sha256'], key
    assert sha(Path(r['source'])) == r['sha256'], key
(out / 'storage-index.json').write_text(json.dumps(index, indent=2, sort_keys=True) + '\n')
summary = dict(status='checkpoint_verified', full_evaluation_status='running_not_verified',
    payloads=len(index), original_bytes=sum(r['bytes'] for r in index.values()),
    stored_bytes=sum((out / r['stored']).stat().st_size for r in index.values()),
    focused_tests=7, benchmark_tests=645, optional_skips=13, preliminary_real_model_comparisons=60,
    cases=1320, primary_cases=1200, plans=7920, planned_production_calls=39600,
    slurm_array=10150, whole_item_accepted=False)
(out / 'archive-verification.json').write_text(json.dumps(summary, indent=2) + '\n')
(out / 'README.md').write_text('''# Nutrient coculture execution checkpoint

This is a verified implementation/input checkpoint, not a completed empirical
evaluation. Slurm array 10150 is running the full production and independent
GLPK replay. Submission, checksum and running-state receipts are retained.

The outcome-blind plan contains all 1,320 cocultures (1,200 primary), six
medium/oxygen arms and five LPs per case/arm: 7,920 plans, 39,600 production calls.
Separate genotypes and organism models preserve full GPR effects and native
maintenance. The three community policies cover permissive uptake, native
unsupplied uptake and closure of unsupplied uptake. Bulk endpoint OD does not
prove both members grew or identify exchanged metabolites/transporters.

Seven focused tests pass, including independent synthetic reciprocal-rescue
checks and shard/serial identity. The 645-test benchmark gate passes with 13
optional skips. Sixty deterministic preliminary real-model comparisons agree
with GLPK; no claim is made for the full catalogue yet. Initial resource-warning
logs are retained; the final focused gate is warning-free.

`package/` is the exact immutable remote package: source, models and final
plan-v3 bindings. No observation file is in that package. Development plan
versions have identical plan rows but older implementation hashes; their hash
comparison is recorded separately. Use only the final frozen package to run.

Current periodic production messages report original interleaved plan positions,
not per-shard completion counts; some shards emit no periodic messages. The
checkpoint read actual files through a read-only Slurm step after worker SSH
was unavailable: 243 and 234 written records, still awaiting full verification.

Sources, commands, package versions, preliminary results and assumptions are
retained. Paths in receipts preserve the original environment and need remapping
elsewhere. The full source report documents the protocol and limitations. Large
files use deterministic gzip; `storage-index.json` records uncompressed sizes
and hashes. Every payload was expanded and verified against its original file.
Run `sha256sum -c SHA256SUMS` to check stored files.

No full numerical or biological accuracy result, transporter validation,
gene-attribution result or whole-item acceptance is claimed.
''')
files = sorted(p for p in out.rglob('*') if p.is_file() and p.name != 'SHA256SUMS')
(out / 'SHA256SUMS').write_text(''.join(sha(p) + '  ' + str(p.relative_to(out)) + '\n' for p in files))
summary['checksummed_files'] = len(files)
print(json.dumps(summary))
