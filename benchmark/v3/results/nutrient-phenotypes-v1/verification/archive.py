from pathlib import Path
import gzip
import hashlib
import json
import shutil

repo = Path('/home/leechuck/Public/software/gspa')
data = Path.home() / '.cache/gspa-v3-data/nutrient-phenotypes-v1'
source_data = Path.home() / '.cache/gspa-v3-data/nutrient-environments-v1'
scratch = Path(__file__).parent
out = repo / 'benchmark/v3/results/nutrient-phenotypes-v1'


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for b in iter(lambda: stream.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


index = {}


def add(path, key):
    assert key not in index
    size = path.stat().st_size
    stored = key + '.gz' if size > 512 * 1024 else key
    target = out / stored
    target.parent.mkdir(parents=True, exist_ok=True)
    if stored != key:
        with path.open('rb') as src, target.open('wb') as dst, gzip.GzipFile(filename='', mode='wb', mtime=0, compresslevel=1, fileobj=dst) as z:
            shutil.copyfileobj(src, z, 1 << 20)
    else:
        shutil.copyfile(path, target)
    index[key] = dict(original=str(path), bytes=size, sha256=digest(path), stored=stored)


plan = json.loads((data / 'plan/protocol.json').read_text())
assert digest(data / 'plan/plans.jsonl') == plan['plans_sha256']
for name, expected in plan['implementation'].items():
    assert digest(repo / name) == expected
run = data / 'remote/jobs/10144'
assert (run / 'exit-code.txt').read_text().strip() == '1'
receipt = json.loads((run / 'predictions/run.json').read_text())
assert receipt['predictions_sha256'] == digest(run / 'predictions/predictions.jsonl')
assert receipt['protocol_sha256'] == digest(data / 'plan/protocol.json')
oracle = json.loads((run / 'oracle.json').read_text())
assert not oracle['agrees'] and oracle['comparisons'] == 13200 and oracle['disagreements'] == 3
assert oracle['implementation_sha256'] == digest(run / 'benchmark/v3/check_nutrient_cobra.py')
assert all(c['gpr_agrees'] for r in oracle['results'] for c in r['comparisons'])
final_root = data / 'remote/oracle-v2/10148'
assert (final_root / 'exit-code.txt').read_text().strip() == '1'
final_oracle = json.loads((final_root / 'oracle-v2.json').read_text())
assert final_oracle['implementation_sha256'] == digest(repo / 'benchmark/v3/check_nutrient_cobra.py')
assert final_oracle['implementation_sha256'] == digest(final_root / 'check_nutrient_cobra.py')
assert final_oracle['run_sha256'] == digest(run / 'predictions/run.json')
assert final_oracle['comparisons'] == 13200
assert not final_oracle['agrees'] and final_oracle['disagreements'] == 2
assert all(c['gpr_agrees'] for r in final_oracle['results'] for c in r['comparisons'])
unverified = [r for r in final_oracle['results'] if not r['agrees']]
assert [(r['case_id'], r['arm']) for r in unverified] == [('NUTRIENT2026:Bs:GNAc + NoC:M+Met', 'additional_inorganic:oxygen=20'), ('NUTRIENT2026:Bs:Glcn + Gly:R+R', 'printed:oxygen=0')]
assert sum(not c['flux_agrees'] for r in final_oracle['results'] for c in r['comparisons']) == 2
evaluation = json.loads((data / 'evaluation.json').read_text())
assert evaluation['run_sha256'] == digest(run / 'predictions/run.json')
assert evaluation['observations_sha256'] == digest(source_data / 'derived-v3/observations.jsonl')

out.mkdir(exist_ok=False)
for p in sorted(data.rglob('*')):
    if p.is_file():
        add(p, 'data/' + str(p.relative_to(data)))
for name in ['derived-v3/inputs.jsonl', 'derived-v3/observations.jsonl', 'derived-v3/report.json', 'models/registered-preflight.json']:
    add(source_data / name, 'observations/' + name)
for p in sorted(scratch.rglob('*')):
    if p.is_file() and 'package' not in p.relative_to(scratch).parts and p.name != 'archive.log':
        add(p, 'verification/' + str(p.relative_to(scratch)))
for name in ['benchmark/v3/run_nutrient_phenotypes.py', 'benchmark/v3/check_nutrient_cobra.py',
             'benchmark/v3/nutrient_phenotypes.sbatch', 'benchmark/v3/tests/test_nutrient_phenotypes.py',
             'scripts/metabolic_experiments.py', 'docs/benchmarks/v3-nutrient-phenotypes.md',
             'docs/benchmarks/v3-nutrient-environments.md', 'tasks/todo.md', 'tasks/v3-acceptance.md']:
    add(repo / name, 'source/' + name)
(out / 'storage-index.json').write_text(json.dumps(index, indent=2, sort_keys=True) + '\n')
for key, record in index.items():
    h = hashlib.sha256()
    opener = gzip.open if key != record['stored'] else open
    with opener(out / record['stored'], 'rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    assert h.hexdigest() == record['sha256'], key
report = dict(status='verified', payloads=len(index), original_bytes=sum(r['bytes'] for r in index.values()),
              stored_payload_bytes=sum((out / r['stored']).stat().st_size for r in index.values()))
(out / 'archive-verification.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report))
