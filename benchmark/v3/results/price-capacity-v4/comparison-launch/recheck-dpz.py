import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile

sys.path.insert(0, '/data/hohndor/gspa-v3/price-index-capacity-v1/code')
import verify_candidate_universe as U


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


part = Path('/data/hohndor/gspa-v3/price-full-matrix-v2/parts/dpz/GCF_000236665.1')
runs = {
    'before': Path('/data/hohndor/gspa-v3/price-prior-capacity-v1/dpz-GCF_000236665.1-stage1/stage1'),
    'after': Path('/data/hohndor/gspa-v3/price-index-capacity-v1/dpz-GCF_000236665.1-stage1/stage1')}
out = Path('/data/hohndor/gspa-v3/price-index-comparison-dpz-v1/verification.json')
if out.exists():
    raise SystemExit('comparison already exists')
publication = json.loads((part / 'job-receipt.json').read_text())
if publication['status'] != 'verified':
    raise ValueError('source part is unverified')
expected = {str(part / name): sha for name, sha in publication['artifacts'].items()}
for run in runs.values():
    receipt = json.loads((run / 'verified_run.json').read_text())
    if receipt['returncode'] != 0:
        raise ValueError('inference failed')
    for path, sha in {**receipt['spec']['inputs'], **receipt['artifacts']}.items():
        if path in expected and expected[path] != sha:
            raise ValueError('conflicting recorded hashes')
        expected[path] = sha
    for name in ['verified_run.json', 'resources.txt', 'command.json']:
        path = run / name
        expected[str(path)] = digest(path)
for path in [part / 'job-receipt.json', Path(__file__), Path(U.__file__), U.sort_executable()]:
    expected[str(path)] = digest(path)
if {p: digest(p) for p in expected} != expected:
    raise ValueError('recorded input/output changed')
counts = {}
with tempfile.TemporaryDirectory(dir=os.environ['TMPDIR']) as tmp:
    baseline, count = U.expected_candidates(part / 'claims.jsonl', Path(tmp))
    for arm, run in runs.items():
        counts[arm] = U.compare_candidates(baseline, run / 'GCF_000236665.1.integrated.tsv', Path(tmp))
if count != publication['claims']:
    raise ValueError('source universe changed')
predictions = {arm: digest(run / 'GCF_000236665.1.integrated.tsv') for arm, run in runs.items()}
if len(set(predictions.values())) != 1:
    raise ValueError('full production predictions changed')
if {p: digest(p) for p in expected} != expected:
    raise ValueError('input/output changed during comparison')
report = {'status': 'verified', 'input_candidates': count, 'arm_counts': counts,
          'prediction_sha256': predictions, 'predictions_byte_identical': True, 'hashes': expected,
          'method': 'Every candidate externally sorted and compared; full TSV bytes hashed; all recorded inputs/outputs rehashed twice.'}
out.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
print(json.dumps({k: v for k, v in report.items() if k != 'hashes'}))
