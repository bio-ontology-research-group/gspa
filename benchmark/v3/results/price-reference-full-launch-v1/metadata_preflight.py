"""Read small completed-cell metadata only; do not score partial panels."""
import hashlib
import json
from pathlib import Path
import sys

root = Path(__file__).resolve().parent
sys.path.insert(0, str(root/'source/benchmark/v3'))
import evaluate_price_reference as evaluator

collection = Path('/data/hohndor/gspa-v3/price-production-v3/collection')
plan = evaluator.load(collection/'plan.json')
registration_sha = evaluator.sha(collection/'registration.json')
assert registration_sha == '454f466cb635b70a9bbb65cbb70d1022065568d1931ac944be95e27a754b02dd'
assert evaluator.sha(collection/'plan.json') == '2407147c51dd0fe25cea0ea8193feefcac75d8cdbe46f0009d9001474339c4a4'
assert evaluator.sha(root/'ec-audit.json') == '5302494b16ecd3f76fd96d24e48e27c298206b27608be72c5c66f1d197f4f79f'
refs, exclusions = evaluator.references(root/'source/benchmark/manifests/suite6_price149_manifest.tsv',
                                      {c['assembly'] for c in plan['cells']})
completed = []
for cell in plan['cells']:
    matches = list((collection/'attempts'/str(cell['index'])).glob('*/completed.json'))
    assert len(matches) <= 1
    if not matches:
        continue
    directory = matches[0].parent
    assert not (directory/'failed.json').exists()
    done = evaluator.load(matches[0])
    start = evaluator.load(directory/'started.json')
    assert done['cell'] == start['cell'] == cell and done['status'] == 'verified'
    assert start['registration_sha256'] == registration_sha
    path = directory/'production/verification.json'
    assert evaluator.sha(path) == done['verification_sha256']
    verification = evaluator.load(path)
    assert verification['status'] == 'verified'
    specs = {}
    for arm in ('stage0', 'stage1'):
        path = directory/'production'/arm/'verified_run.json'
        assert evaluator.sha(path) == verification['runs'][arm]['receipt_sha256']
        receipt = evaluator.load(path)
        assert receipt['returncode'] == 0
        specs[arm] = receipt['spec']
        assert all(Path(p).is_file() for p in receipt['artifacts'])
        assert all(Path(p).is_file() for p in receipt['spec']['inputs'])
    evaluator.comparable_specs(specs['stage0'], specs['stage1'], directory)
    completed.append(cell['index'])
print(json.dumps(dict(status='verified', completed_cell_metadata_checked=completed,
    expected_cells=len(plan['cells']), reference_proteins=len(refs),
    original_reference_pairs=sum(map(len, refs.values())),
    scope='Metadata compatibility only. No full-file hash verification or accuracy scoring of this incomplete collection.'), indent=2))
