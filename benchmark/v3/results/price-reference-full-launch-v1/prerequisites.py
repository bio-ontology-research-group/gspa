"""Check the completed, pinned full collection before read-only scoring."""
import hashlib
import json
from pathlib import Path
import sys


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def check():
    collection = Path('/data/hohndor/gspa-v3/price-production-v3/collection')
    finalizer = Path('/data/hohndor/gspa-v3/price-finalize-v2/jobs/10410')
    expected = {
        collection/'registration.json': '454f466cb635b70a9bbb65cbb70d1022065568d1931ac944be95e27a754b02dd',
        collection/'plan.json': '2407147c51dd0fe25cea0ea8193feefcac75d8cdbe46f0009d9001474339c4a4',
    }
    for path, digest in expected.items():
        if sha(path) != digest:
            raise ValueError(f'Pinned collection changed: {path}')
    plan = json.loads((collection/'plan.json').read_text())
    path = collection/'collection-verification.json'
    result = json.loads(path.read_text())
    if (finalizer/'exit-code.txt').read_text().strip() != '0':
        raise ValueError('Full collection finalizer did not succeed')
    if sha(path) != sha(finalizer/'collection-verification.json'):
        raise ValueError('Finalizer and collection verification differ')
    if (result['status'] != 'verified' or result['cells'] != 132 or result['arms'] != 2
            or result['expected_claims'] != plan['expected_claims']):
        raise ValueError('Incomplete collection verification')
    if result['expected_claims'] != 1526919081 or plan['expected_claims'] != 1526919081:
        raise ValueError('Claim denominator differs')
    if len(result['reports']) != 132 or [r['cell'] for r in result['reports']] != plan['cells']:
        raise ValueError('Verification does not cover the exact registered cells')
    for p, digest in expected.items():
        if result['artifacts'].get(str(p)) != digest:
            raise ValueError('Verification is not bound to the pinned collection')
    bound = dict(expected)
    for p in (path, finalizer/'collection-verification.json', finalizer/'exit-code.txt'):
        bound[p] = sha(p)
    return dict(status='verified', cells=132, expected_claims=1526919081,
                inputs={str(p): digest for p, digest in bound.items()},
                interpretation='Completed collection gate; evaluator rechecks full scoring inputs and artifacts.')


if __name__ == '__main__':
    with Path(sys.argv[1]).open('x') as stream:
        json.dump(check(), stream, indent=2, sort_keys=True)
        stream.write('\n')
