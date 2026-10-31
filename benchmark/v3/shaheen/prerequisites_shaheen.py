"""Check the completed, pinned release-candidate Price collection on Shaheen III before read-only scoring.

Mirrors benchmark/v3/results/price-reference-full-launch-v1/prerequisites.py with parameters
for collection root, finalizer job directory, register receipt, and expected plan/registration hashes.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def find_finalizer_dir(collection_dir, expected_verification_sha=None):
    """Locate finalizer job directory under collection's parent jobs directory."""
    jobs_dir = collection_dir.parent / 'jobs'
    if not jobs_dir.is_dir():
        return None
    candidates = []
    for job_path in sorted(jobs_dir.iterdir(), key=lambda p: p.stat().st_mtime, reverse=True):
        if not job_path.is_dir():
            continue
        vfile = job_path / 'collection-verification.json'
        efile = job_path / 'exit-code.txt'
        if vfile.is_file() and efile.is_file():
            if expected_verification_sha is None or sha(vfile) == expected_verification_sha:
                candidates.append(job_path)
    return candidates[0] if candidates else None


def check(collection_dir, finalizer_dir=None, receipt_file=None,
          expected_registration=None, expected_plan=None, allow_incomplete=False):
    collection = Path(collection_dir).resolve()
    receipt_path = Path(receipt_file).resolve() if receipt_file else collection.parent / 'register-receipt.json'

    if receipt_path.is_file():
        receipt = json.loads(receipt_path.read_text())
        if expected_registration is None:
            expected_registration = receipt.get('registration_sha256')
        if expected_plan is None:
            expected_plan = receipt.get('plan_sha256')

    # Known fallback for Shaheen release-candidate collection if receipt is absent
    if expected_registration is None:
        expected_registration = '90a960f8a5c1b7c89307a1cb67a27c3c668320ba333c4a87a02f3b874f9337c9'
    if expected_plan is None:
        expected_plan = '6ea6903399ec633c722fc485d8eda5391efd1765ee190f678a8faa2ef34a67e7'

    expected = {
        collection / 'registration.json': expected_registration,
        collection / 'plan.json': expected_plan,
    }
    for path, digest in expected.items():
        if not path.is_file():
            raise FileNotFoundError(f'Missing collection file: {path}')
        actual = sha(path)
        if actual != digest:
            raise ValueError(f'Pinned collection file changed: {path} (expected {digest}, got {actual})')

    plan = json.loads((collection / 'plan.json').read_text())
    path = collection / 'collection-verification.json'

    if not path.is_file():
        if allow_incomplete:
            completed_cells = []
            for cell in plan['cells']:
                if list((collection / 'attempts' / str(cell['index'])).glob('*/completed.json')):
                    completed_cells.append(cell['index'])
            return dict(
                status='partial',
                completed_cells=len(completed_cells),
                expected_cells=len(plan['cells']),
                completed_cell_indices=completed_cells,
                interpretation='Partial collection check; collection-verification.json not yet created.'
            )
        raise ValueError(f'Full collection verification missing (collection not finalized): {path}')

    result = json.loads(path.read_text())
    verification_sha = sha(path)

    if finalizer_dir is None:
        finalizer = find_finalizer_dir(collection, verification_sha)
        if finalizer is None:
            raise ValueError(f'Could not find successful finalizer job directory for {path}')
    else:
        finalizer = Path(finalizer_dir).resolve()

    if not (finalizer / 'exit-code.txt').is_file():
        raise ValueError(f'Finalizer exit-code.txt missing: {finalizer}')
    if (finalizer / 'exit-code.txt').read_text().strip() != '0':
        raise ValueError('Full collection finalizer did not succeed (exit-code != 0)')
    if not (finalizer / 'collection-verification.json').is_file():
        raise ValueError(f'Finalizer collection-verification.json missing: {finalizer}')
    if verification_sha != sha(finalizer / 'collection-verification.json'):
        raise ValueError('Finalizer and collection verification digests differ')

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
    for p in (path, finalizer / 'collection-verification.json', finalizer / 'exit-code.txt'):
        bound[p] = sha(p)

    return dict(
        status='verified',
        cells=132,
        expected_claims=1526919081,
        inputs={str(p): digest for p, digest in bound.items()},
        interpretation='Completed Shaheen RC collection gate; evaluator rechecks full scoring inputs and artifacts.'
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('positional_out', nargs='?', type=Path, default=None,
                        help='Path to write prerequisites JSON output')
    parser.add_argument('--collection', type=Path,
                        default=Path('/scratch/project/k10018/gspa-v3/out/price-production-rc1/collection'),
                        help='Path to collection directory')
    parser.add_argument('--finalizer', type=Path, default=None,
                        help='Path to finalizer job directory')
    parser.add_argument('--receipt', type=Path, default=None,
                        help='Path to register-receipt.json')
    parser.add_argument('--expected-registration', default=None,
                        help='Expected sha256 of collection/registration.json')
    parser.add_argument('--expected-plan', default=None,
                        help='Expected sha256 of collection/plan.json')
    parser.add_argument('--allow-incomplete', action='store_true',
                        help='Allow unfinalized collection (returns partial status instead of raising)')
    parser.add_argument('--out', type=Path, default=None,
                        help='Path to write prerequisites JSON output')
    args = parser.parse_args()

    out_path = args.out or args.positional_out
    if not out_path:
        parser.error('An output path is required (pass positional arg or --out)')

    res = check(
        collection_dir=args.collection,
        finalizer_dir=args.finalizer,
        receipt_file=args.receipt,
        expected_registration=args.expected_registration,
        expected_plan=args.expected_plan,
        allow_incomplete=args.allow_incomplete
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open('w') as stream:
        json.dump(res, stream, indent=2, sort_keys=True)
        stream.write('\n')
    print(f"Prerequisites written to {out_path} (status: {res.get('status')})")


if __name__ == '__main__':
    main()
