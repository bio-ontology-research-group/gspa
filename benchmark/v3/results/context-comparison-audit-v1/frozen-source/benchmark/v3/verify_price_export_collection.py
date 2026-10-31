#!/usr/bin/env python3
"""Verify complete plan coverage and rehash every published full-matrix part.

The per-part verifier compares all original matrix cells with exported claims.
This independent collection gate binds those receipts to the frozen plan and
checks the persistent files, including every byte of every claims file. It does
not rerun the matrix comparison or measure biological accuracy.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def verify(plan_path, progress=None):
    plan_path = Path(plan_path).resolve()
    plan = json.loads(plan_path.read_text())
    root = Path(plan['out']) / 'parts'
    jobs = plan['jobs']
    keys = [(j['predictor'], j['assembly']) for j in jobs]
    if not jobs or len(set(keys)) != len(keys):
        raise ValueError('empty plan or duplicate planned matrix')
    if set(keys) != {(p, a) for p in {p for p, _ in keys} for a in {a for _, a in keys}}:
        raise ValueError('plan does not cover each predictor/assembly combination')
    actual = {(p.parent.name, p.name) for p in root.glob('*/*') if p.is_dir()}
    if actual != set(keys):
        raise ValueError('published part coverage differs from plan')
    expected = {}

    def bind(path, sha):
        path = str(Path(path).resolve())
        if path in expected and expected[path] != sha:
            raise ValueError(f'conflicting input identity: {path}')
        expected[path] = sha

    bind(plan_path, digest(plan_path))
    bind(__file__, digest(__file__))
    for path, sha in plan['inputs'].items():
        bind(path, sha)
    parts = []
    required = {'claims.jsonl', 'join_map.tsv', 'part.json', 'matrix-verification.json'}
    for index, job in enumerate(jobs):
        part = root / job['predictor'] / job['assembly']
        receipt_path = part / 'job-receipt.json'
        receipt = json.loads(receipt_path.read_text())
        if (receipt.get('status') != 'verified' or receipt.get('index') != index or
                receipt.get('schema') != 'gspa.price-full-export-job/1' or receipt.get('job') != job):
            raise ValueError(f'job receipt does not match plan: {part}')
        inputs = {**plan['inputs'], str(plan_path): expected[str(plan_path)], job['matrix']: job['matrix_sha256']}
        if receipt['inputs'] != inputs:
            raise ValueError(f'job input identity differs from plan: {part}')
        for path, sha in inputs.items():
            bind(path, sha)
        if receipt['claims'] != job['expected_claims'] or receipt['matrix_rows'] != job['matrix_rows']:
            raise ValueError(f'job counts differ from plan: {part}')
        if set(receipt['artifacts']) != required:
            raise ValueError(f'incomplete or unexpected artifact inventory: {part}')
        verification = json.loads((part / 'matrix-verification.json').read_text())
        if (verification.get('status') != 'verified' or
                verification.get('schema') != 'gspa.price-part-matrix-verification/1' or
                verification.get('predictor') != job['predictor'] or verification.get('assembly') != job['assembly'] or
                verification.get('claims_compared') != job['expected_claims'] or
                verification.get('matrix_rows') != job['matrix_rows']):
            raise ValueError(f'matrix comparison does not match plan: {part}')
        if verification['inputs'].get(job['matrix']) != job['matrix_sha256']:
            raise ValueError(f'matrix source identity differs: {part}')
        for name in required - {'matrix-verification.json'}:
            bound = [sha for path, sha in verification['inputs'].items() if Path(path).name == name]
            if bound != [receipt['artifacts'][name]]:
                raise ValueError(f'matrix comparison and published artifact differ: {part / name}')
        for name, sha in receipt['artifacts'].items():
            bind(part / name, sha)
        bind(receipt_path, digest(receipt_path))
        parts.append(dict(index=index, predictor=job['predictor'], assembly=job['assembly'],
                          claims=receipt['claims'], matrix_rows=receipt['matrix_rows'],
                          receipt_sha256=expected[str(receipt_path.resolve())]))
    total = sum(p['claims'] for p in parts)
    if total != plan['expected_claims']:
        raise ValueError('collection total differs from plan')
    # Reading small receipts is not proof that the bulk files still match them.
    for index, (path, sha) in enumerate(sorted(expected.items())):
        if digest(path) != sha:
            raise ValueError(f'published/input hash changed: {path}')
        if progress:
            progress(dict(files_verified=index + 1, files_total=len(expected), path=path))
    return dict(schema='gspa.price-export-collection-verification/1', status='verified',
                parts=parts, n_parts=len(parts), claims=total, hashes=expected,
                scope=__doc__.strip())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    if args.out.exists():
        raise SystemExit('refusing to overwrite an existing collection receipt')
    start = time.monotonic()
    report = verify(args.plan, lambda r: print(json.dumps(r), flush=True))
    report['elapsed_seconds'] = time.monotonic() - start
    with args.out.open('x') as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write('\n')
    print(json.dumps({k: report[k] for k in ('status', 'n_parts', 'claims', 'elapsed_seconds')}), flush=True)
