#!/usr/bin/env python3
"""Replay a sample of absence comparisons through the committed independent COBRA oracle.

The oracle is ``benchmark/v3/check_absence_cobra.py`` exactly as committed: it rebuilds every cell
from the manifest with COBRApy on GLPK and compares status and growth value, not just the class.
This wrapper only chooses the sample (by file order, never by outcome), re-executes the production
engine so a full result document exists for the oracle to check, and collects the verdicts.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing as mp
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(ROOT / 'benchmark' / 'v3'))

import absence_experiments as AE  # noqa: E402
import metabolic_experiments as ME  # noqa: E402
import check_absence_cobra as ORACLE  # noqa: E402


def sha(path) -> str:
    digest = hashlib.sha256()
    with open(path, 'rb') as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()


def one(task: tuple) -> dict:
    case, work = task
    directory = Path(work) / case['case_id'].replace('/', '_').replace(':', '_')
    directory.mkdir(parents=True, exist_ok=True)
    manifest = dict(case['manifest'])
    manifest_path = directory / 'manifest.json'
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True))
    try:
        model = ME.load_sbml_model(Path(manifest['model']))
        result = AE.run_absence_experiments(manifest, model)
        (directory / 'result.json').write_text(json.dumps(result, indent=2, sort_keys=True))
        verdict = ORACLE.compare(manifest_path, directory / 'result.json')
        return {'case_id': case['case_id'], 'status': 'checked', 'agrees': verdict['agrees'],
                'cells': len(verdict['comparisons']),
                'disagreements': [c for c in verdict['comparisons'] if not c['agrees']]}
    except Exception as exc:
        return {'case_id': case['case_id'], 'status': 'error', 'agrees': False,
                'failure': f'{type(exc).__name__}: {exc}'}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cases', type=Path, required=True)
    parser.add_argument('--every', type=int, default=40)
    parser.add_argument('--workers', type=int, default=6)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    cases = []
    with args.cases.open() as stream:
        for index, line in enumerate(stream):
            if line.strip() and index % args.every == 0:
                cases.append(json.loads(line))
    args.work.mkdir(parents=True, exist_ok=True)
    started = time.time()
    tasks = [(case, str(args.work)) for case in cases]
    if args.workers > 1:
        with mp.get_context('spawn').Pool(args.workers) as pool:
            rows = pool.map(one, tasks, chunksize=1)
    else:
        rows = [one(task) for task in tasks]
    report = {'schema': 'gspa.absence_truth/1:oracle_sample',
              'oracle': 'benchmark/v3/check_absence_cobra.py as committed (COBRApy on GLPK)',
              'oracle_sha256': sha(ROOT / 'benchmark/v3/check_absence_cobra.py'),
              'cases_sha256': sha(args.cases), 'sampling': f'every {args.every}th case by file order',
              'cases_checked': len(rows), 'cases_agreed': sum(1 for r in rows if r['agrees']),
              'cells_checked': sum(r.get('cells', 0) for r in rows),
              'errors': [r for r in rows if r['status'] == 'error'],
              'disagreeing_cases': [r for r in rows if r['status'] == 'checked' and not r['agrees']],
              'elapsed_s': time.time() - started}
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({k: report[k] for k in ('cases_checked', 'cases_agreed', 'cells_checked')}
                     | {'errors': len(report['errors']), 'disagreeing': len(report['disagreeing_cases'])}))
    return 0 if not report['errors'] and not report['disagreeing_cases'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
