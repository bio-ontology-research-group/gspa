#!/usr/bin/env python3
"""One bounded, immutable full-matrix export/verification job for a frozen plan.

Writes bulk output to node-local scratch. Publishes only a verified, hash-matching
copy to persistent storage. Failed jobs write a failure receipt and retain logs;
they never create a completed destination or fall back to compact inputs.
"""
import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import time

try:
    from . import prepare_price_context as C, verify_price_part as V
except ImportError:
    import prepare_price_context as C
    import verify_price_part as V


def save(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True)+'\n')


def run(plan_path, index, scratch):
    plan_path, scratch = Path(plan_path).resolve(), Path(scratch).resolve()
    plan = json.loads(plan_path.read_text()); job = plan['jobs'][index]
    predictor, assembly = job['predictor'], job['assembly']
    root = Path(plan['out']); root.mkdir(parents=True, exist_ok=True)
    final = root/'parts'/predictor/assembly
    # Refuse stale source, layout, target or code before executing imported logic.
    identity = dict(plan['inputs'])
    identity[str(plan_path)] = V.digest(plan_path)
    identity[job['matrix']] = job['matrix_sha256']
    jobdir = scratch/f'price-export-{index}'
    started = time.time()
    receipt = {'schema': 'gspa.price-full-export-job/1', 'job': job, 'inputs': identity,
               'slurm_job_id': os.environ.get('SLURM_JOB_ID'), 'index': index,
               'packages': {p: importlib.metadata.version(p) for p in ('numpy', 'pandas')}}
    try:
        if final.exists(): raise ValueError(f'completed destination already exists: {final}')
        for filename, expected in identity.items():
            if V.digest(filename) != expected: raise ValueError(f'frozen input changed: {filename}')
        jobdir.mkdir(parents=True, exist_ok=False)
        rec = C.export_part(predictor, assembly, job['matrix'], job['matrix_sha256'],
                            plan['layout'], plan['targets'], str(jobdir), trust_pickle=True)
        receipt['target_coverage'] = rec['target_coverage']
        receipt['export_counts'] = rec['counts']
        part = Path(rec['part_dir'])
        if not rec['target_coverage']['verified']: raise ValueError('target statistics do not match original matrix')
        verification = V.verify(part, plan['layout'], plan['targets'], trust_pickle=True)
        if job.get('expected_claims') is not None and verification['claims_compared'] != job['expected_claims']:
            raise ValueError('full claim count differs from the original matrix audit')
        save(part/'matrix-verification.json', verification)
        expected = {p.name: V.digest(p) for p in part.iterdir() if p.is_file()}
        final.parent.mkdir(parents=True, exist_ok=True)
        staging = final.with_name('.'+assembly+'.publishing-'+str(os.getpid()))
        shutil.copytree(part, staging)
        copied = {p.name: V.digest(p) for p in staging.iterdir() if p.is_file()}
        if copied != expected: raise ValueError('published copy does not match verified scratch output')
        for filename, digest in identity.items():
            if V.digest(filename) != digest: raise ValueError(f'input changed during export: {filename}')
        receipt.update(status='verified', artifacts=expected, elapsed_seconds=time.time()-started,
                       claims=verification['claims_compared'], matrix_rows=verification['matrix_rows'])
        save(staging/'job-receipt.json', receipt)
        if final.exists(): raise ValueError('another job published this part')
        staging.rename(final)
        print(json.dumps({k: receipt[k] for k in ('status','index','claims','matrix_rows','elapsed_seconds')}), flush=True)
    except BaseException as exc:
        receipt.update(status='failed', error=f'{type(exc).__name__}: {exc}', elapsed_seconds=time.time()-started)
        failures = root/'failures'; failures.mkdir(exist_ok=True)
        save(failures/f'{index}-{os.environ.get("SLURM_JOB_ID", os.getpid())}.json', receipt)
        raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', required=True, type=Path)
    parser.add_argument('--index', required=True, type=int)
    parser.add_argument('--scratch', required=True, type=Path)
    args = parser.parse_args()
    run(args.plan, args.index, args.scratch)
