#!/usr/bin/env python3
"""Submit a cohort bundle to the web service and run it with the real worker (plan 2026-11-01, check A6 iii).

The bundle goes through the real submission endpoint (`POST /api/jobs`, kind=cohort) when FastAPI is importable;
otherwise through the same validation function the endpoint calls (`service.cohort.accept`), which is recorded.
The worker is `service.jobs.run` (dispatching to `service.cohort.run_cohort`). Prints the final job state.
"""
import argparse
import json
import os
import sys
import time
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--bundle', required=True)
    ap.add_argument('--jobs', required=True, help='GSPA_JOBS root for this run')
    ap.add_argument('--resources', required=True, help='GSPA_COHORT_RESOURCES JSON')
    ap.add_argument('--gspa', required=True, help='GSPA_BIN launcher')
    ap.add_argument('--threads', default='8')
    ap.add_argument('--go-obo', required=True)
    ap.add_argument('--feature-export', action='store_true')
    ap.add_argument('--job-id-file', required=True)
    args = ap.parse_args()
    os.environ.update(GSPA_JOBS=args.jobs, GSPA_COHORT_RESOURCES=args.resources, GSPA_BIN=args.gspa, GSPA_THREADS=args.threads,
                      GSPA_GO_OBO=args.go_obo, GSPA_FEATURE_EXPORT='1' if args.feature_export else '')
    os.environ.setdefault('GSPA_VIZ_PYTHON', sys.executable)
    from service import cohort, jobs
    jobs.ROOT.mkdir(parents=True, exist_ok=True)
    submitted_via = None
    try:
        from fastapi.testclient import TestClient
        from service import portal
        from unittest.mock import patch
        from types import SimpleNamespace
        client = TestClient(portal.app)
        with patch.object(jobs, 'ready', return_value=True), patch.object(jobs.shutil, 'disk_usage', return_value=SimpleNamespace(free=100 * 1024**3)):
            with open(args.bundle, 'rb') as fh:
                r = client.post('/api/jobs', data={'kind': 'cohort', 'preset': 'bacteria'}, files={'file': ('bundle.zip', fh)})
        if r.status_code != 202:
            raise SystemExit(f'submission refused: {r.status_code} {r.text}')
        job = r.json()
        submitted_via = 'POST /api/jobs (FastAPI TestClient)'
    except ImportError:
        import uuid
        job_id = uuid.uuid4().hex
        folder = jobs.ROOT / job_id
        folder.mkdir()
        (folder / 'cohort.zip').write_bytes(Path(args.bundle).read_bytes())
        accepted = cohort.accept(folder)
        job = jobs.state(folder, id=job_id, status='queued', stage='Waiting for a worker', created=time.time(), kind='cohort', preset='bacteria',
                         genomes=accepted['genomes'], sequences=accepted['proteins'], residues=0, analyses={}, gram='neg', filename='bundle.zip')
        submitted_via = 'service.cohort.accept + jobs.state (FastAPI not importable here)'
    folder = jobs.ROOT / job['id']
    Path(args.job_id_file).write_text(job['id'] + '\n')
    t = time.time()
    try:
        jobs.run(folder)
    except Exception as exc:
        jobs.state(folder, status='failed', stage='Job stopped', error=str(exc), finished=time.time())
    state = jobs.state(folder)
    state['submitted_via'] = submitted_via
    state['worker_seconds'] = round(time.time() - t, 1)
    (folder / 'route.json').write_text(json.dumps(state, indent=2) + '\n')
    print(json.dumps(state, indent=1))
    return 0 if state.get('status') == 'complete' else 1


if __name__ == '__main__':
    sys.exit(main())
