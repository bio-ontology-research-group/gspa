"""Execute one label-independent, hash-bound cohort search inside a Slurm array."""
import argparse
import json
import os
from pathlib import Path

import observe_gapsmith as observer


def selected(manifest, expected_hash, index):
    if observer.digest(manifest) != expected_hash:
        raise ValueError('cohort manifest changed')
    data = json.loads(manifest.read_text())
    if data.get('schema') != 'gspa.cohort-pathway-search/1':
        raise ValueError('unsupported cohort manifest')
    rows = data['genomes']
    ids = [r['genome_id'] for r in rows]
    if len(ids) != len(set(ids)) or not 0 <= index < len(rows):
        raise ValueError('duplicate genome IDs or invalid array index')
    row = rows[index]
    if row['taxonomy'] not in ('Bacteria', 'Archaea'):
        raise ValueError('unsupported taxonomy')
    for path, expected in data['sources'].items():
        if observer.digest(Path(path)) != expected:
            raise ValueError('frozen source changed: ' + path)
    if observer.digest(Path(row['proteome'])) != row['sha256']:
        raise ValueError('frozen proteome changed')
    return data, row


def run(manifest, expected_hash, index):
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('cohort observation requires Slurm')
    data, row = selected(manifest, expected_hash, index)
    root = Path(data['output_root']) / row['genome_id']
    root.mkdir(parents=True, exist_ok=False)
    record = dict(schema='gspa.cohort-pathway-execution/1', genome_id=row['genome_id'],
                  index=index, manifest=str(manifest), manifest_sha256=expected_hash,
                  job_id=os.environ['SLURM_JOB_ID'], array_job_id=os.environ.get('SLURM_ARRAY_JOB_ID'),
                  sources=data['sources'], status='running')
    receipt = root/'execution.json'
    receipt.write_text(json.dumps(record, indent=2, sort_keys=True)+'\n')
    try:
        observer.run(root, Path(data['reference']), Path(data['binary']), Path(data['diamond']),
                     row['sha256'], Path(row['proteome']), row['taxonomy'], data['threads'])
        selected(manifest, expected_hash, index)  # bind code, manifest and target after execution too
        record['status'] = 'observed'
        record['observation_sha256'] = observer.digest(root/'observation.json')
    except Exception as exc:
        record.update(status='failed', failure_type=type(exc).__name__, failure=str(exc))
        raise
    finally:
        receipt.write_text(json.dumps(record, indent=2, sort_keys=True)+'\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--manifest-sha256', required=True)
    p.add_argument('--index', type=int, required=True)
    a = p.parse_args()
    run(a.manifest.resolve(), a.manifest_sha256, a.index)
