#!/usr/bin/env python3
"""Compare a CLI batch cohort run with a Nextflow cohort run genome by genome.

For every genome: integrated TSV bytes, cohort report, learned report (model hash and every candidate
record), stage ledger (names, statuses, reasons and change lists; elapsed seconds excluded), provenance
and suggestions. Feature exports are compared file by file when both routes wrote them. The result is a
JSON record with one entry per genome and an overall ``identical`` flag; nothing is summarised away.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b''):
            h.update(chunk)
    return h.hexdigest()


def ledger_view(path):
    d = json.loads(Path(path).read_text())
    stages = []
    for s in d['stages']:
        stages.append({k: v for k, v in s.items() if k != 'elapsed_seconds'})
    top = {k: v for k, v in d.items() if k not in ('stages',)}
    return {'top': top, 'stages': stages}


def rows(path):
    with open(path) as fh:
        return sum(1 for _ in fh) - 1


def compare_genome(gid, direct_dir, nf_dir):
    d_base = Path(direct_dir) / f'{gid}.integrated.tsv'
    n_base = Path(nf_dir) / gid / f'{gid}_integrated.tsv'
    record = {'genome_id': gid, 'files': {}, 'identical': True}
    for suffix in ('', '.cohort.json', '.learned.json', '.stages.json', '.quality.json', '.rescues.tsv'):
        d, n = Path(str(d_base) + suffix), Path(str(n_base) + suffix)
        entry = {'direct': d.exists(), 'nextflow': n.exists()}
        if d.exists() and n.exists():
            if suffix == '.stages.json':
                dv, nv = ledger_view(d), ledger_view(n)
                entry['stage_states_equal'] = [(s['name'], s['status'], s.get('reason')) for s in dv['stages']] == \
                                              [(s['name'], s['status'], s.get('reason')) for s in nv['stages']]
                entry['stage_changes_equal'] = [s.get('changes') for s in dv['stages']] == [s.get('changes') for s in nv['stages']]
                entry['top_level_equal'] = dv['top'] == nv['top']
                entry['stage_states'] = [(s['name'], s['status']) for s in dv['stages']]
                entry['equal'] = entry['stage_states_equal'] and entry['stage_changes_equal'] and entry['top_level_equal']
            elif suffix == '.learned.json':
                dj, nj = json.loads(d.read_text()), json.loads(n.read_text())
                entry['model_sha256_equal'] = dj.get('model', {}).get('sha256') == nj.get('model', {}).get('sha256')
                entry['candidates'] = len(dj.get('candidates', []))
                entry['candidates_equal'] = dj.get('candidates') == nj.get('candidates')
                entry['counts_equal'] = dj.get('counts') == nj.get('counts')
                entry['equal'] = sha(d) == sha(n)
                entry['bytes_equal'] = entry['equal']
                if not entry['equal']:
                    entry['equal'] = entry['model_sha256_equal'] and entry['candidates_equal'] and entry['counts_equal']
            else:
                entry['sha256'] = sha(d)
                entry['equal'] = sha(d) == sha(n)
                if suffix == '':
                    entry['rows'] = rows(d)
            record['identical'] &= bool(entry['equal'])
        elif d.exists() != n.exists():
            entry['equal'] = False
            record['identical'] = False
        record['files'][suffix or '.tsv'] = entry
    # provenance and suggestions: direct writes next to the TSV with 'enabled' naming; Nextflow per genome dir
    for name, d_candidates, n_candidates in (
            ('provenance', [Path(direct_dir) / f'{gid}.integrated.tsv.provenance.json', Path(direct_dir) / f'{gid}.provenance.json'],
             [Path(nf_dir) / gid / 'provenance.json', Path(nf_dir) / gid / f'{gid}_integrated.tsv.provenance.json']),
            ('suggestions', [Path(direct_dir) / f'{gid}.integrated.tsv.suggestions.tsv', Path(direct_dir) / f'{gid}.suggestions.tsv'],
             [Path(nf_dir) / gid / 'suggestions.tsv', Path(nf_dir) / gid / f'{gid}_integrated.tsv.suggestions.tsv'])):
        d = next((p for p in d_candidates if p.exists()), None)
        n = next((p for p in n_candidates if n_candidates and p.exists()), None)
        entry = {'direct': str(d) if d else None, 'nextflow': str(n) if n else None}
        if d and n:
            entry['equal'] = sha(d) == sha(n)
            record['identical'] &= entry['equal']
        elif bool(d) != bool(n):
            entry['equal'] = False
            record['identical'] = False
        record['files'][name] = entry
    # feature exports
    d_feat, n_feat = Path(direct_dir) / 'features' / gid, Path(nf_dir) / 'features' / gid
    if d_feat.is_dir() or n_feat.is_dir():
        d_files = sorted(p.name for p in d_feat.glob('*')) if d_feat.is_dir() else []
        n_files = sorted(p.name for p in n_feat.glob('*')) if n_feat.is_dir() else []
        entry = {'direct_files': d_files, 'nextflow_files': n_files, 'equal': d_files == n_files and bool(d_files)}
        if entry['equal']:
            entry['equal'] = all(sha(d_feat / f) == sha(n_feat / f) for f in d_files)
        record['identical'] &= entry['equal']
        record['files']['features'] = entry
    return record


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--direct', required=True)
    ap.add_argument('--nextflow', required=True, help='Nextflow published results directory (contains one directory per genome)')
    ap.add_argument('--manifest', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    with open(args.manifest) as fh:
        ids = [r['genome_id'] for r in csv.DictReader(fh, delimiter='\t')]
    records = [compare_genome(g, args.direct, args.nextflow) for g in ids]
    nf_dir = Path(args.nextflow)
    listing = sorted(str(p.relative_to(nf_dir)) for p in nf_dir.rglob('*') if p.is_file())[:400]
    result = {'schema': 'gspa.routes-closure-cohort-comparison/1', 'genomes': len(ids),
              'identical': all(r['identical'] for r in records), 'records': records,
              'nextflow_files_sample': listing}
    Path(args.out).write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'genomes': len(ids), 'identical': result['identical'],
                      'per_genome': {r['genome_id']: r['identical'] for r in records}}))


if __name__ == '__main__':
    main()
