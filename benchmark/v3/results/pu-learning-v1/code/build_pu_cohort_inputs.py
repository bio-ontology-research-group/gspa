#!/usr/bin/env python3
"""Batch manifest for a production replay of the 18 label-rich cohort genomes on Shaheen with a given learned
model and claim regime (plan `tasks/v3-row4-pu-learning-2026-11-02.md`, P1 and metric 9).

Mirrors the archived family-excluded verification command (`family-excluded-training-v1/default/command.json`:
`gspa integrate --batch-manifest ... --cohort --cohort-metadata ... --cohort-orthogroups ... --cohort-holdout CLADE
--feature-export-dir ...`) with Shaheen paths; the only additions are `--learned-model <artifact>` and the claim
sub-directory (`homology` or `family_excluded`). No label is read. Python 3.7 compatible (Shaheen login node).
"""
import argparse
import csv
import json
import os


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cohort', required=True, help='training-cohort bundle (batch_manifest.tsv, baseline/claims/<regime>, layout, cohort_inputs, reasoner-cache)')
    ap.add_argument('--claims-subdir', required=True, choices=('homology', 'family_excluded'))
    ap.add_argument('--claims-dir', default=None, help='override: directory holding <genome>.claims.jsonl (e.g. contract-bearing claims)')
    ap.add_argument('--genomes', nargs='*', default=None, help='subset of genome ids (default: all 18)')
    ap.add_argument('--out', required=True, help='job directory (manifest.tsv, inference/, features/)')
    args = ap.parse_args()
    cohort, out = args.cohort, args.out
    os.makedirs(os.path.join(out, 'inference'), exist_ok=True)
    os.makedirs(os.path.join(out, 'features'), exist_ok=True)
    with open(os.path.join(cohort, 'batch_manifest.tsv')) as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))
    columns = list(rows[0])
    keep = []
    for row in rows:
        gid = row['genome_id']
        if args.genomes and gid not in args.genomes:
            continue
        base = dict(row)
        base['claims_path'] = os.path.join(args.claims_dir or os.path.join(cohort, 'baseline', 'claims', args.claims_subdir), gid + '.claims.jsonl')
        base['genome_layout_path'] = os.path.join(cohort, 'layout', gid + '.layout.tsv')
        base['out_path'] = os.path.join(out, 'inference', gid + '.integrated.tsv')
        for key in ('claims_path', 'genome_layout_path'):
            assert os.path.isfile(base[key]), base[key]
        keep.append(base)
    if args.genomes and len(keep) != len(args.genomes):
        raise SystemExit('unknown genome id among %s' % args.genomes)
    with open(os.path.join(out, 'manifest.tsv'), 'w') as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter='\t', restval='-', lineterminator='\n')
        w.writeheader(); w.writerows(keep)
    with open(os.path.join(out, 'manifest.json'), 'w') as fh:
        json.dump({'genomes': [r['genome_id'] for r in keep], 'claims_subdir': args.claims_subdir, 'claims_dir': args.claims_dir,
                   'cohort': cohort}, fh, indent=1, sort_keys=True)
    print('%d genomes; manifest written to %s' % (len(keep), out))


if __name__ == '__main__':
    main()
