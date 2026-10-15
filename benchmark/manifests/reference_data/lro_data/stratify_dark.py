#!/usr/bin/env python3
"""Stratify LRO results by max_pident_annot_homolog column in the
case file. Prints per-bucket hit@k / MRR for each results file.
"""
import argparse
import collections
from pathlib import Path


BUCKETS = [
    ('no_hit_pident=0',   0.0, 0.0),
    ('dark_(0,30]',       0.01, 30.0),
    ('twilight_(30,50]',  30.01, 50.0),
    ('near_(50,70]',      50.01, 70.0),
    ('close_(70,90]',     70.01, 90.0),
    ('very_close_(90,]',  90.01, 101.0),
]


def bucket_of(pid):
    for name, lo, hi in BUCKETS:
        if lo <= pid <= hi:
            return name
    return 'unk'


def load_tsv(path):
    rows = []
    with open(path) as f:
        header = f.readline().rstrip('\n').split('\t')
        for line in f:
            parts = line.rstrip('\n').split('\t')
            rows.append(dict(zip(header, parts)))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cases', required=True)
    ap.add_argument('--results', nargs='+', required=True)
    args = ap.parse_args()

    cases = load_tsv(args.cases)
    depth = {(r['protein_id'], r['reaction_id']):
              float(r['max_pident_annot_homolog'])
              for r in cases}
    print(f'[info] loaded {len(depth)} cases', file=sys.stderr) \
        if False else None

    header = (f'{"variant":36s}  {"bucket":22s}  {"n":>5}  '
              f'{"hit@1":>6}  {"hit@3":>6}  {"hit@5":>6}  {"MRR":>6}')
    print(header)

    for rp in args.results:
        rs = load_tsv(rp)
        buckets = collections.defaultdict(list)
        for r in rs:
            k = (r['protein_id'], r['reaction_id'])
            if k not in depth:
                continue
            b = bucket_of(depth[k])
            buckets[b].append(r)
            buckets['ALL'].append(r)
        name = Path(rp).stem.replace('results_', '')
        for b in ['ALL'] + [x[0] for x in BUCKETS]:
            v = buckets.get(b, [])
            if not v:
                continue
            n = len(v)
            ranks = [int(r.get('rank_of_p', '0') or 0) for r in v]
            h1 = sum(1 for x in ranks if x == 1) / n
            h3 = sum(1 for x in ranks if 1 <= x <= 3) / n
            h5 = sum(1 for x in ranks if 1 <= x <= 5) / n
            mrr = sum(1.0 / x for x in ranks if x > 0) / n
            print(f'{name:36s}  {b:22s}  {n:5d}  '
                  f'{h1:6.3f}  {h3:6.3f}  {h5:6.3f}  {mrr:6.3f}')
        print('')


if __name__ == '__main__':
    import sys
    main()
