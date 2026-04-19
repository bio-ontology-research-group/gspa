#!/usr/bin/env python3
"""Score how "novel" each dark-matter prediction is relative to an ESM2
homology baseline.

For each culture, for each top-1 prediction (gap_ec, candidate):
  plm_cos = cos(ESM2(candidate), EC-centroid(gap_ec))

Also: per-candidate specificity = 1 / (distinct gaps predicted for).

Novelty is a prediction where:
  - context thinks this protein catalyses the reaction (density high,
    low n_gaps-per-candidate)
  - PLM cos to EC-centroid is LOW (ESM2 would not have found it)

Annotated output: culture gap_pathway gap_rxn gap_ec candidate density
n_anchors n_gaps_for_candidate plm_cos flag
  flag = {PLM_would_find | BORDERLINE | context_only_novel | no_plm_embed}
"""
import argparse
import collections
import sys
from pathlib import Path

import numpy as np


def load_npy_index(plm_dir, tag):
    npy = Path(plm_dir) / f'{tag}_esm2t30.npy'
    idx = Path(plm_dir) / f'{tag}_esm2t30.index.tsv'
    if not npy.exists() or not idx.exists():
        return None, {}
    arr = np.load(npy, mmap_mode='r')
    pid_to_row = {}
    with open(idx) as f:
        f.readline()
        for line in f:
            pid, r = line.rstrip('\n').split('\t')
            pid_to_row[pid] = int(r)
    return arr, pid_to_row


def load_ec_centroids(path):
    arr = np.load(Path(path) / 'ec_centroids.npy')
    ec_to_row = {}
    with open(Path(path) / 'ec_centroids.index.tsv') as f:
        f.readline()
        for line in f:
            ec, r, n = line.rstrip('\n').split('\t')
            ec_to_row[ec] = int(r)
    return arr, ec_to_row


def cos(a, b):
    na = np.linalg.norm(a) + 1e-9
    nb = np.linalg.norm(b) + 1e-9
    return float(np.dot(a, b) / (na * nb))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--predictions', required=True)
    ap.add_argument('--plm-dir', required=True,
                    help='dir with {tag}_esm2t30.npy / index.tsv')
    ap.add_argument('--tag', required=True)
    ap.add_argument('--centroids-dir', required=True)
    ap.add_argument('--plm-high', type=float, default=0.7,
                    help='cos > this: PLM would find it')
    ap.add_argument('--plm-low', type=float, default=0.4,
                    help='cos < this: novel (PLM would miss)')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    panel, pid_to_row = load_npy_index(args.plm_dir, args.tag)
    if panel is None:
        print(f'[error] no PLM for {args.tag}', file=sys.stderr)
        sys.exit(1)
    centroids, ec_to_row = load_ec_centroids(args.centroids_dir)
    print(f'[info] PLM shape {panel.shape}, '
          f'{len(ec_to_row)} EC centroids', file=sys.stderr)

    # First pass: count distinct gaps per candidate (for specificity).
    rows = []
    cand_gaps = collections.defaultdict(set)
    with open(args.predictions) as f:
        h = f.readline().rstrip('\n').split('\t')
        for line in f:
            p = line.rstrip('\n').split('\t')
            d = dict(zip(h, p))
            rows.append(d)
            if d.get('rank') == '1':
                cand_gaps[d['candidate']].add((d['gap_pathway'],
                                                 d['gap_rxn']))

    # Write annotated output
    cols = ['culture', 'gap_pathway', 'gap_rxn', 'gap_ec', 'rank',
             'candidate', 'density', 'n_anchors', 'n_gaps_for_candidate',
             'plm_cos', 'specificity', 'flag']
    counts = collections.Counter()
    with open(args.out, 'w') as fout:
        fout.write('\t'.join(cols) + '\n')
        for d in rows:
            if d.get('rank') != '1':
                continue
            cand = d['candidate']
            ec = d['gap_ec']
            if ec.startswith('EC:'):
                ec = ec[3:]
            plm_cos_val = None
            crow = pid_to_row.get(cand)
            erow = ec_to_row.get(ec)
            if crow is not None and erow is not None:
                plm_cos_val = cos(panel[crow], centroids[erow])
            n_gaps = len(cand_gaps.get(cand, set()))
            spec = 1.0 / max(n_gaps, 1)

            if plm_cos_val is None:
                flag = 'no_plm_embed'
            elif plm_cos_val >= args.plm_high:
                flag = 'PLM_would_find'
            elif plm_cos_val <= args.plm_low:
                flag = 'context_only_novel'
            else:
                flag = 'BORDERLINE'
            counts[flag] += 1
            fout.write('\t'.join([
                d['culture'], d['gap_pathway'], d['gap_rxn'], ec,
                d['rank'], cand, d['density'], d['n_anchors'],
                str(n_gaps),
                'NA' if plm_cos_val is None else f'{plm_cos_val:.4f}',
                f'{spec:.4f}', flag,
            ]) + '\n')
    for f, n in counts.most_common():
        print(f'  {f}: {n}', file=sys.stderr)
    print(f'[info] wrote {args.out}', file=sys.stderr)


if __name__ == '__main__':
    main()
