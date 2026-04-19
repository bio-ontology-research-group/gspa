#!/usr/bin/env python3
"""Sharpened novelty scoring for culture dark-matter predictions.

For each top-K (K=5) prediction, compute cos over ALL panel EC centroids
(not just gap_ec). Bin predictions by how well ESM2 would handle them:

  truly_dark      max_cos_any_ec < 0.3       — no ESM2 signal anywhere
  remote_any      0.3 <= max < 0.5           — only remote ESM2 signal
  misleading_esm  cos_to_gap_ec < 0.4 AND max >= 0.7 AND argmax != gap_ec
                                              — ESM2 would predict wrong EC
  weak_partial    cos_to_gap_ec 0.4-0.7       — ESM2 weakly supports target
  esm_consistent  cos_to_gap_ec >= 0.7        — ESM2 would find it

Filters for candidate specificity:
  n_gaps_for_candidate   # distinct gap-reactions this candidate is top-K for

Output: scored TSV + per-bin summary + "validation candidates" shortlist.
"""
import argparse
import collections
import sys
from pathlib import Path

import numpy as np


def load_plm(plm_dir, tag):
    npy = Path(plm_dir) / f'{tag}_esm2t30.npy'
    idx = Path(plm_dir) / f'{tag}_esm2t30.index.tsv'
    arr = np.load(npy)
    pid_to_row = {}
    with open(idx) as f:
        f.readline()
        for line in f:
            pid, r = line.rstrip('\n').split('\t')
            pid_to_row[pid] = int(r)
    # L2 normalise
    norms = np.linalg.norm(arr, axis=1, keepdims=True) + 1e-9
    return arr / norms, pid_to_row


def load_centroids(path):
    arr = np.load(Path(path) / 'ec_centroids.npy')
    ec_to_row = {}
    row_to_ec = {}
    with open(Path(path) / 'ec_centroids.index.tsv') as f:
        f.readline()
        for line in f:
            ec, r, n = line.rstrip('\n').split('\t')
            r = int(r)
            ec_to_row[ec] = r
            row_to_ec[r] = ec
    # Already unit-norm from build script but re-norm for safety
    norms = np.linalg.norm(arr, axis=1, keepdims=True) + 1e-9
    return arr / norms, ec_to_row, row_to_ec


def classify(cos_gap, max_cos_any, argmax_ec, gap_ec):
    if max_cos_any < 0.3:
        return 'truly_dark'
    if max_cos_any < 0.5:
        return 'remote_any'
    # max_cos_any >= 0.5
    if cos_gap is not None and cos_gap >= 0.7:
        return 'esm_consistent'
    if (cos_gap is not None and cos_gap < 0.4
            and max_cos_any >= 0.7 and argmax_ec != gap_ec):
        return 'misleading_esm'
    if cos_gap is not None and 0.4 <= cos_gap < 0.7:
        return 'weak_partial'
    # gap_ec not in centroid set or cos_gap < 0.4 and argmax == gap_ec
    if cos_gap is None:
        # gap_ec unknown in centroid set; protein itself has moderate
        # ESM2 match to some other EC
        return 'gap_ec_uncovered_protein_known'
    return 'other'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--predictions', required=True)
    ap.add_argument('--plm-dir', required=True)
    ap.add_argument('--tag', required=True)
    ap.add_argument('--centroids-dir', required=True)
    ap.add_argument('--top-k', type=int, default=5,
                    help='Only score predictions with rank <= this.')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    panel, pid_to_row = load_plm(args.plm_dir, args.tag)
    centroids, ec_to_row, row_to_ec = load_centroids(args.centroids_dir)
    print(f'[info] panel {panel.shape}, '
          f'centroids {centroids.shape}', file=sys.stderr)

    # Full cos matrix (n_proteins, n_ecs) — small enough to compute once.
    sim_matrix = panel @ centroids.T  # (n_prot, n_ec)
    print(f'[info] sim matrix {sim_matrix.shape}', file=sys.stderr)

    # First pass: count distinct gaps per candidate (at rank<=K)
    rows = []
    cand_gaps = collections.defaultdict(set)
    with open(args.predictions) as f:
        h = f.readline().rstrip('\n').split('\t')
        for line in f:
            p = line.rstrip('\n').split('\t')
            d = dict(zip(h, p))
            if int(d.get('rank', 99)) > args.top_k:
                continue
            rows.append(d)
            cand_gaps[d['candidate']].add((d['gap_pathway'], d['gap_rxn']))

    # Score each prediction
    cols = ['culture', 'gap_pathway', 'gap_rxn', 'gap_ec', 'rank',
             'candidate', 'density', 'n_anchors',
             'n_gaps_for_candidate',
             'cos_to_gap_ec', 'max_cos_any_ec', 'argmax_ec',
             'bin']
    bin_counts = collections.Counter()
    rank_bin_counts = collections.Counter()
    per_bin_examples = collections.defaultdict(list)

    with open(args.out, 'w') as fout:
        fout.write('\t'.join(cols) + '\n')
        for d in rows:
            cand = d['candidate']
            ec = d['gap_ec']
            if ec.startswith('EC:'):
                ec = ec[3:]
            crow = pid_to_row.get(cand)
            if crow is None:
                continue
            sim_row = sim_matrix[crow]
            argmax_idx = int(np.argmax(sim_row))
            max_cos = float(sim_row[argmax_idx])
            argmax_ec = row_to_ec.get(argmax_idx, '')
            erow = ec_to_row.get(ec)
            cos_gap = (float(sim_row[erow]) if erow is not None else None)
            n_gaps = len(cand_gaps.get(cand, set()))
            bin_ = classify(cos_gap, max_cos, argmax_ec, ec)
            bin_counts[bin_] += 1
            rank_bin_counts[(bin_, d['rank'])] += 1
            if len(per_bin_examples[bin_]) < 3 and d['rank'] == '1':
                per_bin_examples[bin_].append(d)
            fout.write('\t'.join([
                d['culture'], d['gap_pathway'], d['gap_rxn'], ec,
                d['rank'], cand, d['density'], d['n_anchors'],
                str(n_gaps),
                'NA' if cos_gap is None else f'{cos_gap:.4f}',
                f'{max_cos:.4f}', argmax_ec,
                bin_,
            ]) + '\n')

    print(f'\n[{args.tag}] bin counts (top-{args.top_k} predictions):',
          file=sys.stderr)
    for b, n in sorted(bin_counts.items(), key=lambda x: -x[1]):
        print(f'  {b:35s} {n}', file=sys.stderr)
    print(f'[{args.tag}] rank distribution within bins:',
          file=sys.stderr)
    for b in bin_counts:
        ranks_in_bin = {r: rank_bin_counts.get((b, str(r)), 0)
                         for r in range(1, args.top_k + 1)}
        print(f'  {b:35s} rank1={ranks_in_bin[1]:>5d} '
              f'rank2={ranks_in_bin[2]:>5d} rank3={ranks_in_bin[3]:>5d} '
              f'rank4={ranks_in_bin[4]:>5d} rank5={ranks_in_bin[5]:>5d}',
              file=sys.stderr)


if __name__ == '__main__':
    main()
