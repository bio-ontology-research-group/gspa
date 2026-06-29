#!/usr/bin/env python3
"""DeepGO-PlusPlus-Light — pick the best **no-GPU** component panel.

DeepGO-PlusPlus's strongest single components are the PLM heads (`mlp` ESM2-650M,
`prostt5`, `esm2_3b`), which need a GPU transformer at inference. DG++-Light drops
all of them and keeps only the components that run on CPU — alignment/search/
lookup signals: DIAMOND BLAST-KNN (`diam`), FoldSeek structure-KNN (`foldseek`,
over precomputed AlphaFold-DB structures — no GPU folding), InterProScan domain
hits both raw (`interpro`) and as a logistic model (`interpro_lr`), STRING Net-KNN
(`net`), and the BM25 literature channel (`lit`, optional). These are exactly the
components the CAFA6 winner found *generalise best on the private test set*.

This script evaluates several candidate no-GPU panels on the no-knowledge GT
(novel proteins) with the same OOF-by-protein + official-cafaeval protocol as
`ablation.py`, so the best CPU panel can be frozen with `train_integrator.py
--save-model` as `models/deepgo_plusplus_light.json` (the Groovy predictor reads
its component list from the JSON, so no code change is needed to ship it).

Usage:
  python eval_light.py --components components --gt gt/gt_no.tsv \
      --dag ~/Public/software/cafa6/go-dag.tsv \
      --ia  ~/Public/software/cafa6/kaggle-official/IA.tsv \
      --obo ~/Public/software/cafa6/go.obo \
      --train-terms ~/Public/software/cafa6/kaggle-official/train_terms.tsv \
      --out light_out
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from collections import defaultdict

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
from train_integrator import load_dag, load_aspect_from_dag_ia, load_component, load_gt  # noqa: E402
from ablation import oof_predictions, write_preds, score_dir_with_cafaeval, log  # noqa: E402

# Candidate no-GPU panels (label -> component list). All components here run on
# CPU only; the GPU PLM heads (mlp/prostt5/esm2_3b) are deliberately excluded and
# *replaced* by `cnn` — a CPU 1D-CNN over sequence (build_cnn_component.py).
LIGHT_PANELS = [
    # without the CNN (the "drop-only" baseline — why a sequence model is needed)
    ('diam+net',                       ['diam', 'net']),
    ('diam+foldseek+interpro_lr+net',  ['diam', 'foldseek', 'interpro_lr', 'net']),
    ('diam+foldseek+interpro+interpro_lr+net',
     ['diam', 'foldseek', 'interpro', 'interpro_lr', 'net']),
    ('diam+foldseek+interpro_lr+net+lit',
     ['diam', 'foldseek', 'interpro_lr', 'net', 'lit']),
    # with the CPU 1D-CNN replacing the PLM heads
    ('cnn (standalone)',               ['cnn']),
    ('cnn+diam+net',                   ['cnn', 'diam', 'net']),
    ('cnn+diam+foldseek+interpro+net', ['cnn', 'diam', 'foldseek', 'interpro', 'net']),
    ('cnn+diam+foldseek+interpro_lr+net',
     ['cnn', 'diam', 'foldseek', 'interpro_lr', 'net']),
    ('cnn+diam+foldseek+interpro+interpro_lr+net',
     ['cnn', 'diam', 'foldseek', 'interpro', 'interpro_lr', 'net']),
    ('cnn+diam+foldseek+interpro+interpro_lr+net+lit',
     ['cnn', 'diam', 'foldseek', 'interpro', 'interpro_lr', 'net', 'lit']),
]
CPU_COMPONENTS = ['cnn', 'diam', 'foldseek', 'interpro', 'interpro_lr', 'net', 'lit']


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--components', default='components')
    ap.add_argument('--gt', default='gt/gt_no.tsv')
    ap.add_argument('--dag', required=True)
    ap.add_argument('--ia', required=True)
    ap.add_argument('--obo', required=True)
    ap.add_argument('--train-terms', required=True)
    ap.add_argument('--out', default='light_out')
    ap.add_argument('--folds', type=int, default=5)
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    log('loading DAG / aspect / IA / freq ...')
    anc = load_dag(args.dag)
    aspect_of = load_aspect_from_dag_ia(args.ia, anc)
    ia = {}
    with open(args.ia) as fh:
        for line in fh:
            g, _, v = line.rstrip('\n').partition('\t')
            try:
                ia[g] = float(v)
            except ValueError:
                pass
    freq = defaultdict(int)
    with open(args.train_terms) as fh:
        next(fh, None)
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) >= 2:
                freq[p[1]] += 1
    logfreq = {t: np.log10(c + 1) for t, c in freq.items()}

    gt = load_gt(args.gt, anc, aspect_of)
    proteins = sorted(gt.keys())
    keep = set(proteins)
    log(f'no-knowledge GT proteins: {len(proteins):,}')

    log('loading CPU components once ...')
    comp = {}
    for c in CPU_COMPONENTS:
        path = os.path.join(args.components, f'{c}.tsv.gz')
        if not os.path.exists(path):
            path = os.path.join(args.components, f'{c}.tsv')
        if not os.path.exists(path):
            log(f'  WARN missing {c} -> skipped')
            continue
        comp[c] = load_component(path, keep, anc, aspect_of)
        log(f'  {c}: {len(comp[c]):,} proteins')

    preds_dir = os.path.join(args.out, 'preds')
    res_dir = os.path.join(args.out, 'results')
    subprocess.run(['rm', '-rf', preds_dir, res_dir], check=False)
    os.makedirs(preds_dir)
    fname_to_label = {}
    for i, (label, clist) in enumerate(LIGHT_PANELS):
        clist = [c for c in clist if c in comp]
        if not clist:
            continue
        rows = oof_predictions(comp, proteins, gt, aspect_of, clist, 'logreg',
                               'scores', ia, logfreq, folds=args.folds)
        fn = f'{i:02d}.tsv'
        write_preds(rows, os.path.join(preds_dir, fn))
        fname_to_label[fn] = label
        log(f'  {label}: preds written')

    log('scoring all panels with one cafaeval pass ...')
    scores = score_dir_with_cafaeval(preds_dir, res_dir, args.obo, args.gt, args.ia)

    rows = []
    for fn, label in fname_to_label.items():
        mf, bp, cc, mean = scores.get(fn, (0.0, 0.0, 0.0, 0.0))
        rows.append((label, mf, bp, cc, mean))
    rows.sort(key=lambda r: -r[4])

    tsv = os.path.join(args.out, 'light_results.tsv')
    with open(tsv, 'w') as fh:
        fh.write('panel\tMF\tBP\tCC\tmean_fw\n')
        for label, mf, bp, cc, mean in rows:
            fh.write(f'{label}\t{mf:.3f}\t{bp:.3f}\t{cc:.3f}\t{mean:.4f}\n')
    log(f'wrote {tsv}')

    print('\n## DeepGO-PlusPlus-Light — no-GPU panels (no-knowledge IA-weighted f_w)\n')
    print('| panel | MF | BP | CC | mean |')
    print('|---|---|---|---|---|')
    for label, mf, bp, cc, mean in rows:
        print(f'| {label} | {mf:.3f} | {bp:.3f} | {cc:.3f} | **{mean:.4f}** |')
    print(f'\nBest no-GPU panel: **{rows[0][0]}**  (mean f_w {rows[0][4]:.4f})')


if __name__ == '__main__':
    main()
