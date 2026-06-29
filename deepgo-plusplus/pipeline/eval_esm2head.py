#!/usr/bin/env python3
"""Does adding an ESM2-35M hierarchy-aware (MCM) head to DG++-Light cpu_lean help?

Leak-free no-knowledge eval, reusing ablation.py's OOF-by-protein + official
cafaeval machinery. Evaluates:
  * each component standalone (incl. the new esm2_head vs the existing esm2_knn),
  * cpu_lean (the 6 shipped components),
  * cpu_lean + esm2_head (the addition),
  * cpu_lean with esm2_knn -> esm2_head (parametric head replacing the kNN).

Run on ws (cafaeval is CPU-heavy and lives there). `cafaeval` must be on PATH.
"""
import argparse, os, sys, time
from collections import defaultdict
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
from train_integrator import load_dag, load_aspect_from_dag_ia, load_component, load_gt
from ablation import oof_predictions, write_preds, score_dir_with_cafaeval, log

LEAN = ['diam', 'interpro', 'cnn', 'net_union', 'esm2_knn', 'proteinfer']


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--components', required=True)
    ap.add_argument('--gt', required=True)
    ap.add_argument('--dag', required=True)
    ap.add_argument('--ia', required=True)
    ap.add_argument('--obo', required=True)
    ap.add_argument('--train-terms', required=True)
    ap.add_argument('--out', default='esm2head_out')
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

    want = LEAN + ['esm2_head']
    comp = {}
    for c in want:
        path = os.path.join(args.components, f'{c}.tsv.gz')
        if not os.path.exists(path):
            path = os.path.join(args.components, f'{c}.tsv')
        if not os.path.exists(path):
            log(f'  WARN missing component {c} -> skipped')
            continue
        comp[c] = load_component(path, keep, anc, aspect_of)
        log(f'  {c}: {len(comp[c]):,} proteins')

    lean = [c for c in LEAN if c in comp]
    configs = []
    for c in lean + (['esm2_head'] if 'esm2_head' in comp else []):
        configs.append(('standalone', c, [c]))
    configs.append(('panel', 'cpu_lean (6)', lean))
    if 'esm2_head' in comp:
        configs.append(('panel', 'cpu_lean + esm2_head (7)', lean + ['esm2_head']))
        swap = [c for c in lean if c != 'esm2_knn'] + ['esm2_head']
        configs.append(('panel', 'cpu_lean knn->head', swap))

    preds_dir = os.path.join(args.out, 'P')
    os.makedirs(preds_dir, exist_ok=True)
    fname_to_cfg = {}
    for i, (grp, label, clist) in enumerate(configs):
        clist = [c for c in clist if c in comp]
        if not clist:
            continue
        log(f'OOF [{i+1}/{len(configs)}] {label}  ({"+".join(clist)})')
        rows = oof_predictions(comp, proteins, gt, aspect_of, clist, 'logreg',
                               'scores', ia, logfreq, folds=args.folds)
        fn = f'cfg{i:02d}.tsv'
        write_preds(rows, os.path.join(preds_dir, fn))
        fname_to_cfg[fn[:-4]] = (grp, label)

    log('scoring all configs with one cafaeval pass ...')
    res_dir = os.path.join(args.out, 'res')
    scores = score_dir_with_cafaeval(preds_dir, res_dir, args.obo, args.gt, args.ia)

    results = []
    for fn, (grp, label) in fname_to_cfg.items():
        mf, bp, cc, mean = scores.get(fn, (0, 0, 0, 0))
        results.append((grp, label, mf, bp, cc, mean))
    order = {'standalone': 0, 'panel': 1}
    results.sort(key=lambda r: (order.get(r[0], 9), -r[5]))
    print(f'\n===== ESM2-35M MCM head as a DG++-Light addition '
          f'(no-knowledge {os.path.basename(args.gt)}, IA-weighted f_w) =====')
    print('%-10s %-26s %6s  %5s %5s %5s' % ('group', 'config', 'mean', 'MF', 'BP', 'CC'))
    last = None
    for grp, label, mf, bp, cc, mean in results:
        if grp != last:
            print('-' * 62); last = grp
        print('%-10s %-26s %6.3f  %5.3f %5.3f %5.3f' % (grp, label, mean, mf, bp, cc))


if __name__ == '__main__':
    main()
