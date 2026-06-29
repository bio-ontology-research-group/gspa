#!/usr/bin/env python3
"""DeepGO-PlusPlus ablation harness — which component / aggregator actually helps.

Reproduces (offline, on the faithful CAFA6 reconstruction) the kind of ablation
the CAFA6 winner reported: it quantifies the *marginal* contribution of every
component and compares the **logistic** aggregator against **gradient-boosted
trees** (the GOAlpha Learning-to-Rank choice).

All configs are scored identically: per-aspect candidate rows -> GroupKFold-by-
protein out-of-fold scores -> official `cafaeval` IA-weighted f_w on the SAME
no-knowledge ground truth. No-knowledge (novel proteins) is the metric the whole
project optimises — it is where CAFA6 is decided and where naive max-merge fails.

The component loaders + OOF logic are imported from / mirror `train_integrator.py`
(single source of truth for the integration math); this script only loops the
configs and loads every component once so ~30 configs cost one parse each.

Usage:
  python ablation.py --components components --gt gt/gt_no.tsv \
      --dag ~/Public/software/cafa6/go-dag.tsv \
      --ia  ~/Public/software/cafa6/kaggle-official/IA.tsv \
      --obo ~/Public/software/cafa6/go.obo \
      --train-terms ~/Public/software/cafa6/kaggle-official/train_terms.tsv \
      --out ablation_out
Outputs: <out>/ablation_results.tsv and a Markdown table on stdout.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from collections import defaultdict

import numpy as np

# Reuse the exact loaders the production trainer uses (no math drift).
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
from train_integrator import load_dag, load_aspect_from_dag_ia, load_component, load_gt  # noqa: E402


def log(m):
    print(f'[{time.strftime("%H:%M:%S")}] {m}', file=sys.stderr, flush=True)


def oof_predictions(comp, proteins, gt, aspect_of, comp_list, model, features,
                    ia, logfreq, folds=5):
    """Return list of (prot, term, oof_score) over all three aspects.

    Mirrors the per-aspect candidate build + GroupKFold OOF loop in
    train_integrator.main(), parameterised by component subset / model / features.
    """
    from sklearn.model_selection import GroupKFold
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    xgb = None
    if model == 'xgb':
        import xgboost as xgb

    out_rows = []
    for aspect in ('MF', 'BP', 'CC'):
        rows_prot, rows_term, X, y = [], [], [], []
        for prot in proteins:
            cand = set()
            for c in comp_list:
                for t in comp[c].get(prot, ()):
                    if aspect_of.get(t) == aspect:
                        cand.add(t)
            truth = gt[prot].get(aspect, set())
            for t in cand:
                feat = [comp[c].get(prot, {}).get(t, 0.0) for c in comp_list]
                if features == 'all':
                    feat.append(ia.get(t, 0.0))
                    feat.append(logfreq.get(t, 0.0))
                X.append(feat)
                y.append(1 if t in truth else 0)
                rows_prot.append(prot)
                rows_term.append(t)
        if not y:
            continue
        X = np.asarray(X, dtype=np.float32)
        y = np.asarray(y, dtype=np.int8)
        groups = np.asarray(rows_prot)
        oof = np.zeros(len(y), dtype=np.float32)
        gkf = GroupKFold(n_splits=folds)
        for tr, te in gkf.split(X, y, groups):
            # A weak standalone component can leave a fold all-negative; the
            # linear solver needs 2 classes. Fall back to the train base rate.
            if y[tr].min() == y[tr].max():
                oof[te] = float(y[tr].mean())
                continue
            if model == 'xgb':
                dtr = xgb.DMatrix(X[tr], label=y[tr])
                dte = xgb.DMatrix(X[te])
                spw = max(1.0, (len(y[tr]) - y[tr].sum()) / max(1, int(y[tr].sum())))
                params = dict(objective='binary:logistic', eval_metric='logloss',
                              max_depth=6, eta=0.1, subsample=0.8,
                              colsample_bytree=0.8, scale_pos_weight=spw, nthread=8)
                bst = xgb.train(params, dtr, num_boost_round=200)
                oof[te] = bst.predict(dte)
            else:
                sc = StandardScaler().fit(X[tr])
                lr = LogisticRegression(max_iter=1000, class_weight='balanced', C=1.0)
                lr.fit(sc.transform(X[tr]), y[tr])
                oof[te] = lr.predict_proba(sc.transform(X[te]))[:, 1]
        out_rows.extend(zip(rows_prot, rows_term, oof))
    return out_rows


def write_preds(rows, path):
    with open(path, 'w') as fh:
        for prot, term, s in rows:
            if s >= 0.001:
                fh.write(f'{prot}\t{term}\t{s:.4f}\n')


def score_dir_with_cafaeval(preds_dir, res_dir, obo, gt_path, ia_path):
    """One cafaeval call over a dir of prediction files (one obo parse for all).

    Returns {filename: (MF, BP, CC, mean)} from evaluation_best_f_w.tsv, where
    cafaeval reports the per-(file, namespace) threshold-max IA-weighted f_w.
    """
    subprocess.run(['cafaeval', obo, preds_dir, gt_path, '-ia', ia_path,
                    '-norm', 'cafa', '-prop', 'max', '-th_step', '0.01',
                    '-out_dir', res_dir],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    best = os.path.join(res_dir, 'evaluation_best_f_w.tsv')
    lines = open(best).read().splitlines()
    hdr = lines[0].split('\t')
    fni, ni, fi = hdr.index('filename'), hdr.index('ns'), hdr.index('f_w')
    per = defaultdict(dict)
    for r in lines[1:]:
        c = r.split('\t')
        per[c[fni]][c[ni]] = float(c[fi])
    out = {}
    for fn, m in per.items():
        mf = m.get('molecular_function', 0.0)
        bp = m.get('biological_process', 0.0)
        cc = m.get('cellular_component', 0.0)
        out[fn] = (mf, bp, cc, (mf + bp + cc) / 3.0)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--components', default='components')
    ap.add_argument('--gt', default='gt/gt_no.tsv')
    ap.add_argument('--dag', required=True)
    ap.add_argument('--ia', required=True)
    ap.add_argument('--obo', required=True)
    ap.add_argument('--train-terms', required=True)
    ap.add_argument('--out', default='ablation_out')
    ap.add_argument('--folds', type=int, default=5)
    ap.add_argument('--all-components',
                    default='diam,foldseek,clean,interpro,mlp,prostt5,net,lit,esm2_3b,interpro_lr')
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    all_comps = [c for c in args.all_components.split(',') if c]

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

    log('loading components once (propagated, filtered to GT) ...')
    comp = {}
    for c in all_comps:
        path = os.path.join(args.components, f'{c}.tsv.gz')
        if not os.path.exists(path):
            path = os.path.join(args.components, f'{c}.tsv')
        if not os.path.exists(path):
            log(f'  WARN missing component {c} -> skipped')
            continue
        comp[c] = load_component(path, keep, anc, aspect_of)
        log(f'  {c}: {len(comp[c]):,} proteins')
    have = [c for c in all_comps if c in comp]

    SHIP6 = [c for c in ['diam', 'foldseek', 'clean', 'interpro', 'mlp', 'prostt5'] if c in comp]
    PANEL7 = SHIP6 + (['net'] if 'net' in comp else [])

    # --- build the config battery ---
    configs = []  # (group, label, comp_list, model, features)

    # 1. per-component standalone (1-feature logreg == raw component f_w)
    for c in have:
        configs.append(('standalone', c, [c], 'logreg', 'scores'))

    # 2. leave-one-out from the recommended 7-comp panel (6 + net)
    configs.append(('loo', 'full-7 (6+net)', PANEL7, 'logreg', 'scores'))
    for c in PANEL7:
        configs.append(('loo', f'-{c}', [x for x in PANEL7 if x != c], 'logreg', 'scores'))

    # 3. cumulative build-up (best-first ordering)
    order = [c for c in ['mlp', 'prostt5', 'diam', 'foldseek', 'interpro', 'clean', 'net']
             if c in comp]
    for i in range(1, len(order) + 1):
        configs.append(('cumulative', '+'.join(order[:i]), order[:i], 'logreg', 'scores'))

    # 4. aggregator comparison on the shipped 6-comp panel
    configs.append(('aggregator', 'logreg / scores (shipped)', SHIP6, 'logreg', 'scores'))
    configs.append(('aggregator', 'xgb / scores (leak-free GBT)', SHIP6, 'xgb', 'scores'))
    configs.append(('aggregator', 'logreg / +IA+freq', SHIP6, 'logreg', 'all'))
    configs.append(('aggregator', 'xgb / +IA+freq (LEAK)', SHIP6, 'xgb', 'all'))

    # Generate OOF preds for every config into one dir, then score all in a
    # single cafaeval pass (one obo parse instead of ~30).
    preds_dir = os.path.join(args.out, 'preds')
    res_dir = os.path.join(args.out, 'results')
    subprocess.run(['rm', '-rf', preds_dir, res_dir], check=False)
    os.makedirs(preds_dir)
    fname_to_cfg = {}
    for i, (grp, label, clist, model, feats) in enumerate(configs):
        clist = [c for c in clist if c in comp]
        if not clist:
            continue
        t0 = time.time()
        rows = oof_predictions(comp, proteins, gt, aspect_of, clist, model, feats,
                               ia, logfreq, folds=args.folds)
        fn = f'{i:03d}.tsv'
        write_preds(rows, os.path.join(preds_dir, fn))
        fname_to_cfg[fn] = (grp, label)
        log(f'[{grp:11s}] {label:32s} preds written  ({time.time()-t0:.0f}s)')

    log('scoring all configs with one cafaeval pass ...')
    scores = score_dir_with_cafaeval(preds_dir, res_dir, args.obo, args.gt, args.ia)

    results = []
    for fn, (grp, label) in fname_to_cfg.items():
        mf, bp, cc, mean = scores.get(fn, (0.0, 0.0, 0.0, 0.0))
        results.append((grp, label, mf, bp, cc, mean))
        log(f'[{grp:11s}] {label:32s} MF={mf:.3f} BP={bp:.3f} CC={cc:.3f} mean={mean:.4f}')

    # --- write TSV + markdown ---
    tsv = os.path.join(args.out, 'ablation_results.tsv')
    with open(tsv, 'w') as fh:
        fh.write('group\tconfig\tMF\tBP\tCC\tmean_fw\n')
        for grp, label, mf, bp, cc, mean in results:
            fh.write(f'{grp}\t{label}\t{mf:.3f}\t{bp:.3f}\t{cc:.3f}\t{mean:.4f}\n')
    log(f'wrote {tsv}')

    print('\n## DeepGO-PlusPlus ablation (no-knowledge, IA-weighted f_w)\n')
    last_grp = None
    for grp, label, mf, bp, cc, mean in results:
        if grp != last_grp:
            print(f'\n### {grp}\n')
            print('| config | MF | BP | CC | mean |')
            print('|---|---|---|---|---|')
            last_grp = grp
        print(f'| {label} | {mf:.3f} | {bp:.3f} | {cc:.3f} | **{mean:.4f}** |')


if __name__ == '__main__':
    main()
