#!/usr/bin/env python3
"""GOAlpha-style Learning-to-Rank integrator over GSPA/CAFA6 component predictions.

Reproduces the core of the CAFA6 1st-place recipe: instead of max/avg combining
component predictors (which on our submission DESTROYED the MLP signal -> 0.359
vs MLP-alone 0.447), train a per-aspect gradient-boosted ranker over the
component scores + GO-term frequency + species, selected on a leak-free
no-knowledge validation.

Features per (protein, GO-term) candidate, per aspect:
  * one score column per component (diam, foldseek, clean, interpro, mlp);
    each component is first propagated to ancestors (max) so columns are
    comparable; absent -> 0.
  * ia(term)              (information accretion, the eval weight)
  * log10 train frequency of the term (from train_terms.tsv)
Label: 1 if term in the protein's propagated experimental GT for that aspect.

Evaluation is leak-free: GroupKFold by protein over the no-knowledge proteins,
out-of-fold LTR scores written per fold; final scoring is delegated to the
official `cafaeval` by the caller (this script writes a prediction TSV).

Usage:
  python train_ltr_integrator.py \
      --components components \
      --gt gt/gt_no.tsv --aspect-of-gt gt \
      --dag ~/Public/software/cafa6/go-dag.tsv \
      --ia ~/Public/software/cafa6/kaggle-official/IA.tsv \
      --train-terms ~/Public/software/cafa6/kaggle-official/train_terms.tsv \
      --taxon Test/testsuperset-taxon-list.tsv \
      --out ltr_out
"""
from __future__ import annotations

import argparse
import gzip
import os
import sys
import time
from collections import defaultdict

import numpy as np

COMPONENTS = ['diam', 'foldseek', 'clean', 'interpro', 'mlp']
ASP = {'F': 'MF', 'P': 'BP', 'C': 'CC'}
NS_OF_ROOT = {'GO:0003674': 'MF', 'GO:0008150': 'BP', 'GO:0005575': 'CC'}


def log(m):
    print(f'[{time.strftime("%H:%M:%S")}] {m}', file=sys.stderr, flush=True)


def opn(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def load_dag(path):
    anc = defaultdict(set)
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            c, _, a = line.rstrip('\n').partition('\t')
            if c and a:
                anc[c].add(a)
    for c in list(anc):
        anc[c].add(c)
    return anc


def load_aspect_from_dag_ia(ia_path, anc):
    """Aspect per term via ancestor containing a known root."""
    asp = {}
    for t in anc:
        for root, ns in NS_OF_ROOT.items():
            if root in anc[t]:
                asp[t] = ns
                break
    return asp


def load_component(path, keep, anc, aspect_of):
    """prot -> {term: score}, propagated to ancestors (max), only kept proteins."""
    out = defaultdict(dict)
    with opn(path) as fh:
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) < 3 or p[0] not in keep:
                continue
            try:
                s = float(p[2])
            except ValueError:
                continue
            d = out[p[0]]
            for t in anc.get(p[1], (p[1],)):
                if s > d.get(t, -1.0):
                    d[t] = s
    return out


def load_gt(path, anc, aspect_of):
    gt = defaultdict(lambda: defaultdict(set))   # prot -> aspect -> {terms}
    with open(path) as fh:
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) < 2:
                continue
            prot, term = p[0], p[1]
            for t in anc.get(term, (term,)):
                a = aspect_of.get(t)
                if a:
                    gt[prot][a].add(t)
    return gt


def main():
    global COMPONENTS
    ap = argparse.ArgumentParser()
    ap.add_argument('--components', default='components')
    ap.add_argument('--gt', required=True, help='no-knowledge GT (prot\\tGO)')
    ap.add_argument('--dag', required=True)
    ap.add_argument('--ia', required=True)
    ap.add_argument('--train-terms', required=True)
    ap.add_argument('--taxon', required=True, help='testsuperset-taxon-list.tsv (ID\\tSpecies) + per-protein taxon')
    ap.add_argument('--protein-taxon', help='optional prot\\ttaxon map; else species feature skipped')
    ap.add_argument('--out', default='ltr_out')
    ap.add_argument('--folds', type=int, default=5)
    ap.add_argument('--component-list', default='diam,foldseek,clean,interpro,mlp',
                    help='comma-separated component names = filenames (without .tsv.gz) '
                         'in --components dir; e.g. diam,foldseek,clean,interpro_lr,mlp')
    ap.add_argument('--model', choices=['xgb', 'logreg'], default='logreg',
                    help='logreg = linear over component scores only (cannot memorize '
                         'term base-rates -> honest integration); xgb = gradient boosting '
                         '(can leak term identity via IA/freq on in-distribution CV)')
    ap.add_argument('--save-model', default=None,
                    help='write the frozen per-aspect logreg (StandardScaler + coefficients) '
                         'to this JSON so it can be applied to novel proteins at inference '
                         '(requires --model logreg --features scores)')
    ap.add_argument('--features', choices=['scores', 'all'], default='scores',
                    help='scores = 5 component scores only; all = + IA + log-freq')
    args = ap.parse_args()

    from sklearn.model_selection import GroupKFold
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    xgb = None
    if args.model == 'xgb':
        import xgboost as xgb

    COMPONENTS = [c for c in args.component_list.split(',') if c]

    os.makedirs(args.out, exist_ok=True)
    log('loading DAG + aspects + IA + term freq ...')
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
    # term train-frequency
    freq = defaultdict(int)
    ntrain = 0
    with open(args.train_terms) as fh:
        next(fh, None)
        seen = set()
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) < 2:
                continue
            freq[p[1]] += 1
            if p[0] not in seen:
                seen.add(p[0]); ntrain += 1
    logfreq = {t: np.log10(c + 1) for t, c in freq.items()}

    # proteins of interest = GT proteins
    gt = load_gt(args.gt, anc, aspect_of)
    proteins = sorted(gt.keys())
    keep = set(proteins)
    log(f'GT proteins: {len(proteins):,}')

    log('loading components (propagated) ...')
    comp = {}
    for c in COMPONENTS:
        path = os.path.join(args.components, f'{c}.tsv.gz')
        if not os.path.exists(path):
            path = os.path.join(args.components, f'{c}.tsv')
        comp[c] = load_component(path, keep, anc, aspect_of)
        log(f'  {c}: {len(comp[c]):,} proteins')

    if args.save_model and (args.model != 'logreg' or args.features != 'scores'):
        raise SystemExit('--save-model requires --model logreg --features scores')

    # Build candidate (protein, term) rows per aspect
    results = {}
    frozen = {}  # aspect -> {mean, scale, coef, intercept} for the deployable model
    for aspect in ('MF', 'BP', 'CC'):
        log(f'=== aspect {aspect}: building candidates ===')
        rows_prot, rows_term, X, y = [], [], [], []
        for prot in proteins:
            # candidate terms = union of all component terms for this protein in this aspect
            cand = set()
            for c in COMPONENTS:
                for t in comp[c].get(prot, ()):
                    if aspect_of.get(t) == aspect:
                        cand.add(t)
            truth = gt[prot].get(aspect, set())
            for t in cand:
                feat = [comp[c].get(prot, {}).get(t, 0.0) for c in COMPONENTS]
                if args.features == 'all':
                    feat.append(ia.get(t, 0.0))
                    feat.append(logfreq.get(t, 0.0))
                X.append(feat)
                y.append(1 if t in truth else 0)
                rows_prot.append(prot)
                rows_term.append(t)
        X = np.asarray(X, dtype=np.float32)
        y = np.asarray(y, dtype=np.int8)
        groups = np.asarray(rows_prot)
        log(f'  rows={len(y):,} pos={int(y.sum()):,} ({100*y.mean():.1f}%)')

        # leak-free OOF predictions via GroupKFold by protein
        oof = np.zeros(len(y), dtype=np.float32)
        gkf = GroupKFold(n_splits=args.folds)
        for k, (tr, te) in enumerate(gkf.split(X, y, groups)):
            if args.model == 'xgb':
                dtr = xgb.DMatrix(X[tr], label=y[tr])
                dte = xgb.DMatrix(X[te])
                spw = max(1.0, (len(y[tr]) - y[tr].sum()) / max(1, y[tr].sum()))
                params = dict(objective='binary:logistic', eval_metric='logloss',
                              max_depth=6, eta=0.1, subsample=0.8, colsample_bytree=0.8,
                              scale_pos_weight=spw, nthread=8)
                bst = xgb.train(params, dtr, num_boost_round=200)
                oof[te] = bst.predict(dte)
            else:  # logreg over component scores only (cannot memorize term identity)
                sc = StandardScaler().fit(X[tr])
                lr = LogisticRegression(max_iter=1000, class_weight='balanced', C=1.0)
                lr.fit(sc.transform(X[tr]), y[tr])
                oof[te] = lr.predict_proba(sc.transform(X[te]))[:, 1]
        if args.model == 'logreg':
            # report learned component weights (avg over a full-data fit) for interpretability
            sc = StandardScaler().fit(X)
            lr = LogisticRegression(max_iter=1000, class_weight='balanced').fit(sc.transform(X), y)
            names = COMPONENTS + (['ia', 'logfreq'] if args.features == 'all' else [])
            log('  weights: ' + '  '.join(f'{n}={w:+.2f}' for n, w in zip(names, lr.coef_[0])))
            if args.save_model:
                frozen[aspect] = {
                    'mean': sc.mean_.astype(float).tolist(),
                    'scale': sc.scale_.astype(float).tolist(),
                    'coef': lr.coef_[0].astype(float).tolist(),
                    'intercept': float(lr.intercept_[0]),
                }
        results[aspect] = (rows_prot, rows_term, oof)

    if args.save_model:
        import json
        model = {'components': COMPONENTS, 'features': args.features, 'aspects': frozen}
        with open(args.save_model, 'w') as fh:
            json.dump(model, fh, indent=2)
        log(f'wrote frozen integrator -> {args.save_model}')

    # write LTR prediction TSV (per aspect concatenated) for cafaeval
    pred_dir = os.path.join(args.out, 'preds')
    os.makedirs(pred_dir, exist_ok=True)
    with open(os.path.join(pred_dir, 'ltr.tsv'), 'w') as out:
        for aspect in ('MF', 'BP', 'CC'):
            rp, rt, oof = results[aspect]
            for prot, term, s in zip(rp, rt, oof):
                if s >= 0.001:
                    out.write(f'{prot}\t{term}\t{s:.4f}\n')
    log(f'wrote {pred_dir}/ltr.tsv  -> score with cafaeval against the same GT')


if __name__ == '__main__':
    main()
