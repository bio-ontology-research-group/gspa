#!/usr/bin/env python3
"""LR-InterPro (GOAlpha component): per-GO-term logistic regression over the
binary InterPro feature vector. Trained on the pre-t0 SwissProt-exp population
(no leakage), applied to test proteins.

Replaces our weak raw InterPro->GO mapping (no-knowledge f_w 0.165) with the
GOAlpha recipe (~0.385).

Inputs (TSV: protein \\t IPR1;IPR2;... \\t GO1;GO2;... ; GO column optional for test):
  --train  train_ipr.tsv   (interpros + propagated GO labels)
  --test   test_ipr.tsv    (interpros; GO col ignored)
Outputs:
  --out    interpro_lr.tsv  (protein \\t GO \\t score) for the test proteins

Per-aspect, per-term liblinear LR over a shared sparse binary InterPro matrix;
only terms with >= --min-pos positive training proteins are modelled.
"""
from __future__ import annotations

import argparse
import sys
import time
from collections import defaultdict

import numpy as np
from scipy import sparse


def log(m):
    print(f'[{time.strftime("%H:%M:%S")}] {m}', file=sys.stderr, flush=True)


def read_rows(path):
    prots, iprs, gos = [], [], []
    with open(path) as fh:
        for line in fh:
            p = line.rstrip('\n').split('\t')
            prots.append(p[0])
            iprs.append([x for x in (p[1].split(';') if len(p) > 1 and p[1] else []) if x])
            gos.append([x for x in (p[2].split(';') if len(p) > 2 and p[2] else []) if x])
    return prots, iprs, gos


def aspect_of(obo):
    ns = {}
    cur = None
    NS = {'molecular_function': 'MF', 'biological_process': 'BP', 'cellular_component': 'CC'}
    for line in open(obo):
        line = line.rstrip('\n')
        if line == '[Term]':
            cur = None
        elif line.startswith('id: GO:'):
            cur = line[4:]
        elif line.startswith('namespace:') and cur:
            ns[cur] = NS.get(line.split(': ', 1)[1])
    return ns


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--train', required=True)
    ap.add_argument('--test', required=True)
    ap.add_argument('--obo', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--min-pos', type=int, default=10)
    ap.add_argument('--threads', type=int, default=8)
    args = ap.parse_args()

    from sklearn.linear_model import LogisticRegression
    from joblib import Parallel, delayed

    asp = aspect_of(args.obo)
    log('reading train ...')
    tr_p, tr_ipr, tr_go = read_rows(args.train)
    log('reading test ...')
    te_p, te_ipr, _ = read_rows(args.test)

    # IPR feature vocab from train
    feat = {}
    for iprs in tr_ipr:
        for f in iprs:
            if f not in feat:
                feat[f] = len(feat)
    nF = len(feat)
    log(f'IPR features: {nF:,}; train {len(tr_p):,}; test {len(te_p):,}')

    def build_X(ipr_lists):
        rows, cols = [], []
        for i, iprs in enumerate(ipr_lists):
            for f in iprs:
                j = feat.get(f)
                if j is not None:
                    rows.append(i); cols.append(j)
        data = np.ones(len(rows), dtype=np.float32)
        return sparse.csr_matrix((data, (rows, cols)), shape=(len(ipr_lists), nF))

    Xtr = build_X(tr_ipr)
    Xte = build_X(te_ipr)

    # term -> positive train indices
    term_pos = defaultdict(list)
    for i, gos in enumerate(tr_go):
        for g in gos:
            term_pos[g].append(i)
    terms = [t for t, pos in term_pos.items() if len(pos) >= args.min_pos and asp.get(t)]
    log(f'terms with >= {args.min_pos} pos: {len(terms):,}')

    ntr = len(tr_p)

    def fit_predict(t):
        y = np.zeros(ntr, dtype=np.int8)
        y[term_pos[t]] = 1
        try:
            lr = LogisticRegression(solver='liblinear', C=1.0, class_weight='balanced')
            lr.fit(Xtr, y)
            p = lr.predict_proba(Xte)[:, 1]
        except Exception:
            return t, None
        return t, p

    log('training per-term LR (parallel) ...')
    res = Parallel(n_jobs=args.threads, batch_size=16)(delayed(fit_predict)(t) for t in terms)

    log('writing predictions ...')
    n = 0
    with open(args.out, 'w') as out:
        for t, p in res:
            if p is None:
                continue
            for i in np.nonzero(p >= 0.01)[0]:
                out.write(f'{te_p[i]}\t{t}\t{p[i]:.4f}\n')
                n += 1
    log(f'wrote {n:,} predictions -> {args.out}')


if __name__ == '__main__':
    main()
