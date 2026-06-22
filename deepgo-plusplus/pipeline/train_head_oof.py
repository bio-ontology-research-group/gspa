#!/usr/bin/env python3
"""K-fold OUT-OF-FOLD PLM head scores over the pre-t0 train population.

For the pre-t0-validation integrator (CAFA6 (b)): the shipped integrator was
tuned on no-knowledge TEST labels; to make it blind-CAFA-faithful we instead
train it on a pre-t0 validation set. This produces the clean (out-of-sample)
PLM head scores for that set: K-fold over the train proteins, each fold's
predictions come from a head trained on the other folds. Mirrors train_head.py
(per-aspect MLP, frequency-weighted BCE) but writes OOF train predictions
instead of test predictions.

  --train-emb train_<model>.npz   --labels train_ipr.tsv   --dag go-dag.tsv
  --kfold 5   --out head_<model>_trainoof.tsv   (protein\\tGO\\tscore)
"""
import argparse, sys, time
from collections import defaultdict
import numpy as np
import torch, torch.nn as nn


def log(m): print(f'[{time.strftime("%H:%M:%S")}] {m}', file=sys.stderr, flush=True)

ASP_ROOT = {'MF': 'GO:0003674', 'BP': 'GO:0008150', 'CC': 'GO:0005575'}


def load_labels(path, aspect_of, aspect):
    y = defaultdict(set)
    with open(path) as fh:
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) < 3 or not p[2]:
                continue
            for g in p[2].split(';'):
                if aspect_of.get(g) == aspect:
                    y[p[0]].add(g)
    return y


def aspect_map(dag_path):
    anc = defaultdict(set)
    with open(dag_path) as fh:
        for line in fh:
            if line.startswith('#'): continue
            c, _, a = line.rstrip('\n').partition('\t')
            if c and a: anc[c].add(a)
    inv = {v: k for k, v in ASP_ROOT.items()}
    out = {}
    for t, A in anc.items():
        for root, asp in inv.items():
            if root in A: out[t] = asp; break
    return out


class MLP(nn.Module):
    def __init__(self, d, n, hidden=1024):
        super().__init__()
        self.net = nn.Sequential(nn.Linear(d, hidden), nn.GELU(), nn.Dropout(0.3),
                                 nn.Linear(hidden, n))
    def forward(self, x): return self.net(x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--train-emb', required=True)
    ap.add_argument('--labels', required=True)
    ap.add_argument('--dag', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--kfold', type=int, default=5)
    ap.add_argument('--min-pos', type=int, default=10)
    ap.add_argument('--epochs', type=int, default=30)
    ap.add_argument('--batch', type=int, default=512)
    ap.add_argument('--seed', type=int, default=0)
    args = ap.parse_args()
    dev = 'cuda' if torch.cuda.is_available() else 'cpu'
    rng = np.random.RandomState(args.seed)

    aspect_of = aspect_map(args.dag)
    tr = np.load(args.train_emb, allow_pickle=True)
    tr_ids = list(tr['ids']); tr_emb = tr['emb'].astype(np.float32)
    log(f'train {tr_emb.shape}  dev={dev}  kfold={args.kfold}')

    fout = open(args.out, 'w')
    for aspect in ('MF', 'BP', 'CC'):
        ylab = load_labels(args.labels, aspect_of, aspect)
        freq = defaultdict(int)
        for p, gs in ylab.items():
            for g in gs: freq[g] += 1
        terms = sorted([t for t, c in freq.items() if c >= args.min_pos])
        if not terms:
            continue
        tix = {t: j for j, t in enumerate(terms)}
        # rows = train proteins that have any label in this aspect
        rows = [i for i, p in enumerate(tr_ids) if p in ylab]
        rows = np.array(rows)
        Xall = tr_emb[rows]
        Y = np.zeros((len(rows), len(terms)), dtype=np.float32)
        for r, ri in enumerate(rows):
            for g in ylab[tr_ids[ri]]:
                if g in tix: Y[r, tix[g]] = 1.0
        pos = Y.sum(0); w = np.clip(np.log1p(pos.max() / np.maximum(pos, 1)), 0.25, 16)
        posw = torch.tensor(w, device=dev)
        log(f'{aspect}: {len(rows)} prot, {len(terms)} terms')

        # k-fold OOF
        order = rng.permutation(len(rows))
        folds = np.array_split(order, args.kfold)
        for k in range(args.kfold):
            te_idx = folds[k]
            tr_idx = np.concatenate([folds[j] for j in range(args.kfold) if j != k])
            Xtr = torch.tensor(Xall[tr_idx], device=dev)
            Ytr = torch.tensor(Y[tr_idx], device=dev)
            Xte = torch.tensor(Xall[te_idx], device=dev)
            model = MLP(Xtr.shape[1], len(terms)).to(dev)
            opt = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
            lossf = nn.BCEWithLogitsLoss(pos_weight=posw)
            n = Xtr.shape[0]
            for ep in range(args.epochs):
                model.train(); perm = torch.randperm(n, device=dev)
                for b in range(0, n, args.batch):
                    bi = perm[b:b + args.batch]
                    opt.zero_grad()
                    loss = lossf(model(Xtr[bi]), Ytr[bi]); loss.backward(); opt.step()
            model.eval()
            with torch.no_grad():
                preds = torch.sigmoid(model(Xte)).cpu().numpy()
            for r, ridx in enumerate(te_idx):
                pid = tr_ids[rows[ridx]]
                row = preds[r]
                for j in np.nonzero(row >= 0.01)[0]:
                    fout.write(f'{pid}\t{terms[j]}\t{row[j]:.4f}\n')
            log(f'{aspect}: fold {k+1}/{args.kfold} done')
    fout.close(); log(f'done -> {args.out}')


if __name__ == '__main__':
    main()
