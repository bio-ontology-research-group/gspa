#!/usr/bin/env python3
"""Hierarchy-aware per-aspect PLM head — C-HMCNN Max-Constraint Module (+ controls).

A drop-in successor to `train_head.py` / `train_head_oof.py` that swaps the flat
BCE loss for a hierarchy-aware objective over the GO **is_a ∪ part_of** DAG. Same
encoder (per-aspect MLP over PLM embeddings), same per-aspect term vocab and
frequency pos-weights, so the ONLY changed variable vs the shipped heads is the
loss — a clean A/B for "does hierarchy-aware training help each component?".

Three losses (`--loss`):
  bce      flat BCEWithLogits + freq pos-weight                 (baseline = train_head.py)
  mcm      C-HMCNN Max-Constraint Module + MCLoss               (Giunchiglia & Lukasiewicz, NeurIPS'20)
  softreg  flat BCE + λ·mean(relu(h_child − h_parent))          (soft true-path penalty, the cheap control)

The MCM is implemented with `scatter_reduce(amax)` over the per-aspect descendant
edge list (ancestor bucket ← descendant), so it scales to thousands of terms
(no (B,n,n) materialisation as in the reference impl). The DAG is the project's
`go-dag.tsv` (child<TAB>ancestor over is_a ∪ part_of); edges are restricted to the
aspect's own term set (the head is per-aspect, so cross-aspect part_of edges drop
out naturally).

Modes (`--mode`):
  test   train on the full train pop, predict --test-emb        (standalone eval)
  oof    K-fold, write held-out train predictions               (integrator training set)

Usage:
  python train_head_hmcnn.py --loss mcm --mode test \
      --train-emb emb/train_prostt5.npz --test-emb emb/test_prostt5.npz \
      --labels train_ipr.tsv --dag go-dag.tsv --out preds/head_prostt5_mcm.tsv
"""
import argparse, os, sys, time
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


def load_dag(dag_path):
    """child -> set(proper ancestors), over is_a ∪ part_of (go-dag.tsv)."""
    anc = defaultdict(set)
    with open(dag_path) as fh:
        for line in fh:
            if line.startswith('#'): continue
            c, _, a = line.rstrip('\n').partition('\t')
            if c and a: anc[c].add(a)
    return anc


def aspect_map(anc):
    inv = {v: k for k, v in ASP_ROOT.items()}
    out = {}
    for t, A in anc.items():
        for root, asp in inv.items():
            if root in A: out[t] = asp; break
    return out


def build_edges(terms, tix, anc, device):
    """Descendant edge list within this aspect's term set.

    Returns (par, chi) long tensors: for every ordered pair where term i is an
    ancestor-or-self of term j (both in `terms`), an edge bucket i <- source j.
    MCM(h)[:,i] = max over j with edge (i,j) of h[:,j].
    Also returns (par_d, chi_d) = the same minus self-loops (for the softreg penalty).
    """
    par, chi = [], []
    for jt in terms:
        j = tix[jt]
        targets = {jt} | {a for a in anc.get(jt, ()) if a in tix}
        for t in targets:
            par.append(tix[t]); chi.append(j)
    par = torch.tensor(par, device=device); chi = torch.tensor(chi, device=device)
    nonself = par != chi
    return par, chi, par[nonself], chi[nonself]


def mcm(h, par, chi, n):
    """Coherent output: out[:,i] = max over descendants j (incl self) of h[:,j]."""
    B = h.shape[0]
    out = torch.zeros(B, n, device=h.device, dtype=h.dtype)
    idx = par.unsqueeze(0).expand(B, -1)
    out.scatter_reduce_(1, idx, h[:, chi], reduce='amax', include_self=True)
    return out


def mc_loss(logits, y, par, chi, n, posw, eps=1e-7):
    """C-HMCNN MCLoss: positives delegated to the most-confident TRUE descendant,
    negatives suppressed via the coherent (max-propagated) probability."""
    h = torch.sigmoid(logits)
    mcm_neg = mcm(h, par, chi, n)          # coherent prob (any descendant fires)
    mcm_pos = mcm(h * y, par, chi, n)      # max over descendants that are truly positive
    loss_pos = -(y * torch.log(mcm_pos.clamp_min(eps))) * posw
    loss_neg = -((1 - y) * torch.log((1 - mcm_neg).clamp_min(eps)))
    return (loss_pos + loss_neg).mean()


class MLP(nn.Module):
    def __init__(self, d, n, hidden=1024):
        super().__init__()
        self.net = nn.Sequential(nn.Linear(d, hidden), nn.GELU(), nn.Dropout(0.3),
                                 nn.Linear(hidden, n))
    def forward(self, x): return self.net(x)


def train_one(Xtr, Ytr, par, chi, par_d, chi_d, posw, args, dev):
    n_terms = Ytr.shape[1]
    model = MLP(Xtr.shape[1], n_terms).to(dev)
    opt = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
    bce = nn.BCEWithLogitsLoss(pos_weight=posw)
    N = Xtr.shape[0]
    for ep in range(args.epochs):
        model.train(); perm = torch.randperm(N, device=dev)
        for b in range(0, N, args.batch):
            bi = perm[b:b + args.batch]
            opt.zero_grad()
            logits = model(Xtr[bi]); yb = Ytr[bi]
            if args.loss == 'bce':
                loss = bce(logits, yb)
            elif args.loss == 'mcm':
                loss = mc_loss(logits, yb, par, chi, n_terms, posw)
            else:  # softreg
                h = torch.sigmoid(logits)
                viol = torch.relu(h[:, chi_d] - h[:, par_d]).mean()
                loss = bce(logits, yb) + args.lam * viol
            loss.backward(); opt.step()
    return model


def predict(model, X, dev, batch=4096):
    model.eval(); out = []
    with torch.no_grad():
        for b in range(0, X.shape[0], batch):
            out.append(torch.sigmoid(model(X[b:b + batch])).cpu().numpy())
    return np.concatenate(out, 0) if out else np.zeros((0, 0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--loss', choices=['bce', 'mcm', 'softreg'], default='mcm')
    ap.add_argument('--mode', choices=['test', 'oof'], default='test')
    ap.add_argument('--train-emb', help='required unless --load-model')
    ap.add_argument('--test-emb', help='required for --mode test')
    ap.add_argument('--labels', help='required unless --load-model')
    ap.add_argument('--dag', help='required unless --load-model')
    ap.add_argument('--out', required=True)
    ap.add_argument('--emb-model', help='embedding model tag stored in the checkpoint '
                    '(e.g. esm2_650m, prostt5); for the deploy apply path')
    ap.add_argument('--save-model', help='(--mode test) save the per-aspect MLP heads '
                    '+ vocab + emb metadata to this .pt for deployment')
    ap.add_argument('--load-model', help='apply-only: load a saved .pt and predict over '
                    '--apply-emb (skips training)')
    ap.add_argument('--apply-emb', help='embedding .npz (ids, emb) to score with --load-model')
    ap.add_argument('--min-pos', type=int, default=10)
    ap.add_argument('--epochs', type=int, default=30)
    ap.add_argument('--batch', type=int, default=512)
    ap.add_argument('--kfold', type=int, default=5)
    ap.add_argument('--lam', type=float, default=1.0, help='softreg penalty weight')
    ap.add_argument('--seed', type=int, default=0)
    args = ap.parse_args()
    dev = 'cuda' if torch.cuda.is_available() else 'cpu'
    rng = np.random.RandomState(args.seed)

    # ---- apply-only deploy path: load saved heads, score an embedding npz ----
    if args.load_model:
        ckpt = torch.load(args.load_model, map_location=dev, weights_only=False)
        ae = np.load(args.apply_emb or args.test_emb, allow_pickle=True)
        ids = list(ae['ids']); X = torch.tensor(ae['emb'].astype(np.float32), device=dev)
        log(f'apply {args.load_model} (emb_model={ckpt.get("emb_model")}) over {X.shape}')
        with open(args.out, 'w') as fout:
            for aspect in ('MF', 'BP', 'CC'):
                a = ckpt['aspects'].get(aspect)
                if not a: continue
                terms = a['terms']
                model = MLP(ckpt['dim'], len(terms), ckpt.get('hidden', 1024)).to(dev)
                model.load_state_dict(a['state_dict'])
                preds = predict(model, X, dev)
                for i, p in enumerate(ids):
                    row = preds[i]
                    for j in np.nonzero(row >= 0.01)[0]:
                        fout.write(f'{p}\t{terms[j]}\t{row[j]:.4f}\n')
        log(f'done -> {args.out}')
        return

    anc = load_dag(args.dag)
    aspect_of = aspect_map(anc)
    tr = np.load(args.train_emb, allow_pickle=True)
    tr_ids = list(tr['ids']); tr_emb = tr['emb'].astype(np.float32)
    if args.mode == 'test':
        te = np.load(args.test_emb, allow_pickle=True)
        te_ids = list(te['ids']); Xte = torch.tensor(te['emb'].astype(np.float32), device=dev)
    log(f'loss={args.loss} mode={args.mode} train {tr_emb.shape} dev={dev}')

    fout = open(args.out, 'w')
    saved = {}; emb_dim = tr_emb.shape[1]
    for aspect in ('MF', 'BP', 'CC'):
        ylab = load_labels(args.labels, aspect_of, aspect)
        freq = defaultdict(int)
        for gs in ylab.values():
            for g in gs: freq[g] += 1
        terms = sorted([t for t, c in freq.items() if c >= args.min_pos])
        if not terms: continue
        tix = {t: j for j, t in enumerate(terms)}
        par, chi, par_d, chi_d = build_edges(terms, tix, anc, dev)
        rows = np.array([i for i, p in enumerate(tr_ids) if p in ylab])
        Xall = tr_emb[rows]
        Y = np.zeros((len(rows), len(terms)), dtype=np.float32)
        for r, ri in enumerate(rows):
            for g in ylab[tr_ids[ri]]:
                if g in tix: Y[r, tix[g]] = 1.0
        pos = Y.sum(0); w = np.clip(np.log1p(pos.max() / np.maximum(pos, 1)), 0.25, 16)
        posw = torch.tensor(w, device=dev)
        log(f'{aspect}: {len(rows)} prot, {len(terms)} terms, {par.numel()} edges')

        if args.mode == 'test':
            Xtr = torch.tensor(Xall, device=dev); Ytr = torch.tensor(Y, device=dev)
            model = train_one(Xtr, Ytr, par, chi, par_d, chi_d, posw, args, dev)
            if args.save_model:
                saved[aspect] = {'terms': terms,
                                 'state_dict': {k: v.cpu() for k, v in model.state_dict().items()}}
            preds = predict(model, Xte, dev)
            for i, p in enumerate(te_ids):
                row = preds[i]
                for j in np.nonzero(row >= 0.01)[0]:
                    fout.write(f'{p}\t{terms[j]}\t{row[j]:.4f}\n')
            log(f'{aspect}: wrote test predictions')
        else:  # oof
            order = rng.permutation(len(rows)); folds = np.array_split(order, args.kfold)
            for k in range(args.kfold):
                te_idx = folds[k]
                tr_idx = np.concatenate([folds[j] for j in range(args.kfold) if j != k])
                Xtr = torch.tensor(Xall[tr_idx], device=dev); Ytr = torch.tensor(Y[tr_idx], device=dev)
                model = train_one(Xtr, Ytr, par, chi, par_d, chi_d, posw, args, dev)
                preds = predict(model, torch.tensor(Xall[te_idx], device=dev), dev)
                for r, ridx in enumerate(te_idx):
                    pid = tr_ids[rows[ridx]]; row = preds[r]
                    for j in np.nonzero(row >= 0.01)[0]:
                        fout.write(f'{pid}\t{terms[j]}\t{row[j]:.4f}\n')
                log(f'{aspect}: fold {k+1}/{args.kfold} done')
    fout.close()
    if args.mode == 'test' and args.save_model and saved:
        torch.save({'aspects': saved, 'dim': emb_dim, 'hidden': 1024,
                    'loss': args.loss,
                    'emb_model': args.emb_model or os.path.basename(args.train_emb or '')},
                   args.save_model)
        log(f'saved heads -> {args.save_model}  ({len(saved)} aspects, dim={emb_dim})')
    log(f'done -> {args.out}')


if __name__ == '__main__':
    main()
