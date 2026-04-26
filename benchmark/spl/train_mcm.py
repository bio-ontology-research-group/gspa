#!/usr/bin/env python3
"""Train a DeepGO-Plus-style head on frozen ESM2 pooled features with
C-HMCNN's Max-Constrained-Module (MCM) enforcing the GO true-path rule.

MCM is the post-hoc ancestor-propagation layer from Giunchiglia &
Lukasiewicz 2020. For each term i::

    out_i  =  max_{j ∈ {i} ∪ ancestors(i)}  raw_j

applied to sigmoid scores. This guarantees the output is ancestor-closed
(if a specific child scores high, all its ancestors inherit that score),
which is exactly the GO true-path rule. It is differentiable and
O(labels × ancestors_per_label) per forward pass.

Trade-off vs SPL (this is plan B for BP/MF whose SDDs blow up): MCM is
not a proper probabilistic model — the outputs are element-wise
sigmoid scores with max-constraint applied, not a normalised joint
distribution over valid label configurations. It satisfies D3
(consistency) and D6 (efficiency) of the SPL desiderata but drops D1
(probabilistic) and D2 (expressiveness beyond hierarchy). For the
true-path rule specifically, these limitations don't hurt benchmark
F-max.

Loss: BCE on the max-constrained output (C-HMCNN's MCLoss). Using
log-sigmoid + BCEWithLogitsLoss-style reparameterisation for numerical
stability. The effect is that training lifts ancestor scores to match
the strongest descendant, propagating gradients naturally.

Output: checkpoint with the gating MLP state_dict + the ancestor matrix
R (so predict_mcm.py can apply the constraint at inference).
"""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from time import perf_counter

LOG = logging.getLogger("train_mcm")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--train-pooled", type=Path, required=True)
    ap.add_argument("--hierarchy", type=Path, required=True,
                    help="NPZ from build_go_hierarchy.py (ancestors matrix).")
    ap.add_argument("--val-frac", type=float, default=0.05)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--hidden", type=int, default=2048,
                    help="Hidden dim of the MLP head.")
    ap.add_argument("--dropout", type=float, default=0.2)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--wd", type=float, default=1e-5)
    ap.add_argument("--batch-size", type=int, default=128)
    ap.add_argument("--epochs", type=int, default=200)
    ap.add_argument("--patience", type=int, default=15)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--pos-weight", type=float, default=None,
                    help="Scalar pos_weight for BCE. If unset, auto-compute "
                         "from label frequency: (n_neg / n_pos).")
    return ap.parse_args()


def get_constr_out(x, R):
    """MCM: for each term, take max over its ancestors' scores.
    x: (B, T) raw (pre-sigmoid) logits; output: (B, T) logits.
    R: (T, T) boolean — R[i, j] = True iff j ∈ {i} ∪ ancestors(i).
    """
    import torch
    # expand to (B, T, T), mask with R and take max over dim 2
    B, T = x.shape
    c_out = x.unsqueeze(1).expand(B, T, T)  # (B, T, T), each row = all scores
    R_b = R.unsqueeze(0).expand(B, T, T)    # (B, T, T), each row = ancestor mask for term i
    # Where R is False, mask to very negative so it doesn't contribute to max
    masked = torch.where(R_b, c_out, torch.full_like(c_out, -1e9))
    out, _ = masked.max(dim=2)
    return out


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")

    import numpy as np
    import torch
    import torch.nn as nn
    from torch.utils.data import DataLoader, TensorDataset

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    LOG.info("loading pooled features: %s", args.train_pooled)
    data = np.load(args.train_pooled, allow_pickle=True)
    pooled = data["pooled"].astype(np.float32)
    labels = data["labels"].astype(np.float32)
    LOG.info("  pooled %s labels %s", pooled.shape, labels.shape)

    keep = (labels.sum(axis=1) > 0) & (np.abs(pooled).sum(axis=1) > 0)
    pooled = pooled[keep]
    labels = labels[keep]
    N, D = pooled.shape
    T = labels.shape[1]
    LOG.info("  after filter: %d rows, %d features, %d labels", N, D, T)

    # Ancestor matrix R: R[i, j] = True iff j ∈ {i} ∪ ancestors(i)
    LOG.info("loading hierarchy: %s", args.hierarchy)
    h = np.load(args.hierarchy)
    anc = h["ancestors"].astype(bool)  # (T, T) — anc[i, j] = j is ancestor of i
    R = anc.copy()
    np.fill_diagonal(R, True)  # include self
    assert R.shape == (T, T), f"hierarchy size mismatch: {R.shape} vs {T}"
    LOG.info("  R density: %.3f%%", R.sum() / R.size * 100)

    # Train / val split
    rng = np.random.RandomState(args.seed)
    perm = rng.permutation(N)
    val_n = max(1, int(N * args.val_frac))
    val_idx = perm[:val_n]
    train_idx = perm[val_n:]
    X_train = torch.from_numpy(pooled[train_idx])
    Y_train = torch.from_numpy(labels[train_idx])
    X_val = torch.from_numpy(pooled[val_idx])
    Y_val = torch.from_numpy(labels[val_idx])
    LOG.info("train=%d val=%d", len(X_train), len(X_val))

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    LOG.info("device: %s", device)
    R_t = torch.from_numpy(R).to(device)

    # Simple MLP head: D → H → H → T
    head = nn.Sequential(
        nn.Linear(D, args.hidden),
        nn.ReLU(),
        nn.Dropout(args.dropout),
        nn.Linear(args.hidden, args.hidden),
        nn.ReLU(),
        nn.Dropout(args.dropout),
        nn.Linear(args.hidden, T),
    ).to(device)
    LOG.info("head params: %d", sum(p.numel() for p in head.parameters()))

    optimizer = torch.optim.Adam(head.parameters(), lr=args.lr, weight_decay=args.wd)

    # Positive-label upweight: sparse multi-label GO has ~50 positives
    # out of 5707 per protein; a uniform mean BCE is dominated by zeros.
    n_pos = float(labels.sum())
    n_neg = float(labels.size - n_pos)
    auto_pw = n_neg / max(1.0, n_pos)
    pw = args.pos_weight if args.pos_weight is not None else auto_pw
    LOG.info("pos_weight = %.2f  (auto %.2f, n_pos=%d n_neg=%d)",
             pw, auto_pw, int(n_pos), int(n_neg))
    pos_weight_t = torch.tensor(pw, device=device)

    train_loader = DataLoader(TensorDataset(X_train, Y_train),
                              batch_size=args.batch_size, shuffle=True)
    val_loader = DataLoader(TensorDataset(X_val, Y_val),
                            batch_size=args.batch_size, shuffle=False)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    best_val = float("inf")
    best_epoch = -1
    since_best = 0
    history = []

    for epoch in range(args.epochs):
        head.train()
        t0 = perf_counter()
        tot = 0.0
        nb = 0
        for x, y in train_loader:
            x = x.to(device)
            y = y.to(device)
            logits = head(x)
            constrained = get_constr_out(logits, R_t)
            # MCLoss of C-HMCNN (weighted):
            #   y=1: BCE on constrained outputs — wants σ(c_i) → 1
            #   y=0: BCE on raw logits — wants σ(x_i) → 0
            # Positives upweighted by pos_weight so the sparse-positive
            # signal isn't washed out by the 99%+ zeros.
            pos_term = - pos_weight_t * y * nn.functional.logsigmoid(constrained)
            neg_term = - (1 - y) * nn.functional.logsigmoid(-logits)
            loss = (pos_term + neg_term).mean()
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            tot += loss.item()
            nb += 1
        train_loss = tot / max(1, nb)
        dt_train = perf_counter() - t0

        head.eval()
        val_tot = 0.0
        vb = 0
        with torch.no_grad():
            for x, y in val_loader:
                x = x.to(device)
                y = y.to(device)
                logits = head(x)
                constrained = get_constr_out(logits, R_t)
                loss = (- y * nn.functional.logsigmoid(constrained)
                        - (1 - y) * nn.functional.logsigmoid(-logits)).mean()
                val_tot += loss.item()
                vb += 1
        val_loss = val_tot / max(1, vb)

        history.append({"epoch": epoch, "train_loss": train_loss,
                        "val_loss": val_loss, "dt_s": dt_train})
        LOG.info("epoch %d/%d  train=%.4f  val=%.4f  dt=%.1fs",
                 epoch + 1, args.epochs, train_loss, val_loss, dt_train)

        if val_loss < best_val - 1e-4:
            best_val = val_loss
            best_epoch = epoch
            since_best = 0
            torch.save({
                "head_state_dict": head.state_dict(),
                "ancestors": anc,
                "epoch": epoch,
                "val_loss": val_loss,
                "args": {k: (str(v) if isinstance(v, Path) else v)
                         for k, v in vars(args).items()},
            }, args.out_dir / "best.pt")
        else:
            since_best += 1
            if since_best >= args.patience:
                LOG.info("early stopping at epoch %d (best %d, val %.4f)",
                         epoch + 1, best_epoch + 1, best_val)
                break

        with (args.out_dir / "history.json").open("w") as fh:
            json.dump(history, fh, indent=2)

    LOG.info("best val %.4f at epoch %d", best_val, best_epoch + 1)


if __name__ == "__main__":
    main()
