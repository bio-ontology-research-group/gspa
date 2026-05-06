#!/usr/bin/env python3
"""Train an SPL head on top of frozen ESM2 pooled embeddings for the
5,707-term DeepGO-Plus vocabulary with the GO true-path constraint.

Input artifacts (all produced by earlier steps):

- Pooled NPZ from ``embed_train_data.py`` — arrays ``pooled`` (N, D),
  ``labels`` (N, T), ``proteins`` (N,).
- Compiled SDD + vtree from ``compile_sdd.py``.

Architecture (phase 1 of the plan: train gate only, ESM2 frozen via the
precomputed pool). The gating MLP maps the pooled embedding to the
vector of circuit parameters; CircuitMPE computes the exact log-
likelihood of the label configuration under the constraint. Loss:
``-log p(y_true | x)`` via CircuitMPE.cross_entropy.

Training is fast because there is no encoder forward pass — a few
seconds per epoch on an RTX-4090 for the full 79k-protein training set.
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import sys
from pathlib import Path
from time import perf_counter

LOG = logging.getLogger("train_spl")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--train-pooled", type=Path, required=True)
    ap.add_argument("--val-frac", type=float, default=0.05,
                    help="Fraction of train to hold out for validation.")
    ap.add_argument("--sdd", type=Path, required=True)
    ap.add_argument("--vtree", type=Path, required=True)
    ap.add_argument("--hierarchy", type=Path, default=None,
                    help="Aspect hierarchy NPZ. If present + has "
                         "'full_vocab_indices', the pooled-NPZ labels are "
                         "sliced to those columns (per-aspect training).")
    ap.add_argument("--hmc-utils", type=Path, required=True,
                    help="Path to the vendored SPL hmc-utils directory.")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--gate-layers", default="1280,512,512",
                    help="Comma-separated hidden layer sizes of the gating MLP.")
    ap.add_argument("--num-reps", type=int, default=1)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--wd", type=float, default=1e-5)
    ap.add_argument("--batch-size", type=int, default=256)
    ap.add_argument("--epochs", type=int, default=200)
    ap.add_argument("--patience", type=int, default=20,
                    help="Early-stopping patience (epochs without val improvement).")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--eval-every", type=int, default=1,
                    help="Evaluate on validation every N epochs.")
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")

    # Import SPL vendored utilities
    hmc_utils_path = str(args.hmc_utils.resolve())
    if hmc_utils_path not in sys.path:
        sys.path.insert(0, hmc_utils_path)
    pypsdd_path = str((args.hmc_utils / "pypsdd").resolve())
    if pypsdd_path not in sys.path:
        sys.path.insert(0, pypsdd_path)

    import numpy as np
    import torch
    from torch.utils.data import DataLoader, TensorDataset
    from compute_mpe import CircuitMPE
    from GatingFunction import DenseGatingFunction

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    # --- Data
    LOG.info("loading pooled features: %s", args.train_pooled)
    data = np.load(args.train_pooled, allow_pickle=True)
    pooled = data["pooled"].astype(np.float32)
    labels = data["labels"].astype(np.float32)  # uint8 → float for BCE
    LOG.info("  pooled %s, labels %s", pooled.shape, labels.shape)

    # Optional per-aspect slicing
    if args.hierarchy is not None:
        h = np.load(args.hierarchy)
        if "full_vocab_indices" in h.files:
            idx = h["full_vocab_indices"].astype(np.int64)
            labels = labels[:, idx]
            LOG.info("  sliced labels → %s via full_vocab_indices", labels.shape)

    # Filter zero-row labels (rows with no in-vocab annotation — no gradient)
    keep_mask = labels.sum(axis=1) > 0
    pooled = pooled[keep_mask]
    labels = labels[keep_mask]
    # Also filter rows where pooled is all zero (never embedded)
    non_zero = np.abs(pooled).sum(axis=1) > 0
    pooled = pooled[non_zero]
    labels = labels[non_zero]
    LOG.info("  after filter: %d rows", len(pooled))

    # Train / val split
    n = len(pooled)
    val_n = max(1, int(n * args.val_frac))
    rng = np.random.RandomState(args.seed)
    perm = rng.permutation(n)
    val_idx = perm[:val_n]
    train_idx = perm[val_n:]
    X_train = torch.from_numpy(pooled[train_idx])
    Y_train = torch.from_numpy(labels[train_idx])
    X_val = torch.from_numpy(pooled[val_idx])
    Y_val = torch.from_numpy(labels[val_idx])
    LOG.info("train=%d val=%d", X_train.shape[0], X_val.shape[0])

    # --- Circuit
    LOG.info("loading circuit: sdd=%s vtree=%s", args.sdd, args.vtree)
    cmpe = CircuitMPE(str(args.vtree), str(args.sdd))
    LOG.info("  circuit loaded")

    # --- Gating
    gate_layers = [int(x) for x in args.gate_layers.split(",") if x]
    LOG.info("gate layers: %s", gate_layers)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    LOG.info("device: %s", device)

    gate = DenseGatingFunction(cmpe.beta,
                               gate_layers=gate_layers,
                               num_reps=args.num_reps).to(device)
    LOG.info("gate params: %d", sum(p.numel() for p in gate.parameters()))

    optimizer = torch.optim.Adam(gate.parameters(), lr=args.lr, weight_decay=args.wd)

    # --- Loaders
    def loader(X, Y, shuffle: bool) -> DataLoader:
        ds = TensorDataset(X, Y)
        return DataLoader(ds, batch_size=args.batch_size, shuffle=shuffle,
                          num_workers=0, pin_memory=False)

    train_loader = loader(X_train, Y_train, shuffle=True)
    val_loader = loader(X_val, Y_val, shuffle=False)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    best_val_nll = float("inf")
    best_epoch = -1
    since_best = 0
    history = []

    for epoch in range(args.epochs):
        gate.train()
        t0 = perf_counter()
        tot_loss = 0.0
        n_batches = 0
        for x, y in train_loader:
            x = x.to(device)
            y = y.to(device)
            thetas = gate(x)
            cmpe.set_params(thetas)
            loss = cmpe.cross_entropy(y, log_space=True).mean()
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            tot_loss += loss.item()
            n_batches += 1
        train_loss = tot_loss / max(1, n_batches)
        dt_train = perf_counter() - t0

        # Val
        val_loss = None
        if epoch % args.eval_every == 0:
            gate.eval()
            with torch.no_grad():
                tot = 0.0
                nb = 0
                for x, y in val_loader:
                    x = x.to(device)
                    y = y.to(device)
                    thetas = gate(x)
                    cmpe.set_params(thetas)
                    nll = cmpe.cross_entropy(y, log_space=True).mean()
                    tot += nll.item()
                    nb += 1
                val_loss = tot / max(1, nb)

        history.append({"epoch": epoch,
                        "train_loss": train_loss,
                        "val_loss": val_loss,
                        "dt_train_s": dt_train})

        LOG.info("epoch %d/%d  train_nll=%.4f  val_nll=%s  dt=%.1fs",
                 epoch + 1, args.epochs, train_loss,
                 f"{val_loss:.4f}" if val_loss is not None else "--",
                 dt_train)

        if val_loss is not None:
            if val_loss < best_val_nll - 1e-4:
                best_val_nll = val_loss
                best_epoch = epoch
                since_best = 0
                torch.save({
                    "gate_state_dict": gate.state_dict(),
                    "epoch": epoch,
                    "val_loss": val_loss,
                    "args": vars(args) | {k: str(v) for k, v in vars(args).items()
                                          if isinstance(v, Path)},
                    "gate_layers": gate_layers,
                }, args.out_dir / "best.pt")
            else:
                since_best += 1
                if since_best >= args.patience:
                    LOG.info("early stopping at epoch %d (best=%d, val=%.4f)",
                             epoch + 1, best_epoch + 1, best_val_nll)
                    break

        with (args.out_dir / "history.json").open("w") as fh:
            json.dump(history, fh, indent=2)

    LOG.info("best val NLL %.4f at epoch %d", best_val_nll, best_epoch + 1)


if __name__ == "__main__":
    main()
