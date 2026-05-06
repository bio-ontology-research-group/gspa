#!/usr/bin/env python3
"""Run an MCM-trained head on a pooled-feature NPZ and emit a predictions
TSV in the gspa benchmark sidecar format::

    protein_id\\tterm\\tscore\\tannotation_type

Applies the same max-over-ancestors constraint at inference that training
used, so the emitted scores are already ancestor-consistent (true-path
rule satisfied by construction).
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

LOG = logging.getLogger("predict_mcm")


def get_constr_out(x, R):
    import torch
    B, T = x.shape
    c_out = x.unsqueeze(1).expand(B, T, T)
    R_b = R.unsqueeze(0).expand(B, T, T)
    masked = torch.where(R_b, c_out, torch.full_like(c_out, -1e9))
    out, _ = masked.max(dim=2)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--checkpoint", type=Path, required=True)
    ap.add_argument("--pooled", type=Path, required=True,
                    help="Pooled NPZ with 'pooled' and 'proteins' arrays.")
    ap.add_argument("--terms", type=Path, required=True)
    ap.add_argument("--out-tsv", type=Path, required=True)
    ap.add_argument("--batch-size", type=int, default=32)
    ap.add_argument("--min-score", type=float, default=0.01)
    ap.add_argument("--hidden", type=int, default=2048,
                    help="Must match training config.")
    ap.add_argument("--dropout", type=float, default=0.0,
                    help="Set to 0 at inference.")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")

    import numpy as np
    import torch
    import torch.nn as nn

    vocab = [ln.strip() for ln in args.terms.read_text().splitlines() if ln.strip()]
    T = len(vocab)
    LOG.info("vocabulary: %d terms", T)

    data = np.load(args.pooled, allow_pickle=True)
    pooled = data["pooled"].astype(np.float32)
    proteins = [str(p) for p in data["proteins"]]
    N, D = pooled.shape
    LOG.info("pooled: N=%d D=%d", N, D)

    LOG.info("loading checkpoint: %s", args.checkpoint)
    ck = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    anc = ck["ancestors"].astype(bool)
    R = anc.copy()
    np.fill_diagonal(R, True)
    head = nn.Sequential(
        nn.Linear(D, args.hidden), nn.ReLU(), nn.Dropout(args.dropout),
        nn.Linear(args.hidden, args.hidden), nn.ReLU(), nn.Dropout(args.dropout),
        nn.Linear(args.hidden, T),
    )
    head.load_state_dict(ck["head_state_dict"])
    head.eval()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    head = head.to(device)
    R_t = torch.from_numpy(R).to(device)

    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_tsv.open("w") as fh:
        fh.write("protein_id\tterm\tscore\tannotation_type\n")
        with torch.no_grad():
            for start in range(0, N, args.batch_size):
                end = min(N, start + args.batch_size)
                x = torch.from_numpy(pooled[start:end]).to(device)
                logits = head(x)
                constrained = get_constr_out(logits, R_t)
                probs = torch.sigmoid(constrained).cpu().numpy()
                for i, pid in enumerate(proteins[start:end]):
                    p = probs[i]
                    keep = p >= args.min_score
                    for j in np.where(keep)[0]:
                        fh.write(f"{pid}\t{vocab[j]}\t{p[j]:.4f}\tGO\n")
                if (start // args.batch_size) % 20 == 0:
                    LOG.info("  %d / %d", end, N)
    LOG.info("wrote %s", args.out_tsv)


if __name__ == "__main__":
    main()
