#!/usr/bin/env python3
"""Pre-compute frozen ESM2 pooled embeddings + binary GO-label targets
for a DeepGO-Plus train/test parquet, so the SPL training loop can run
as a pure MLP + circuit pass with no encoder overhead.

Reads a parquet with columns ``proteins``, ``sequences``, and either
``prop_annotations`` (for train) or ``exp_annotations`` (for test). Emits
one NPZ with::

    pooled      (N, D)   float32 — L2-**not**-normalized mean-pooled
                                  last-layer representations
    labels      (N, T)   uint8   — binary label tensor
    proteins    (N,)     str     — protein IDs, same order as pooled/labels

Usage::

    embed_train_data.py \\
        --parquet /data/hohndor/gapfix/data/deepgoplus-real/data/train_data.parquet \\
        --terms   /data/hohndor/gspa-neural/work/go_terms_5707.txt \\
        --label-col prop_annotations \\
        --model esm2_t33_650M_UR50D \\
        --batch-size 8 \\
        --out pooled_train.npz
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

LOG = logging.getLogger("embed_train_data")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--parquet", type=Path, required=True)
    ap.add_argument("--terms", type=Path, required=True,
                    help="Vocabulary file, one GO ID per line.")
    ap.add_argument("--label-col", default="prop_annotations",
                    choices=("prop_annotations", "exp_annotations", "annotations"))
    ap.add_argument("--model", default="esm2_t33_650M_UR50D")
    ap.add_argument("--batch-size", type=int, default=8)
    ap.add_argument("--max-seq-len", type=int, default=1022)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--max-rows", type=int, default=None,
                    help="For smoke tests: limit to the first N rows.")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")

    import numpy as np
    import pandas as pd
    import torch
    import esm

    # Vocabulary
    vocab = [ln.strip() for ln in args.terms.read_text().splitlines() if ln.strip()]
    vocab_idx = {t: i for i, t in enumerate(vocab)}
    T = len(vocab)
    LOG.info("vocabulary: %d terms", T)

    # Parquet
    LOG.info("loading parquet: %s", args.parquet)
    df = pd.read_parquet(args.parquet)
    if args.max_rows:
        df = df.head(args.max_rows)
    LOG.info("  %d rows", len(df))

    # Labels
    N = len(df)
    labels = np.zeros((N, T), dtype=np.uint8)
    proteins = np.empty(N, dtype=object)
    kept: list[int] = []
    for i, (_, row) in enumerate(df.iterrows()):
        pid = str(row["proteins"]).strip()
        proteins[i] = pid
        terms = row.get(args.label_col)
        if terms is None or len(terms) == 0:
            # No label => cannot learn from this row; skip
            continue
        hit_any = False
        for t in terms:
            s = str(t).split("|", 1)[0]
            j = vocab_idx.get(s)
            if j is not None:
                labels[i, j] = 1
                hit_any = True
        if hit_any:
            kept.append(i)
    LOG.info("rows with ≥1 in-vocab label: %d / %d", len(kept), N)

    # Encoder
    LOG.info("loading %s", args.model)
    loader = getattr(esm.pretrained, args.model)
    model, alphabet = loader()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = model.to(device).eval()
    batch_converter = alphabet.get_batch_converter()
    pad_idx = alphabet.padding_idx
    num_layers = model.num_layers
    embed_dim = model.embed_dim
    LOG.info("  embed_dim=%d, device=%s", embed_dim, device)

    pooled = np.zeros((N, embed_dim), dtype=np.float32)

    def pool_batch(batch: list[tuple[str, str]]):
        trimmed = [(pid, seq[: args.max_seq_len]) for pid, seq in batch]
        _labels, _strs, tokens = batch_converter(trimmed)
        tokens = tokens.to(device)
        with torch.no_grad():
            out = model(tokens, repr_layers=[num_layers], return_contacts=False)
        reps = out["representations"][num_layers]  # (B, L, D)
        mask = (tokens != pad_idx)
        mask[:, 0] = False
        last_true = mask.long().sum(dim=1) - 1
        for i in range(mask.size(0)):
            mask[i, last_true[i].item()] = False
        w = mask.float().unsqueeze(-1)
        embs = (reps * w).sum(dim=1) / w.sum(dim=1).clamp(min=1.0)
        return embs.cpu().numpy().astype(np.float32)

    # Iterate
    batch: list[tuple[int, str, str]] = []
    for i, (_, row) in enumerate(df.iterrows()):
        pid = str(row["proteins"]).strip()
        seq = str(row["sequences"]).strip()
        if not seq:
            continue
        batch.append((i, pid, seq))
        if len(batch) >= args.batch_size:
            idxs = [b[0] for b in batch]
            embs = pool_batch([(b[1], b[2]) for b in batch])
            pooled[idxs] = embs
            batch = []
        if (i + 1) % (args.batch_size * 50) == 0:
            LOG.info("  embedded %d / %d", i + 1, N)
    if batch:
        idxs = [b[0] for b in batch]
        embs = pool_batch([(b[1], b[2]) for b in batch])
        pooled[idxs] = embs

    # Save
    args.out.parent.mkdir(parents=True, exist_ok=True)
    LOG.info("saving %s", args.out)
    np.savez_compressed(
        args.out,
        pooled=pooled,
        labels=labels,
        proteins=proteins.astype(str),
    )
    LOG.info("done")


if __name__ == "__main__":
    main()
