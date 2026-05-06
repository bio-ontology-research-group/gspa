#!/usr/bin/env python3
"""Fuse multiple predictor TSVs into an ensemble TSV for one genome.

Reads all the given predictor TSVs (normalized sidecar format:
  protein_id<TAB>term<TAB>score<TAB>annotation_type),
and emits a combined TSV where each (protein, term) pair gets a single
fused score. Supported fusion modes:

  max   : max over predictors (default)
  mean  : mean over predictors (NaN-skipped; 0 if absent)
  rank  : per-predictor rank-normalized (1/rank) then max

Usage::
    build_ensemble_preds.py \
        --tag mg1655 \
        --pred /path/to/esm2-deepgoplus.tsv \
        --pred /path/to/proteinfer.tsv \
        --pred /path/to/clean.tsv \
        --out /path/to/mg1655.ensemble.tsv \
        --mode max
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path


HEADER = "protein_id\tterm\tscore\tannotation_type"


def load_one(path: Path):
    """Return {(pid, term, ann): score}; if duplicates within a predictor,
    keep the maximum."""
    out: dict[tuple[str, str, str], float] = {}
    with path.open() as fh:
        hdr = next(fh, "").rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(hdr)}
        required = ("protein_id", "term", "score", "annotation_type")
        for c in required:
            if c not in idx:
                raise SystemExit(f"{path} missing column {c}")
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            try:
                score = float(fields[idx["score"]])
            except ValueError:
                continue
            key = (fields[idx["protein_id"]],
                   fields[idx["term"]],
                   fields[idx["annotation_type"]])
            if key not in out or out[key] < score:
                out[key] = score
    return out


def fuse_max(loaded: list[dict]):
    merged: dict[tuple[str, str, str], float] = {}
    for d in loaded:
        for k, v in d.items():
            if k not in merged or merged[k] < v:
                merged[k] = v
    return merged


def fuse_mean(loaded: list[dict]):
    sums: dict[tuple[str, str, str], float] = defaultdict(float)
    cnts: dict[tuple[str, str, str], int] = defaultdict(int)
    keys: set = set()
    for d in loaded:
        keys |= d.keys()
    # Include a predictor's 0-score for keys it doesn't emit. This makes
    # mean fusion penalize predictors that drop below min_score.
    n = len(loaded)
    for k in keys:
        total = 0.0
        for d in loaded:
            total += d.get(k, 0.0)
        sums[k] = total / n
    return dict(sums)


def fuse_rank(loaded: list[dict]):
    """Per-predictor: convert score to reciprocal rank within that predictor,
    then take max across predictors. Score = 1 / (rank_in_predictor)."""
    per_pred = []
    for d in loaded:
        # rank by descending score; ties keep insertion order
        sorted_items = sorted(d.items(), key=lambda kv: -kv[1])
        rr: dict[tuple[str, str, str], float] = {}
        for i, (k, _) in enumerate(sorted_items):
            rr[k] = 1.0 / (1 + i)
        per_pred.append(rr)
    return fuse_max(per_pred)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred", action="append", required=True,
                    help="repeatable: path to one predictor TSV")
    ap.add_argument("--out", required=True)
    ap.add_argument("--mode", choices=["max", "mean", "rank"], default="max")
    ap.add_argument("--min-score", type=float, default=0.0,
                    help="drop output rows with fused score below this")
    args = ap.parse_args()

    loaded = [load_one(Path(p)) for p in args.pred]
    if args.mode == "max":
        merged = fuse_max(loaded)
    elif args.mode == "mean":
        merged = fuse_mean(loaded)
    else:
        merged = fuse_rank(loaded)

    with open(args.out, "w") as fh:
        fh.write(HEADER + "\n")
        for (pid, term, ann), score in sorted(merged.items()):
            if score < args.min_score:
                continue
            fh.write(f"{pid}\t{term}\t{score:.4f}\t{ann}\n")
    print(f"wrote {args.out}: {len(merged)} rows")


if __name__ == "__main__":
    main()
