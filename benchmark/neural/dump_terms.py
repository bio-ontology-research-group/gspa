#!/usr/bin/env python3
"""Dump a DeepGO-Plus-style terms pickle/parquet to a plain text file.

The GSPA ESM2-DeepGO-Plus predictor's sidecar expects a text file with
one GO term per line (index = FC output column). deepgo-nesy and its
cluster-side twin ``gapfix`` store the same ordered vocabulary as a
1-column pandas DataFrame in ``terms.pkl`` / ``terms.parquet`` with
column name ``terms`` (or ``functions`` for the older ``mf.pkl`` subset).
This tiny helper materializes the list.

Usage::

    dump_terms.py --in /data/hohndor/gapfix/data/deepgoplus-real/data/terms.pkl \\
                  --out go_terms_5707.txt
"""
from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="src", type=Path, required=True,
                    help="terms.pkl or terms.parquet (single-column DataFrame)")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--column", default=None,
                    help="Column to read (default: first column of the DataFrame)")
    args = ap.parse_args()

    import pandas as pd
    if args.src.suffix == ".parquet":
        df = pd.read_parquet(args.src)
    else:
        df = pd.read_pickle(args.src)
    if not hasattr(df, "columns"):
        raise SystemExit(f"{args.src} did not load as a DataFrame")
    col = args.column or df.columns[0]
    if col not in df.columns:
        raise SystemExit(f"column {col!r} not in {df.columns.tolist()}")
    terms = [str(t) for t in df[col].tolist()]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(terms) + "\n")
    print(f"wrote {len(terms)} terms to {args.out}")


if __name__ == "__main__":
    main()
