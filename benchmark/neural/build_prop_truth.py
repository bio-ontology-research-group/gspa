#!/usr/bin/env python3
"""Build a prop_annotations-based truth TSV from the DeepGO-Plus test parquet.

DeepGO-Plus / CAFA evaluation expects the truth set to be ancestor-propagated
(every annotation implies all its GO ancestors). The parquet ships both
``exp_annotations`` (experimental-evidence terms, non-propagated) and
``prop_annotations`` (propagated). Use ``prop_annotations`` for the F-max
comparison.
"""
from __future__ import annotations

import argparse
from pathlib import Path


def load_aspects(obo_path: Path) -> dict[str, str]:
    ns_map = {"molecular_function": "MF",
              "biological_process": "BP",
              "cellular_component": "CC"}
    out: dict[str, str] = {}
    cid = None
    with obo_path.open() as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith("[Term]"):
                cid = None
            elif line.startswith("id: GO:"):
                cid = line.split(" ", 1)[1].strip()
            elif line.startswith("namespace: "):
                ns = line.split(" ", 1)[1].strip()
                if cid and ns in ns_map:
                    out[cid] = ns_map[ns]
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--parquet", type=Path, required=True)
    ap.add_argument("--go-obo", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--column", default="prop_annotations",
                    choices=("prop_annotations", "exp_annotations", "annotations"))
    args = ap.parse_args()

    import pandas as pd
    df = pd.read_parquet(args.parquet)
    aspects = load_aspects(args.go_obo)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with args.out.open("w") as fh:
        fh.write("protein_id\taspect\tterm\n")
        for _, row in df.iterrows():
            pid = str(row["proteins"]).strip()
            terms = row.get(args.column)
            if terms is None:
                continue
            for t in terms:
                s = str(t)
                # annotations column has "GO:xxx|EVIDENCE" — strip evidence
                if "|" in s:
                    s = s.split("|", 1)[0]
                if s.startswith("GO:"):
                    asp = aspects.get(s, "")
                    fh.write(f"{pid}\t{asp}\t{s}\n")
                    n += 1
    print(f"wrote {n} rows to {args.out}")


if __name__ == "__main__":
    main()
