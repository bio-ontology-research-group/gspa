#!/usr/bin/env python3
"""Dump the deepgo-nesy/DeepGO-Plus test_data parquet as FASTA + truth TSVs.

Outputs three files alongside {out_prefix}:

- ``{out_prefix}.faa`` — FASTA; headers use the ``proteins`` column as the ID
  (matches how the sidecar emits predictions).
- ``{out_prefix}_truth_exp.tsv`` — experimental-evidence-only truth,
  compatible with ``benchmark_pgap_v2.py``'s ``load_truth``:
  ``protein_id\\taspect\\tterm``. Aspect is filled in from ``--go-obo``
  if supplied; otherwise left blank (benchmark_pgap_v2 treats blank
  truth aspects as matching predictions in the same '' aspect bucket,
  so if you populate aspects here you must also populate them on the
  prediction side).
- ``{out_prefix}_truth_all.tsv`` — everything in ``annotations`` minus the
  evidence-code suffix.

Usage::

    dump_test_fasta.py \\
        --in /data/hohndor/gapfix/data/deepgoplus-real/data/test_data.parquet \\
        --out-prefix /data/hohndor/gspa-neural/work/dgp_test \\
        --go-obo /data/hohndor/gapfix/data/deepgoplus-real/data/go.obo
"""
from __future__ import annotations

import argparse
from pathlib import Path


def load_go_aspects(obo_path: Path) -> dict[str, str]:
    aspects = {"molecular_function": "MF", "biological_process": "BP",
               "cellular_component": "CC"}
    mapping: dict[str, str] = {}
    current_id = None
    with obo_path.open() as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith("[Term]"):
                current_id = None
            elif line.startswith("id: GO:"):
                current_id = line.split(" ", 1)[1].strip()
            elif line.startswith("namespace: "):
                ns = line.split(" ", 1)[1].strip()
                if current_id and ns in aspects:
                    mapping[current_id] = aspects[ns]
    return mapping


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="src", type=Path, required=True)
    ap.add_argument("--out-prefix", type=Path, required=True)
    ap.add_argument("--go-obo", type=Path, default=None,
                    help="Optional go.obo used to fill in aspect column of the truth TSVs")
    args = ap.parse_args()

    import pandas as pd
    df = pd.read_parquet(args.src)
    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)

    aspects = load_go_aspects(args.go_obo) if args.go_obo else {}

    faa = args.out_prefix.with_suffix(".faa")
    exp_tsv = Path(f"{args.out_prefix}_truth_exp.tsv")
    all_tsv = Path(f"{args.out_prefix}_truth_all.tsv")

    with faa.open("w") as fh_faa, exp_tsv.open("w") as fh_exp, all_tsv.open("w") as fh_all:
        fh_exp.write("protein_id\taspect\tterm\n")
        fh_all.write("protein_id\taspect\tterm\n")
        for _, row in df.iterrows():
            pid = str(row["proteins"]).strip()
            seq = str(row["sequences"]).strip()
            if not pid or not seq:
                continue
            fh_faa.write(f">{pid}\n{seq}\n")
            exp_terms = row.get("exp_annotations")
            if exp_terms is not None:
                for term in exp_terms:
                    t = str(term)
                    if t.startswith("GO:"):
                        asp = aspects.get(t, "")
                        fh_exp.write(f"{pid}\t{asp}\t{t}\n")
            all_terms = row.get("annotations")
            if all_terms is not None:
                for entry in all_terms:
                    t = str(entry).split("|", 1)[0]
                    if t.startswith("GO:"):
                        asp = aspects.get(t, "")
                        fh_all.write(f"{pid}\t{asp}\t{t}\n")

    print(f"wrote {faa}")
    print(f"wrote {exp_tsv}")
    print(f"wrote {all_tsv}")


if __name__ == "__main__":
    main()
