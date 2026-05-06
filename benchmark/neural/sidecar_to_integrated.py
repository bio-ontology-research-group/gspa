#!/usr/bin/env python3
"""Convert a sidecar TSV (protein_id, term, score, annotation_type) into the
integrated-TSV shape that ``benchmark_pgap_v2.py`` consumes. Optionally
fills in GO aspect from a supplied ``go.obo``.

Output columns match ``load_gspa_integrated``'s required fields::

    protein_id\\ttype\\tfunction_id\\tgo_aspect\\tposterior_prob

All other columns of the real integrated TSV are omitted (the loader
picks fields by index from the header).
"""
from __future__ import annotations

import argparse
from pathlib import Path


def load_go_aspects(obo_path: Path) -> dict[str, str]:
    """Parse go.obo to extract {GO_ID: MF|BP|CC}. Minimal parser."""
    aspects = {"molecular_function": "MF", "biological_process": "BP",
               "cellular_component": "CC"}
    mapping: dict[str, str] = {}
    current_id: str | None = None
    current_ns: str | None = None
    with obo_path.open() as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith("[Term]"):
                current_id = None
                current_ns = None
            elif line.startswith("id: GO:"):
                current_id = line.split(" ", 1)[1].strip()
            elif line.startswith("namespace: "):
                current_ns = line.split(" ", 1)[1].strip()
                if current_id and current_ns in aspects:
                    mapping[current_id] = aspects[current_ns]
    return mapping


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="src", type=Path, required=True,
                    help="Sidecar TSV (protein_id, term, score, annotation_type)")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--go-obo", type=Path, default=None,
                    help="Optional go.obo for aspect lookup")
    args = ap.parse_args()

    aspects = load_go_aspects(args.go_obo) if args.go_obo else {}

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.src.open() as fin, args.out.open("w") as fout:
        fout.write("protein_id\ttype\tfunction_id\tgo_aspect\tposterior_prob\n")
        header = fin.readline()  # skip
        n = 0
        for line in fin:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            pid, term, score, ann = parts[0], parts[1], parts[2], parts[3]
            aspect = aspects.get(term, "") if ann == "GO" else ""
            fout.write(f"{pid}\t{ann}\t{term}\t{aspect}\t{score}\n")
            n += 1
    print(f"wrote {n} rows to {args.out}")


if __name__ == "__main__":
    main()
