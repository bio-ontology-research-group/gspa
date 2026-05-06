#!/usr/bin/env python3
"""Extract (UniProt_acc, EC_number) rows from uniprot_sprot.dat.gz.

For benchmark EC truth: we need EC numbers per SwissProt protein. Each
SwissProt flat entry has ``DE   EC=X.X.X.X;`` lines (possibly multiple)
inside ``DE   RecName`` / ``AltName`` / ``SubName`` blocks. Also ``DE
Flags: Fragment;`` etc. We harvest every ``EC=`` occurrence and emit one
row per (acc, ec).

Usage::
    extract_ec_truth.py --dat uniprot_sprot.dat.gz --accessions accs.txt --out ec.tsv
"""
from __future__ import annotations

import argparse
import gzip
import re
import sys
from pathlib import Path


EC_RE = re.compile(r"EC=([0-9n.\-]+)")


def load_accs(path: Path) -> set[str]:
    return {line.strip() for line in path.read_text().splitlines() if line.strip()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dat", required=True, help="uniprot_sprot.dat.gz")
    ap.add_argument("--accessions", required=True,
                    help="One UniProt accession per line; restricts output.")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    targets = load_accs(Path(args.accessions))
    print(f"targets: {len(targets):,}", file=sys.stderr)

    opener = gzip.open if args.dat.endswith(".gz") else open
    accs = []       # list of all AC tokens in current entry
    ecs: set[str] = set()
    rows: list[tuple[str, str]] = []
    n_entries = 0

    def flush():
        nonlocal accs, ecs
        matched = [a for a in accs if a in targets]
        for acc in matched:
            for ec in ecs:
                if ec.endswith("."):
                    continue  # partial / malformed
                rows.append((acc, ec))
        accs = []
        ecs = set()

    with opener(args.dat, "rt") as fh:
        for line in fh:
            if line.startswith("//"):
                flush()
                n_entries += 1
                if n_entries % 100_000 == 0:
                    print(f"  {n_entries:,} entries, rows={len(rows):,}",
                          file=sys.stderr)
                continue
            if line.startswith("AC   "):
                for tok in line[5:].rstrip(";\n").split(";"):
                    tok = tok.strip()
                    if tok:
                        accs.append(tok)
            elif line.startswith("DE   ") and "EC=" in line:
                for m in EC_RE.finditer(line):
                    ec = m.group(1).rstrip(".;")
                    if ec:
                        ecs.add(ec)
        flush()

    # Dedup and write
    out_rows = sorted(set(rows))
    with open(args.out, "w") as fh:
        fh.write("accession\taspect\tfunction_id\n")
        for acc, ec in out_rows:
            fh.write(f"{acc}\t\tEC:{ec}\n")
    print(f"wrote {args.out}: {len(out_rows):,} rows from {n_entries:,} entries",
          file=sys.stderr)


if __name__ == "__main__":
    main()
