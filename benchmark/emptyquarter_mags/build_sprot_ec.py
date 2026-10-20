#!/usr/bin/env python3
"""Swiss-Prot accession -> EC number table from uniprot_sprot.dat(.gz).

Reads the DE lines of every entry ("DE   ...; EC=1.1.1.1 {ECO:...};" or "DE            EC=1.1.1.1;") and
writes one row per (primary accession, complete EC number); partial numbers (with "-") are skipped.

  build_sprot_ec.py --dat uniprot_sprot.dat.gz --out sprot_ec.tsv
"""
import argparse
import gzip
import re

EC_RE = re.compile(r"EC=(\d+\.\d+\.\d+\.\d+)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dat", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    op = gzip.open if a.dat.endswith(".gz") else open
    n_entries = n_rows = 0
    acc, ecs = None, set()
    with op(a.dat, "rt", errors="replace") as f, open(a.out, "w") as out:
        for line in f:
            if line.startswith("AC   ") and acc is None:
                acc = line[5:].split(";")[0].strip()
            elif line.startswith("DE   "):
                ecs.update(EC_RE.findall(line))
            elif line.startswith("//"):
                n_entries += 1
                if acc and ecs:
                    for ec in sorted(ecs):
                        out.write(f"{acc}\tEC:{ec}\n")
                        n_rows += 1
                acc, ecs = None, set()
    print(f"{n_entries} entries, {n_rows} (accession, EC) rows -> {a.out}")


if __name__ == "__main__":
    main()
