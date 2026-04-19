#!/usr/bin/env python3
"""Convert MMseqs2 easy-cluster output (rep <TAB> member) to a
protein_id <TAB> cluster_id TSV suitable for --orthogroups on
`gspa-cli integrate`."""
import sys
with open(sys.argv[1]) as fh, open(sys.argv[2], "w") as out:
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 2:
            out.write(f"{parts[1]}\t{parts[0]}\n")
