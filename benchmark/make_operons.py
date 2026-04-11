#!/usr/bin/env python3
"""Extract operon-like clusters from a GFF: adjacent CDSs on the same strand with intergenic distance < 300bp."""
import re
import sys

gff = sys.argv[1]
out = sys.argv[2]
max_gap = int(sys.argv[3]) if len(sys.argv) > 3 else 300

cds = []
with open(gff) as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 9 or fields[2] != 'CDS':
            continue
        contig = fields[0]
        start = int(fields[3])
        end = int(fields[4])
        strand = fields[6]
        attrs = fields[8]
        m = re.search(r'Name=([^;]+)', attrs)
        if not m:
            m = re.search(r'protein_id=([^;]+)', attrs)
        if not m:
            continue
        pid = m.group(1)
        cds.append((contig, start, end, strand, pid))

cds.sort(key=lambda r: (r[0], r[1]))

operons = []
current = []
prev_end = None
prev_strand = None
prev_contig = None

for contig, start, end, strand, pid in cds:
    if (prev_contig == contig and
            prev_strand == strand and
            prev_end is not None and
            start - prev_end <= max_gap):
        current.append(pid)
    else:
        if len(current) >= 2:
            operons.append(current)
        current = [pid]
    prev_contig = contig
    prev_end = end
    prev_strand = strand

if len(current) >= 2:
    operons.append(current)

with open(out, 'w') as fh:
    for op in operons:
        fh.write('\t'.join(op) + '\n')

print(f"Wrote {len(operons)} operons (avg size {sum(len(o) for o in operons) / len(operons):.1f}) to {out}")
