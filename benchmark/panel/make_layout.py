#!/usr/bin/env python3
"""Emit genome-layout TSV: protein_id(uniprot) contig start end strand.

Parses GFF CDS rows, maps RefSeq protein_id via map TSV to UniProt acc.
Coalesces multi-segment CDS into the outermost span.

Output columns:
  protein_id  contig  start  end  strand
"""
import argparse
import re


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--gff', required=True)
    ap.add_argument('--map', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    m = {}
    with open(args.map) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 2:
                m[parts[0]] = parts[1]

    coords = {}
    with open(args.gff) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9 or parts[2] != 'CDS':
                continue
            contig, start, end, strand = parts[0], int(parts[3]), int(parts[4]), parts[6]
            mm = re.search(r'protein_id=([^;]+)', parts[8])
            if not mm:
                continue
            rs = mm.group(1)
            up = m.get(rs)
            if not up:
                continue
            if up in coords:
                c = coords[up]
                if c['contig'] == contig:
                    c['start'] = min(c['start'], start)
                    c['end'] = max(c['end'], end)
            else:
                coords[up] = {'contig': contig, 'start': start,
                              'end': end, 'strand': strand}

    with open(args.out, 'w') as f:
        f.write('protein_id\tcontig\tstart\tend\tstrand\n')
        for up, c in sorted(coords.items(), key=lambda x: (x[1]['contig'], x[1]['start'])):
            f.write(f'{up}\t{c["contig"]}\t{c["start"]}\t{c["end"]}\t{c["strand"]}\n')
    print(f'[info] wrote {len(coords)} rows -> {args.out}')


if __name__ == '__main__':
    main()
