#!/usr/bin/env python3
"""Relabel panel genome protein FASTAs to UniProt accessions.

Given {tag}_protein.faa (headers use RefSeq IDs) and {tag}_map.tsv
(RefSeq -> UniProt), emit {tag}.faa with UniProt accessions as the sole
header identifier. Proteins without a mapping are dropped (reported).
"""
import argparse
import os
import sys


def load_map(path):
    m = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 2:
                m[parts[0]] = parts[1]
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--in-fasta', required=True)
    ap.add_argument('--map', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    m = load_map(args.map)
    written = dropped = 0
    with open(args.in_fasta) as fin, open(args.out, 'w') as fout:
        keep = False
        for line in fin:
            if line.startswith('>'):
                header = line[1:].split()[0]
                up = m.get(header)
                if up:
                    fout.write(f'>{up}\n')
                    keep = True
                    written += 1
                else:
                    keep = False
                    dropped += 1
            elif keep:
                fout.write(line)
    print(f'[info] {written} written, {dropped} dropped (unmapped) -> {args.out}',
          file=sys.stderr)


if __name__ == '__main__':
    main()
