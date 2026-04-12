#!/usr/bin/env python3
"""
Build per-genome RefSeq→UniProt mappings from the NCBI
gene_refseq_uniprotkb_collab.gz file.

For each genome, reads its RefSeq protein IDs (from the NCBI protein
FASTA or GFF) and its UniProt accessions (from the proteome FASTA),
then extracts matching rows from the collab file. Prefers 'identical'
over 'similar' matches.

This replaces the per-genome UniProt REST API xref_refseq query,
which returns empty results for some proteomes (e.g. Synechocystis).

Usage:
  build_refseq_uniprot_map.py \
    --collab gene_refseq_uniprotkb_collab.gz \
    --genomes tag1:refseq_ids.txt:uniprot_ids.txt \
              tag2:refseq_ids.txt:uniprot_ids.txt \
    --out-dir maps/
"""
import argparse
import gzip
import os
import sys
from collections import defaultdict


def load_ids(path):
    ids = set()
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                ids.add(line)
    return ids


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--collab', required=True,
                   help='gene_refseq_uniprotkb_collab.gz')
    p.add_argument('--genomes', nargs='+', required=True,
                   help='tag:refseq_ids_file:uniprot_ids_file (repeatable)')
    p.add_argument('--out-dir', required=True)
    args = p.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # Build per-genome lookup sets
    genome_refseq = {}   # tag -> set of RefSeq IDs
    genome_uniprot = {}  # tag -> set of UniProt IDs
    all_refseq = set()
    all_uniprot = set()
    for spec in args.genomes:
        parts = spec.split(':')
        tag, refseq_file, uniprot_file = parts[0], parts[1], parts[2]
        rs = load_ids(refseq_file)
        uni = load_ids(uniprot_file)
        genome_refseq[tag] = rs
        genome_uniprot[tag] = uni
        all_refseq |= rs
        all_uniprot |= uni
        print(f'  {tag}: {len(rs):,} RefSeq, {len(uni):,} UniProt', flush=True)
    print(f'  union: {len(all_refseq):,} RefSeq, {len(all_uniprot):,} UniProt', flush=True)

    # Scan collab file: keep rows where BOTH refseq ∈ all_refseq AND uniprot ∈ all_uniprot
    # Store: (refseq, uniprot, method) tuples per genome
    genome_maps = {tag: {} for tag in genome_refseq}  # tag -> {refseq -> (uniprot, method)}

    opener = gzip.open if args.collab.endswith('.gz') else open
    n = 0
    kept = 0
    print(f'Scanning {args.collab}...', flush=True)
    with opener(args.collab, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            n += 1
            if n % 20_000_000 == 0:
                print(f'  {n:,} rows, {kept:,} kept', flush=True)
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 5:
                continue
            refseq = fields[0]
            uniprot = fields[1]
            method = fields[4]

            if refseq not in all_refseq:
                continue
            if uniprot not in all_uniprot:
                continue

            for tag in genome_refseq:
                if refseq in genome_refseq[tag] and uniprot in genome_uniprot[tag]:
                    existing = genome_maps[tag].get(refseq)
                    # Prefer 'identical' over 'similar'
                    if existing is None or (existing[1] != 'identical' and method == 'identical'):
                        genome_maps[tag][refseq] = (uniprot, method)
                        kept += 1

    print(f'Total: {n:,} rows, {kept:,} kept', flush=True)

    for tag in genome_refseq:
        out_path = os.path.join(args.out_dir, f'{tag}.refseq_to_uniprot.tsv')
        with open(out_path, 'w') as out:
            for refseq, (uniprot, method) in sorted(genome_maps[tag].items()):
                out.write(f'{refseq}\t{uniprot}\n')
        n_identical = sum(1 for _, m in genome_maps[tag].values() if m == 'identical')
        print(f'  {tag}: {len(genome_maps[tag]):,} mappings ({n_identical:,} identical) -> {out_path}')


if __name__ == '__main__':
    main()
