#!/usr/bin/env python3
"""
Extract PGAP GO annotations for genomes using NCBI's gene2go file.

Strategy:
  1. From each genome's GFF, extract GeneID → protein_id mapping
     (from the Dbxref=GeneID:NNN attribute)
  2. Map protein_id → UniProt accession (via refseq_to_uniprot)
  3. Scan gene2go.gz for those GeneIDs
  4. Emit: uniprot_acc, aspect, go_term

This gives us PGAP's GO annotations for genomes whose RefSeq GFF
doesn't carry Ontology_term/go_function fields directly.
"""
import argparse
import gzip
import re
import sys
from collections import defaultdict

ASPECT_MAP = {'Function': 'MF', 'Process': 'BP', 'Component': 'CC'}


def extract_geneid_proteinid(gff_path):
    """Return {GeneID (str) -> protein_id} from CDS rows."""
    mapping = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] != 'CDS':
                continue
            attrs = fields[8]
            pid_m = re.search(r'protein_id=([^;]+)', attrs)
            gid_m = re.search(r'GeneID:(\d+)', attrs)
            if pid_m and gid_m:
                mapping[gid_m.group(1)] = pid_m.group(1)
    return mapping


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--gene2go', required=True)
    p.add_argument('--genomes', nargs='+', required=True,
                   help='tag:gff_path:refseq_to_uniprot_tsv')
    p.add_argument('--out-dir', required=True)
    args = p.parse_args()

    import os
    os.makedirs(args.out_dir, exist_ok=True)

    # Build per-genome: GeneID set + GeneID->UniProt
    genome_geneids = {}     # tag -> set of GeneID strings
    geneid_to_uniprot = {}  # global GeneID -> UniProt (via protein_id -> UniProt)
    all_geneids = set()

    for spec in args.genomes:
        parts = spec.split(':')
        tag, gff, refseq_map = parts[0], parts[1], parts[2]

        # Load RefSeq -> UniProt
        rs2uni = {}
        with open(refseq_map) as f:
            for line in f:
                cols = line.rstrip('\n').split('\t')
                if len(cols) >= 2:
                    rs2uni[cols[0]] = cols[1]

        # Extract GeneID -> protein_id from GFF
        gid_pid = extract_geneid_proteinid(gff)

        # Build GeneID -> UniProt
        gids = set()
        for gid, pid in gid_pid.items():
            uni = rs2uni.get(pid)
            if uni:
                geneid_to_uniprot[gid] = uni
                gids.add(gid)
        genome_geneids[tag] = gids
        all_geneids |= gids
        print(f'  {tag}: {len(gid_pid)} GeneIDs in GFF, {len(gids)} with UniProt mapping', flush=True)

    print(f'  union: {len(all_geneids)} GeneIDs', flush=True)

    # Scan gene2go
    genome_rows = {tag: [] for tag in genome_geneids}
    opener = gzip.open if args.gene2go.endswith('.gz') else open
    n = 0
    kept = 0
    with opener(args.gene2go, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            n += 1
            if n % 10_000_000 == 0:
                print(f'  {n:,} rows, {kept:,} kept', flush=True)
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 8:
                continue
            gid = fields[1]
            if gid not in all_geneids:
                continue
            go_id = fields[2]
            category = fields[7]
            aspect = ASPECT_MAP.get(category, category)
            uni = geneid_to_uniprot.get(gid)
            if not uni:
                continue
            for tag, gids in genome_geneids.items():
                if gid in gids:
                    genome_rows[tag].append((uni, aspect, go_id))
                    kept += 1

    print(f'  total: {n:,} rows, {kept:,} kept', flush=True)

    for tag, rows in genome_rows.items():
        out_path = os.path.join(args.out_dir, f'{tag}_pgap.tsv')
        with open(out_path, 'w') as out:
            out.write('accession\taspect\tgo_term\n')
            for row in rows:
                out.write('\t'.join(row) + '\n')
        unique_accs = len({r[0] for r in rows})
        print(f'  {tag}: {len(rows):,} annotations on {unique_accs:,} proteins -> {out_path}')


if __name__ == '__main__':
    main()
