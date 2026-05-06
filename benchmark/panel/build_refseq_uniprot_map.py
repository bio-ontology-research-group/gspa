#!/usr/bin/env python3
"""Build RefSeq_protein_id -> UniProt_accession maps for panel genomes.

Uses the gene_refseq_uniprotkb_collab.gz dump (already available at
/data/hohndor/gspa/reference/). For each genome, parses the GFF to get
the set of RefSeq protein IDs (NP_/WP_/YP_...), filters the collab file
to produce {tag}_map.tsv with format: refseq_id <TAB> uniprot_acc.

If multiple UniProt accessions map to one RefSeq ID we keep the first
(typically Swiss-Prot beats TrEMBL by collab ordering).
"""
import argparse
import gzip
import os
import re
import sys


def refseq_ids_from_gff(gff_path):
    ids = set()
    open_fn = gzip.open if gff_path.endswith('.gz') else open
    with open_fn(gff_path, 'rt') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) < 9 or parts[2] != 'CDS':
                continue
            m = re.search(r'protein_id=([^;]+)', parts[8])
            if m:
                ids.add(m.group(1))
    return ids


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--manifest', required=True)
    ap.add_argument('--genomes-dir', required=True,
                    help='dir containing {tag}_genomic.gff')
    ap.add_argument('--collab', required=True,
                    help='gene_refseq_uniprotkb_collab.gz')
    ap.add_argument('--out-dir', required=True)
    args = ap.parse_args()

    tags = []
    per_tag_ids = {}
    with open(args.manifest) as f:
        header = f.readline()
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 1 or not parts[0]:
                continue
            tag = parts[0]
            gff = os.path.join(args.genomes_dir, f'{tag}_genomic.gff')
            if not os.path.exists(gff):
                print(f'[warn] {gff} missing; skipping {tag}', file=sys.stderr)
                continue
            ids = refseq_ids_from_gff(gff)
            print(f'[info] {tag}: {len(ids)} RefSeq IDs', file=sys.stderr)
            tags.append(tag)
            per_tag_ids[tag] = ids

    # Build one master "wanted" set for efficient collab scan
    wanted = set()
    for ids in per_tag_ids.values():
        wanted.update(ids)
    print(f'[info] total wanted RefSeq IDs: {len(wanted)}', file=sys.stderr)

    # Scan collab in one pass
    refseq2uniprot = {}
    with gzip.open(args.collab, 'rt') as f:
        header = f.readline()
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            rs, up = parts[0], parts[1]
            if rs in wanted and rs not in refseq2uniprot:
                refseq2uniprot[rs] = up
    print(f'[info] mapped {len(refseq2uniprot)}/{len(wanted)} RefSeq IDs',
          file=sys.stderr)

    os.makedirs(args.out_dir, exist_ok=True)
    for tag in tags:
        path = os.path.join(args.out_dir, f'{tag}_map.tsv')
        with open(path, 'w') as f:
            for rs in sorted(per_tag_ids[tag]):
                up = refseq2uniprot.get(rs)
                if up:
                    f.write(f'{rs}\t{up}\n')
        print(f'[info] wrote {path}', file=sys.stderr)


if __name__ == '__main__':
    main()
