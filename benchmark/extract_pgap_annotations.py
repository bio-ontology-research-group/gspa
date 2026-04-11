#!/usr/bin/env python3
"""
Extract PGAP's GO annotations from an NCBI genomic.gff and emit a TSV
keyed on UniProt accession (for comparison against experimental GOA).

PGAP GFFs expose three kinds of functional evidence:
  - Ontology_term=GO:0006260,GO:0003887,...  (direct GO)
  - Dbxref=UniProtKB/Swiss-Prot:P12345      (for matching to GOA ground truth)
  - product="thr operon leader peptide"     (text; not parsed here)

The script also accepts an external mapping from RefSeq protein ID to
UniProt accession for cases where the PGAP GFF doesn't provide a
UniProtKB Dbxref (e.g., H. pylori, M. jannaschii). The mapping is
derived from the UniProt proteome FASTA headers, which carry both
(RefSeq Xref via the UniProt accession lookup).

Output TSV columns:
  uniprot_acc   aspect (MF/BP/CC)   go_term   source_protein_id
"""
import argparse
import re
import sys
from collections import defaultdict


def parse_gff_attrs(s):
    out = {}
    for kv in s.split(';'):
        kv = kv.strip()
        if not kv or '=' not in kv:
            continue
        k, v = kv.split('=', 1)
        out[k] = v
    return out


def uniprot_from_dbxref(dbxref):
    for part in dbxref.split(','):
        part = part.strip()
        if part.startswith('UniProtKB/Swiss-Prot:'):
            return part.split(':', 1)[1]
        if part.startswith('UniProtKB/TrEMBL:'):
            return part.split(':', 1)[1]
        if part.startswith('UniProtKB:'):
            return part.split(':', 1)[1]
    return None


def load_uniprot_fasta_xref(uniprot_fasta):
    """
    UniProt FASTA headers look like:
        >sp|P12345|NAME_SPECIES Description OS=... OX=...
    Each header carries the UniProt accession. Cross-reference to RefSeq
    protein_id comes via `xref` files from UniProt; here we use the
    minimal approach of loading the proteome's UniProt accessions and
    returning a dict {uniprot_acc: None} that the caller augments with a
    mapping file.
    """
    accs = set()
    with open(uniprot_fasta) as fh:
        for line in fh:
            if line.startswith('>'):
                parts = line.split('|', 2)
                if len(parts) >= 3:
                    accs.add(parts[1])
    return accs


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--gff', required=True)
    p.add_argument('--uniprot-proteome-fasta', default=None,
                   help='Used to sanity-filter by UniProt accessions in the target proteome')
    p.add_argument('--refseq-to-uniprot', default=None,
                   help='TSV: refseq_protein_id\\tuniprot_acc')
    p.add_argument('--go-owl-aspect-map', default=None,
                   help='Optional TSV go_term\\taspect to fill missing aspects')
    p.add_argument('--out', required=True)
    args = p.parse_args()

    # 1. Build optional RefSeq -> UniProt mapping
    refseq_to_uniprot = {}
    if args.refseq_to_uniprot:
        with open(args.refseq_to_uniprot) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 2:
                    refseq_to_uniprot[parts[0]] = parts[1]

    # 2. Optional proteome filter (keep only rows whose UniProt acc is in this set)
    target_uniprot = None
    if args.uniprot_proteome_fasta:
        target_uniprot = load_uniprot_fasta_xref(args.uniprot_proteome_fasta)
        print(f'target proteome: {len(target_uniprot):,} UniProt accessions', flush=True)

    # 3. Optional GO -> aspect map (for filling when GFF doesn't tell us)
    go_aspect = {}
    if args.go_owl_aspect_map:
        with open(args.go_owl_aspect_map) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 2:
                    go_aspect[parts[0]] = parts[1]

    rows = []
    missing_uniprot = 0
    total_cds = 0
    cds_with_go = 0
    with open(args.gff) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] != 'CDS':
                continue
            total_cds += 1
            attrs = parse_gff_attrs(fields[8])
            # Resolve UniProt acc
            uniprot = None
            dbxref = attrs.get('Dbxref', '')
            if dbxref:
                uniprot = uniprot_from_dbxref(dbxref)
            if not uniprot:
                protein_id = attrs.get('protein_id') or attrs.get('Name')
                if protein_id and protein_id in refseq_to_uniprot:
                    uniprot = refseq_to_uniprot[protein_id]
            if not uniprot:
                missing_uniprot += 1
                continue
            if target_uniprot and uniprot not in target_uniprot:
                continue

            ontology = attrs.get('Ontology_term', '')
            terms = {}
            if ontology:
                for t in ontology.split(','):
                    t = t.strip()
                    if t.startswith('GO:'):
                        terms[t] = None  # aspect filled below
                cds_with_go += 1

            # Try to derive per-term aspect from go_function / go_process / go_component
            go_function = attrs.get('go_function', '')
            go_process = attrs.get('go_process', '')
            go_component = attrs.get('go_component', '')
            for block, aspect in [(go_function, 'MF'),
                                   (go_process, 'BP'),
                                   (go_component, 'CC')]:
                # Format: "label|GO:NNNNNNN||evidence"
                for chunk in block.split(','):
                    m = re.search(r'GO:(\d+)', chunk)
                    if m:
                        gid = f'GO:{m.group(1)}'
                        terms[gid] = aspect

            # Fill missing aspects from the map
            for gid in list(terms.keys()):
                if terms[gid] is None:
                    terms[gid] = go_aspect.get(gid)

            for gid, aspect in terms.items():
                if aspect is None:
                    continue  # drop unknown-aspect entries; can't evaluate per-aspect F-max
                rows.append((uniprot, aspect, gid))

    with open(args.out, 'w') as out:
        out.write('accession\taspect\tgo_term\n')
        for row in rows:
            out.write('\t'.join(row) + '\n')

    unique_accs = len({r[0] for r in rows})
    print(f'total CDS: {total_cds:,}')
    print(f'CDS with GO: {cds_with_go:,}')
    print(f'CDS missing UniProt cross-ref: {missing_uniprot:,}')
    print(f'annotations written: {len(rows):,} on {unique_accs:,} proteins -> {args.out}')


if __name__ == '__main__':
    main()
