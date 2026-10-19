#!/usr/bin/env python3
"""
Convert an annotation TSV + proteome FASTA into a synthetic GFF3 + GAF pair
that `gspa evaluate` can ingest.

Inputs:
  --fasta     UniProt-style FASTA (bare accessions after header cleanup)
  --tsv       annotation TSV. Two supported schemas:
                 * GSPA integrated (header: ...protein_id...type...function_id...go_aspect...posterior_prob...)
                 * PGAP simple (header: accession, aspect, go_term)
  --threshold Score threshold (only used for GSPA schema; default 0.5)
  --gff-out   Output synthetic GFF3
  --gaf-out   Output GAF

Synthetic GFF3 puts one CDS per protein on a single "synthetic1" contig
with dummy positions. Feature ID = protein accession (matches the GAF keys).
"""

import argparse
import sys
from pathlib import Path

ASPECT_LETTER = {'MF': 'F', 'BP': 'P', 'CC': 'C'}


def read_fasta_ids(path):
    ids = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                tok = line[1:].split()[0]
                if '|' in tok:
                    parts = tok.split('|')
                    if len(parts) >= 2:
                        tok = parts[1]
                ids.append(tok)
    return ids


def load_gspa(tsv, threshold):
    """Return dict protein_id -> set of (aspect, go_term)."""
    out = {}
    with open(tsv) as fh:
        headers = next(fh, '').rstrip('\n').split('\t')
        idx_pid = headers.index('protein_id')
        idx_type = headers.index('type')
        idx_func = headers.index('function_id')
        idx_aspect = headers.index('go_aspect')
        idx_prob = headers.index('posterior_prob')
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) <= idx_prob or f[idx_type] != 'GO':
                continue
            try:
                s = float(f[idx_prob])
            except ValueError:
                continue
            if s < threshold:
                continue
            pid = f[idx_pid]
            aspect = f[idx_aspect] or 'BP'
            out.setdefault(pid, set()).add((aspect, f[idx_func]))
    return out


def load_pgap(tsv):
    out = {}
    with open(tsv) as fh:
        next(fh, None)
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            out.setdefault(parts[0], set()).add((parts[1], parts[2]))
    return out


def detect_schema(tsv):
    with open(tsv) as fh:
        hdr = next(fh, '').rstrip('\n').split('\t')
    if 'posterior_prob' in hdr and 'function_id' in hdr:
        return 'gspa'
    return 'pgap'


def write_gff(fasta_ids, annots, gff_path, acc_mapper=None):
    """Write one 300bp CDS per FASTA protein on a synthetic contig.
    If acc_mapper is provided, GFF feature IDs are the MAPPED (UniProt) ids
    so they match the GAF keys from pgap_annotations.tsv.
    """
    pos = 1
    with open(gff_path, 'w') as out:
        out.write('##gff-version 3\n')
        out.write('##sequence-region synthetic1 1 {}\n'.format(300 * len(fasta_ids)))
        for pid in fasta_ids:
            used_id = acc_mapper[pid] if acc_mapper and pid in acc_mapper else pid
            end = pos + 299
            out.write(f'synthetic1\tsynth\tgene\t{pos}\t{end}\t.\t+\t0\tID=gene-{used_id}\n')
            out.write(f'synthetic1\tsynth\tCDS\t{pos}\t{end}\t.\t+\t0\tID={used_id};Parent=gene-{used_id}\n')
            pos = end + 1


def write_gaf(fasta_ids, annots, gaf_path, acc_mapper=None):
    """GAF 2.2 minimal."""
    with open(gaf_path, 'w') as out:
        out.write('!gaf-version: 2.2\n')
        for pid in fasta_ids:
            used_id = acc_mapper[pid] if acc_mapper and pid in acc_mapper else pid
            terms = annots.get(used_id) or set()
            for aspect, go_term in terms:
                aletter = ASPECT_LETTER.get(aspect, 'P')
                # db, db_obj_id, db_obj_sym, qualifier, go_id, ref, evidence, with, aspect, name, syn, type, taxon, date, src, ext, isoform
                out.write('\t'.join([
                    'SYNTH', used_id, used_id, '', go_term, 'GO_REF:synth',
                    'IEA', '', aletter, '', '', 'protein', 'taxon:0',
                    '20260411', 'GSPA_BENCH', '', ''
                ]) + '\n')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--fasta', required=True)
    p.add_argument('--tsv', required=True)
    p.add_argument('--schema', choices=['auto', 'gspa', 'pgap'], default='auto')
    p.add_argument('--threshold', type=float, default=0.5)
    p.add_argument('--gff-out', required=True)
    p.add_argument('--gaf-out', required=True)
    p.add_argument('--acc-mapper', default=None,
                   help='TSV orig_acc\\tmapped_acc, used to normalize FASTA ids to match TSV keys')
    args = p.parse_args()

    ids = read_fasta_ids(args.fasta)
    print(f'FASTA ids: {len(ids):,}')

    acc_mapper = None
    if args.acc_mapper:
        acc_mapper = {}
        with open(args.acc_mapper) as fh:
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 2:
                    acc_mapper[parts[0]] = parts[1]

    schema = args.schema if args.schema != 'auto' else detect_schema(args.tsv)
    if schema == 'gspa':
        annots = load_gspa(args.tsv, args.threshold)
    else:
        annots = load_pgap(args.tsv)
    print(f'schema={schema}, proteins_with_annotations={len(annots):,}')

    write_gff(ids, annots, args.gff_out, acc_mapper=acc_mapper)
    write_gaf(ids, annots, args.gaf_out, acc_mapper=acc_mapper)
    print(f'wrote {args.gff_out} and {args.gaf_out}')


if __name__ == '__main__':
    main()
