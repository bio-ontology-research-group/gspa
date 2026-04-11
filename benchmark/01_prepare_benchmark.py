#!/usr/bin/env python3
"""
GSPA Benchmark Step 1: Prepare benchmark dataset from Swiss-Prot + GOA.

Creates a benchmark set of proteins with experimentally validated GO annotations.
Splits into train/test sets for weight optimization.

Usage:
    python3 01_prepare_benchmark.py \
        --sprot uniprot_sprot.fasta.gz \
        --goa goa_uniprot_all.gaf.gz \
        --outdir benchmark_data \
        --n-test 5000
"""

import argparse
import gzip
import os
import random
from collections import defaultdict

EXPERIMENTAL_EVIDENCE = {'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP',
                         'HTP', 'HDA', 'HMP', 'HGI', 'HEP', 'TAS'}

def parse_goa(goa_file):
    """Parse GOA GAF file, returning accession -> {aspect: set(GO terms)} for experimental evidence."""
    annotations = defaultdict(lambda: defaultdict(set))
    opener = gzip.open if goa_file.endswith('.gz') else open
    with opener(goa_file, 'rt') as f:
        for line in f:
            if line.startswith('!'): continue
            fields = line.strip().split('\t')
            if len(fields) < 15: continue
            db = fields[0]
            accession = fields[1]
            qualifier = fields[3]
            go_id = fields[4]
            evidence = fields[6]
            aspect = fields[8]  # F=MF, P=BP, C=CC

            if 'NOT' in qualifier.upper(): continue
            if evidence not in EXPERIMENTAL_EVIDENCE: continue
            if db not in ('UniProtKB', 'Swiss-Prot'): continue

            aspect_name = {'F': 'MF', 'P': 'BP', 'C': 'CC'}.get(aspect, aspect)
            annotations[accession][aspect_name].add(go_id)

    return annotations

def parse_sprot_fasta(fasta_file):
    """Parse Swiss-Prot FASTA, return dict of accession -> sequence."""
    sequences = {}
    opener = gzip.open if fasta_file.endswith('.gz') else open
    current_acc = None
    current_seq = []
    with opener(fasta_file, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                if current_acc:
                    sequences[current_acc] = ''.join(current_seq)
                # >sp|P12345|NAME_ORGANISM Description
                parts = line.split('|')
                current_acc = parts[1] if len(parts) >= 2 else line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.strip())
    if current_acc:
        sequences[current_acc] = ''.join(current_seq)
    return sequences

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sprot', required=True)
    parser.add_argument('--goa', required=True)
    parser.add_argument('--outdir', default='benchmark_data')
    parser.add_argument('--n-test', type=int, default=5000,
                        help='Number of test proteins')
    parser.add_argument('--min-annotations', type=int, default=3,
                        help='Minimum GO annotations per protein')
    parser.add_argument('--max-length', type=int, default=1000,
                        help='Max protein length (for speed)')
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    random.seed(args.seed)

    # Parse GOA
    print("Parsing GOA annotations...")
    annotations = parse_goa(args.goa)
    print(f"  {len(annotations)} proteins with experimental GO annotations")

    # Parse Swiss-Prot sequences
    print("Parsing Swiss-Prot FASTA...")
    sequences = parse_sprot_fasta(args.sprot)
    print(f"  {len(sequences)} sequences")

    # Filter: proteins with both sequence and experimental annotations
    eligible = []
    for acc in annotations:
        if acc not in sequences: continue
        total_anns = sum(len(v) for v in annotations[acc].values())
        if total_anns < args.min_annotations: continue
        if len(sequences[acc]) > args.max_length: continue
        eligible.append(acc)

    print(f"  {len(eligible)} proteins with >= {args.min_annotations} experimental annotations and length <= {args.max_length}")

    # Shuffle and split
    random.shuffle(eligible)
    test_accs = set(eligible[:args.n_test])
    train_accs = set(eligible[args.n_test:args.n_test*3])  # 2x for training

    print(f"  Test set: {len(test_accs)} proteins")
    print(f"  Train set: {len(train_accs)} proteins")

    # Write test FASTA
    with open(os.path.join(args.outdir, 'test_proteins.fasta'), 'w') as f:
        for acc in sorted(test_accs):
            f.write(f'>{acc}\n{sequences[acc]}\n')

    # Write train FASTA
    with open(os.path.join(args.outdir, 'train_proteins.fasta'), 'w') as f:
        for acc in sorted(train_accs):
            f.write(f'>{acc}\n{sequences[acc]}\n')

    # Write ground truth annotations
    for split_name, split_accs in [('test', test_accs), ('train', train_accs)]:
        with open(os.path.join(args.outdir, f'{split_name}_go_annotations.tsv'), 'w') as f:
            f.write('accession\taspect\tgo_term\n')
            for acc in sorted(split_accs):
                for aspect, terms in annotations[acc].items():
                    for term in sorted(terms):
                        f.write(f'{acc}\t{aspect}\t{term}\n')

    # Write Swiss-Prot FASTA excluding test proteins (for building reference DB)
    print("Writing reference database (Swiss-Prot minus test proteins)...")
    with open(os.path.join(args.outdir, 'reference_db.fasta'), 'w') as f:
        for acc in sorted(sequences.keys()):
            if acc in test_accs: continue  # Exclude test proteins from reference
            f.write(f'>{acc}\n{sequences[acc]}\n')

    # Stats
    for aspect in ['MF', 'BP', 'CC']:
        test_with = sum(1 for a in test_accs if aspect in annotations[a])
        test_terms = sum(len(annotations[a].get(aspect, set())) for a in test_accs)
        print(f"  {aspect}: {test_with} test proteins, {test_terms} annotations")

    print("Done. Files written to:", args.outdir)

if __name__ == '__main__':
    main()
