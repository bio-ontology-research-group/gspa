#!/usr/bin/env python3
"""Extract a FASTA for a set of accessions from a UniProt SwissProt flat file.

Used to build the pre-t0 training FASTA for the DG++-Light 1D-CNN component from
`uniprot_sprot.dat.gz` (the same flat file `build_text_string_index.py` parses),
restricted to the accessions in `train_terms.tsv`.

Usage:
  python extract_sprot_fasta.py uniprot_sprot.dat.gz accessions.txt out.fasta
where accessions.txt is one accession per line (or a TSV whose first column is the
accession; a header is skipped if its first field is 'EntryID'/'accession').
"""
from __future__ import annotations

import gzip
import sys


def opn(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def load_accessions(path):
    acc = set()
    with open(path) as fh:
        for i, line in enumerate(fh):
            f = line.rstrip('\n').split('\t')[0].strip()
            if i == 0 and f.lower() in ('entryid', 'accession', 'protein'):
                continue
            if f:
                acc.add(f)
    return acc


def main():
    dat, acc_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
    want = load_accessions(acc_path)
    print(f'want {len(want):,} accessions', file=sys.stderr)

    n = 0
    with opn(dat) as fh, open(out_path, 'w') as out:
        accs = []
        seq_lines = []
        in_seq = False
        for line in fh:
            if line.startswith('AC '):
                for a in line[5:].replace(' ', '').split(';'):
                    if a:
                        accs.append(a)
            elif line.startswith('SQ '):
                in_seq = True
                seq_lines = []
            elif line.startswith('//'):
                primary = accs[0] if accs else None
                hit = next((a for a in accs if a in want), None)
                if hit:
                    seq = ''.join(seq_lines).replace(' ', '')
                    out.write(f'>{hit}\n{seq}\n')
                    n += 1
                    if n % 10000 == 0:
                        print(f'  wrote {n:,}', file=sys.stderr)
                accs = []
                seq_lines = []
                in_seq = False
            elif in_seq:
                seq_lines.append(line.strip())
    print(f'wrote {n:,} sequences -> {out_path}', file=sys.stderr)


if __name__ == '__main__':
    main()
