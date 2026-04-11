#!/usr/bin/env python3
"""Filter a UniProt-style FASTA to exclude specific accessions."""
import argparse
import sys


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--fasta', required=True)
    p.add_argument('--exclude', required=True, help='File with UniProt accessions to drop')
    p.add_argument('--output', required=True)
    args = p.parse_args()

    with open(args.exclude) as f:
        exclude = {line.strip() for line in f if line.strip()}
    print(f'excluding {len(exclude):,} accessions', flush=True)

    kept = 0
    dropped = 0
    header = None
    seq_lines = []
    keep_current = True

    def flush(header, seq_lines, keep, out):
        if header is None:
            return 0, 0
        if keep:
            out.write(header)
            out.writelines(seq_lines)
            return 1, 0
        return 0, 1

    with open(args.fasta) as fh, open(args.output, 'w') as out:
        for line in fh:
            if line.startswith('>'):
                # flush prev
                k, d = flush(header, seq_lines, keep_current, out)
                kept += k
                dropped += d
                header = line
                seq_lines = []
                # parse accession from >sp|ACC|NAME or >tr|ACC|NAME
                parts = line.split('|', 2)
                acc = parts[1] if len(parts) >= 3 else line[1:].split()[0]
                keep_current = acc not in exclude
            else:
                seq_lines.append(line)
        # flush last
        k, d = flush(header, seq_lines, keep_current, out)
        kept += k
        dropped += d

    print(f'kept {kept:,}, dropped {dropped:,}')


if __name__ == '__main__':
    main()
