#!/usr/bin/env python3
"""Merge per-genome M3 training TSVs, split train/valid by genome."""
import argparse
import os
import random
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--in-dir', required=True,
                    help='parent dir containing {genome}/training.tsv')
    ap.add_argument('--out-train', required=True)
    ap.add_argument('--out-valid', required=True)
    ap.add_argument('--valid-frac', type=float, default=0.2)
    ap.add_argument('--seed', type=int, default=42)
    args = ap.parse_args()

    genomes = sorted(os.listdir(args.in_dir))
    random.seed(args.seed)
    random.shuffle(genomes)
    n_valid = max(1, int(len(genomes) * args.valid_frac))
    valid_tags = set(genomes[:n_valid])
    print(f'[info] valid genomes ({n_valid}): {sorted(valid_tags)}', file=sys.stderr)
    print(f'[info] train genomes ({len(genomes) - n_valid}): '
          f'{sorted(set(genomes) - valid_tags)}', file=sys.stderr)

    header = None
    train_rows = 0
    valid_rows = 0
    with open(args.out_train, 'w') as ft, open(args.out_valid, 'w') as fv:
        for tag in genomes:
            path = os.path.join(args.in_dir, tag, 'training.tsv')
            if not os.path.exists(path):
                continue
            with open(path) as fin:
                h = fin.readline()
                if header is None:
                    header = h
                    ft.write(h)
                    fv.write(h)
                out = fv if tag in valid_tags else ft
                for line in fin:
                    out.write(line)
                    if tag in valid_tags:
                        valid_rows += 1
                    else:
                        train_rows += 1
    print(f'[info] train rows: {train_rows}, valid rows: {valid_rows}',
          file=sys.stderr)


if __name__ == '__main__':
    main()
