#!/usr/bin/env python3
"""Build a comparison table of ablation F-max results."""
import argparse
import json
import os


def load(path):
    try:
        return json.load(open(path))
    except Exception:
        return None


def fmt(f, ci):
    if not ci:
        return f'{f:.3f}'
    return f'{f:.3f} [{ci[0]:.3f},{ci[1]:.3f}]'


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--dir', required=True)
    p.add_argument('--genomes', nargs='+', required=True)
    p.add_argument('--configs', nargs='+', default=['diamond', 'pfam', 'combined', 'priors'])
    args = p.parse_args()

    for truth in ('exp', 'all'):
        rows = [['genome'] + args.configs]
        for tag in args.genomes:
            row = [tag]
            for cfg in args.configs:
                d = load(os.path.join(args.dir, f'{tag}_{cfg}_fmax.json'))
                if not d:
                    row.append('-')
                    continue
                blk = d.get('by_truth', {}).get(truth)
                if not blk or not blk.get('results'):
                    row.append('-')
                    continue
                r = blk['results'][0]
                row.append(fmt(r['fmax_overall'], r.get('fmax_ci')))
            rows.append(row)

        label = 'Experimental-only truth' if truth == 'exp' else 'Full-GOA truth (all evidence)'
        print(f'\n### Ablation — F-max overall ({label})\n')
        widths = [max(len(str(r[i])) for r in rows) for i in range(len(rows[0]))]
        for i, r in enumerate(rows):
            print('  '.join(str(v).ljust(widths[j]) for j, v in enumerate(r)))
            if i == 0:
                print('-' * (sum(widths) + 2 * (len(widths) - 1)))


if __name__ == '__main__':
    main()
