#!/usr/bin/env python3
"""
Scan GOA once, emit TWO truth sets per proteome:
  <tag>_truth_exp.tsv  -- experimental evidence codes only (EXP, IDA, IMP, ...)
  <tag>_truth_all.tsv  -- all evidence codes (including IEA) except NOT
"""
import argparse
import gzip
import os

EXPERIMENTAL = {
    'EXP', 'IDA', 'IMP', 'IPI', 'IGI', 'IEP',
    'HTP', 'HDA', 'HMP', 'HGI', 'HEP', 'TAS', 'IC',
}
ASPECT_MAP = {'F': 'MF', 'P': 'BP', 'C': 'CC'}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--goa', required=True)
    p.add_argument('--accessions', action='append', required=True,
                   help='Repeatable: tag:path/to/accessions.txt')
    p.add_argument('--out-dir', required=True)
    args = p.parse_args()

    tag_to_accs = {}
    all_accs = set()
    for spec in args.accessions:
        tag, path = spec.split(':', 1)
        with open(path) as f:
            accs = {line.strip() for line in f if line.strip()}
        tag_to_accs[tag] = accs
        all_accs |= accs
        print(f'  {tag}: {len(accs):,} accessions', flush=True)
    print(f'  union: {len(all_accs):,}', flush=True)

    os.makedirs(args.out_dir, exist_ok=True)
    tag_to_rows_exp = {t: [] for t in tag_to_accs}
    tag_to_rows_all = {t: [] for t in tag_to_accs}

    opener = gzip.open if args.goa.endswith('.gz') else open
    n = 0
    kept_exp = 0
    kept_all = 0
    print(f'Scanning {args.goa}...', flush=True)
    with opener(args.goa, 'rt') as f:
        for line in f:
            n += 1
            if n % 50_000_000 == 0:
                print(f'  {n:,} lines, exp={kept_exp:,}, all={kept_all:,}', flush=True)
            if line.startswith('!'):
                continue
            fields = line.split('\t', 10)
            if len(fields) < 9:
                continue
            acc = fields[1]
            if acc not in all_accs:
                continue
            qualifier = fields[3]
            if 'NOT' in qualifier.upper():
                continue
            evidence = fields[6]
            go_id = fields[4]
            aspect = ASPECT_MAP.get(fields[8], fields[8])
            row = (acc, aspect, go_id)
            for tag, accs in tag_to_accs.items():
                if acc in accs:
                    tag_to_rows_all[tag].append(row)
                    kept_all += 1
                    if evidence in EXPERIMENTAL:
                        tag_to_rows_exp[tag].append(row)
                        kept_exp += 1
    print(f'Total: {n:,} lines, exp={kept_exp:,}, all={kept_all:,}', flush=True)

    for tag in tag_to_accs:
        for suffix, rows in [('exp', tag_to_rows_exp[tag]),
                             ('all', tag_to_rows_all[tag])]:
            out_path = os.path.join(args.out_dir, f'{tag}_truth_{suffix}.tsv')
            with open(out_path, 'w') as out:
                out.write('accession\taspect\tgo_term\n')
                for row in rows:
                    out.write('\t'.join(row) + '\n')
            unique = len({r[0] for r in rows})
            print(f'  {tag} {suffix}: {len(rows):,} on {unique:,} -> {out_path}')


if __name__ == '__main__':
    main()
