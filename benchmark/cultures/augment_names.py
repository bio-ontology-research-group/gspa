"""Add reaction_name and ec_name columns to prediction TSVs.

Reaction name lookup order:
  1. Per-culture gapsmith Reactions.tbl (has MetaCyc + SEED IDs)
  2. SEED reactions.tsv (for SEED IDs not in culture's tbl)
EC name: extracted from ec2go.txt.
"""
import collections
import csv
import re
import sys
from pathlib import Path


def load_ec_names(path):
    out = {}
    with open(path) as f:
        for line in f:
            if line.startswith('!') or not line.strip():
                continue
            # Format:  EC:1.1.1.1 > GO:alcohol ... activity ; GO:0004022
            m = re.match(r'^EC:(\S+)\s+>\s+GO:(.+?)\s*;\s*GO:\d+',
                          line.strip())
            if not m:
                continue
            ec = m.group(1)
            name = m.group(2).strip()
            # Prefer the entry with a specific EC (keep shortest name
            # since partial ECs like 1.-.-.- have generic names)
            if ec not in out or len(name) < len(out[ec]):
                out[ec] = name
    return out


def load_seed_reaction_names(path):
    out = {}
    with open(path) as f:
        f.readline()
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 3:
                out[parts[0]] = parts[2]  # id -> name
    return out


def load_culture_rxn_names(reactions_tbl):
    out = {}
    with open(reactions_tbl) as f:
        rdr = csv.DictReader(f, delimiter='\t')
        for row in rdr:
            rxn = row['rxn']
            name = (row.get('name') or '').strip()
            if rxn and name and rxn not in out:
                out[rxn] = name
    return out


def augment(input_path, output_path, ec_names, seed_names,
             culture_names, ec_col='gap_ec', rxn_col='gap_rxn'):
    with open(input_path) as fin, open(output_path, 'w') as fout:
        h = fin.readline().rstrip('\n').split('\t')
        new_h = h + ['reaction_name', 'ec_name']
        fout.write('\t'.join(new_h) + '\n')
        ec_i = h.index(ec_col)
        rxn_i = h.index(rxn_col)
        for line in fin:
            parts = line.rstrip('\n').split('\t')
            rxn = parts[rxn_i]
            ec = parts[ec_i]
            rname = culture_names.get(rxn) or seed_names.get(rxn) or ''
            ename = ec_names.get(ec, '')
            parts.extend([rname, ename])
            fout.write('\t'.join(parts) + '\n')


def main():
    ec_path = '/data/hohndor/gspa/reference/ec2go.txt'
    seed_path = '/data/hohndor/gspa/bin/gapsmith/data_merged/seed_reactions.tsv'
    ec_names = load_ec_names(ec_path)
    seed_names = load_seed_reaction_names(seed_path)
    print(f'[info] {len(ec_names)} EC names, '
          f'{len(seed_names)} SEED reaction names', file=sys.stderr)

    # Per culture: augment dark_matter_predictions.tsv
    for tag in ['MR59-1', 'MR60-1', 'C-1.1', 'C-1.3']:
        tbl = f'/data/hohndor/gspa/proteomes/cultures/{tag}/gapsmith/{tag}-all-Reactions.tbl'
        culture_names = load_culture_rxn_names(tbl)
        print(f'[{tag}] {len(culture_names)} culture-specific '
              f'reaction names', file=sys.stderr)

        src = f'/data/hohndor/gspa/proteomes/cultures/{tag}/dark_matter_predictions.tsv'
        dst = f'/data/hohndor/gspa/proteomes/cultures/{tag}/dark_matter_predictions_named.tsv'
        augment(src, dst, ec_names, seed_names, culture_names)
        print(f'  wrote {dst}', file=sys.stderr)

    # Validation shortlist — union of culture-name tables
    all_culture_names = {}
    for tag in ['MR59-1', 'MR60-1', 'C-1.1', 'C-1.3']:
        tbl = f'/data/hohndor/gspa/proteomes/cultures/{tag}/gapsmith/{tag}-all-Reactions.tbl'
        all_culture_names.update(load_culture_rxn_names(tbl))
    for val_src in ['validation_candidates_pident.tsv',
                     'validation_candidates.tsv']:
        src = f'/data/hohndor/gspa/proteomes/cultures/{val_src}'
        dst = f'/data/hohndor/gspa/proteomes/cultures/' \
              f'{Path(val_src).stem}_named.tsv'
        augment(src, dst, ec_names, seed_names, all_culture_names)
        print(f'[validation] wrote {dst}', file=sys.stderr)


if __name__ == '__main__':
    main()
