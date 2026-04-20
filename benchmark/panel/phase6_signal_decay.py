#!/usr/bin/env python3
"""Phase 6b analysis — signal-decay for shortlist entries.

For each (orthogroup, reaction_id) pair that appears in the Phase 5
master shortlist, compute log_lr at panel sizes N=29, N=60, N=97 using
the three catalogs built by run_signal_decay.sh. Report:

  - monotonic_increase: log_lr strictly grows with N (suggests
    sample-size artifact)
  - plateau: log_lr ~stable from N=60 onward (suggests real signal)
  - noisy: neither trend (fluctuates)

The test is meaningful only for pairs that actually appear in all three
catalogs with both n_sig_total > 0 and n_base_with > 0.
"""
import argparse
import collections
import csv
import math
import sys
from pathlib import Path


def load_catalog_subset(path, keys_of_interest, eps=1e-3):
    """{(og, rxn): log_lr} for the (og, rxn) pairs requested.
    Streams the catalog to avoid RAM blow-up."""
    wanted = set(keys_of_interest)
    out = {}
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if parts[0] == 'orthogroup_id' or len(parts) < 6:
                continue
            k = (parts[0], parts[1])
            if k not in wanted:
                continue
            try:
                n_sig = int(parts[2])
                n_sig_tot = int(parts[3])
                n_base = int(parts[4])
                n_base_tot = int(parts[5])
            except ValueError:
                continue
            if n_sig_tot == 0 or n_base_tot == 0:
                continue
            p_sig = n_sig / n_sig_tot
            p_base = n_base / n_base_tot
            out[k] = math.log10((p_sig + eps) / (p_base + eps))
    return out


def load_ec_aliases(path):
    """{EC -> {SEED reactions}}"""
    out = collections.defaultdict(set)
    with open(path) as f:
        f.readline()
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            ec = parts[2].strip()
            if not ec or '-' in ec:
                continue
            for rxn in parts[0].split('|'):
                rxn = rxn.strip()
                if rxn:
                    out[ec].add(rxn)
    return dict(out)


def load_shortlist_keys(shortlist_path, ec_to_rxns):
    """For each shortlist row expand gap_ec → SEED equivalents and
    emit (og, seed_rxn) keys. These are the keys the catalog uses.
    Returns (keys, meta, gap_key_to_seed_keys).
    """
    keys = set()
    meta = {}  # (og, seed_rxn) -> (reaction_name, ec_name, best_log_lr)
    gap_to_seeds = collections.defaultdict(set)  # (og, gap_rxn, ec) -> {seed}
    with open(shortlist_path) as f:
        rdr = csv.DictReader(f, delimiter='\t')
        for r in rdr:
            og = r.get('orthogroup', '')
            gap_rxn = r.get('gap_rxn', '')
            ec = r.get('gap_ec', '')
            if not og or og == 'unclustered':
                continue
            try:
                lr = float(r.get('log_lr', '0'))
            except ValueError:
                lr = 0.0

            seed_set = set(ec_to_rxns.get(ec, set())) if ec else set()
            if not seed_set:
                continue
            for sr in seed_set:
                k = (og, sr)
                keys.add(k)
                prev = meta.get(k, ('', '', -10))
                if lr > prev[2]:
                    meta[k] = (r.get('reaction_name', ''),
                                r.get('ec_name', ''), lr)
                gap_to_seeds[(og, gap_rxn, ec)].add(sr)
    return keys, meta, gap_to_seeds


def classify(lrs):
    """lrs = (lr_29, lr_60, lr_97). Return a label."""
    a, b, c = lrs
    # Treat differences < 0.1 as flat
    d1 = b - a
    d2 = c - b
    if abs(d1) < 0.1 and abs(d2) < 0.1:
        return 'flat'
    if d1 > 0 and d2 > 0:
        return 'monotonic_up' if c - a > 0.2 else 'flat'
    if d1 < 0 and d2 < 0:
        return 'monotonic_down'
    if abs(d2) < 0.1:
        return 'plateau'
    return 'noisy'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--shortlist', required=True)
    ap.add_argument('--n29', required=True)
    ap.add_argument('--n60', required=True)
    ap.add_argument('--n97', required=True)
    ap.add_argument('--ec-aliases', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    print(f'[0/4] loading EC aliases...', file=sys.stderr)
    ec_to_rxns = load_ec_aliases(args.ec_aliases)
    print(f'  {len(ec_to_rxns)} ECs', file=sys.stderr)

    print(f'[1/4] collecting shortlist keys...', file=sys.stderr)
    keys, meta, _ = load_shortlist_keys(args.shortlist, ec_to_rxns)
    print(f'  {len(keys)} unique (og, rxn) pairs in shortlist',
          file=sys.stderr)

    print(f'[2/4] streaming N=29 catalog for those keys...',
          file=sys.stderr)
    lr29 = load_catalog_subset(args.n29, keys)
    print(f'  {len(lr29)} matched', file=sys.stderr)

    print(f'[3/4] streaming N=60...', file=sys.stderr)
    lr60 = load_catalog_subset(args.n60, keys)
    print(f'  {len(lr60)} matched', file=sys.stderr)

    print(f'[4/4] streaming N=97...', file=sys.stderr)
    lr97 = load_catalog_subset(args.n97, keys)
    print(f'  {len(lr97)} matched', file=sys.stderr)

    # Join
    rows = []
    for k in keys:
        og, rxn = k
        r_name, ec_name, best_lr = meta[k]
        a = lr29.get(k)
        b = lr60.get(k)
        c = lr97.get(k)
        if a is None and b is None and c is None:
            continue
        lrs = (a if a is not None else 0.0,
               b if b is not None else 0.0,
               c if c is not None else 0.0)
        n_present = sum(1 for x in (a, b, c) if x is not None)
        label = classify(lrs) if n_present == 3 else 'sparse'
        rows.append([
            og, rxn, r_name, ec_name, n_present,
            f'{lrs[0]:.3f}', f'{lrs[1]:.3f}', f'{lrs[2]:.3f}',
            label,
        ])

    labels = collections.Counter(r[-1] for r in rows)
    print(f'\n[summary] {len(rows)} pairs with data in >= 1 subpanel:',
          file=sys.stderr)
    for k, v in labels.most_common():
        print(f'  {k}: {v}', file=sys.stderr)

    rows.sort(key=lambda r: (-float(r[7]),  # N=97 log_lr desc
                              r[-1]))
    with open(args.out, 'w') as f:
        f.write('orthogroup\treaction_id\treaction_name\tec_name\t'
                'n_present\tlog_lr_n29\tlog_lr_n60\tlog_lr_n97\t'
                'decay_label\n')
        for r in rows:
            f.write('\t'.join(str(x) for x in r) + '\n')
    print(f'  wrote {args.out}', file=sys.stderr)


if __name__ == '__main__':
    main()
