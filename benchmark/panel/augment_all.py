#!/usr/bin/env python3
"""Phase 4B: augment every query genome's dark_matter.tsv with
cross-genome log_lr from the panel non-anchor catalog, using a single
catalog load.

Strategy: pre-filter the catalog to entries with log_lr >= <threshold>
(default 0.3, the "meaningful" enrichment floor) — this cuts memory by
~1000x (128M rows → ~127K rows). Entries below threshold are treated
as log_lr=0.

For each genome:
  - Query candidate = protein ID (no mmseqs needed; panel == query)
  - Orthogroup = direct lookup of f"{tag}:{candidate}" in the ortho map
  - SEED equivalents of gap_rxn via EC
  - Sum catalog entries across SEED equivalents; emit log_lr
"""
import argparse
import collections
import csv
import math
import sys
from pathlib import Path


def load_ortho_map(path):
    out = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            # We want: full tag-prefixed member -> rep (stripped or not).
            # The catalog's orthogroup_id column uses the MMseqs rep's
            # plain pid (tag stripped). So normalize rep the same way.
            member, rep = parts[0], parts[1]
            if ':' in rep:
                rep = rep.split(':', 1)[1]
            out[member] = rep
    return out


def load_ec_aliases(path):
    ec_to_rxns = collections.defaultdict(set)
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
                    ec_to_rxns[ec].add(rxn)
    return dict(ec_to_rxns)


def load_catalog_filtered(path, min_log_lr, eps=1e-3):
    """Load only rows where log_lr >= min_log_lr. Returns
    {(og, rxn): (n_sig, n_sig_tot, n_base, n_base_tot, log_lr)}."""
    out = {}
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if parts[0] == 'orthogroup_id':
                continue
            if len(parts) < 6:
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
            lr = math.log10((p_sig + eps) / (p_base + eps))
            if lr < min_log_lr:
                continue
            out[(parts[0], parts[1])] = (n_sig, n_sig_tot, n_base,
                                          n_base_tot, lr)
    return out


def augment_one(pred_path, out_path, tag, ortho, catalog, ec_to_rxns):
    rxn_cov = {r for (_, r) in catalog.keys()}
    n_in = 0
    n_orthogrouped = 0
    n_match = 0
    with open(pred_path) as fin, open(out_path, 'w') as fout:
        hdr = fin.readline().rstrip('\n').split('\t')
        new_hdr = hdr + ['orthogroup', 'seed_rxns', 'n_sig_nonanc',
                          'n_sig_total', 'n_base_with', 'n_base_total',
                          'log_lr', 'n_hits']
        fout.write('\t'.join(new_hdr) + '\n')
        try:
            cand_i = hdr.index('candidate')
        except ValueError:
            cand_i = hdr.index('gene') if 'gene' in hdr else 5
        rxn_i = hdr.index('gap_rxn')
        ec_i = hdr.index('gap_ec')

        for line in fin:
            parts = line.rstrip('\n').split('\t')
            n_in += 1
            cand = parts[cand_i]
            gap_rxn = parts[rxn_i]
            gap_ec = parts[ec_i]

            og = ortho.get(f'{tag}:{cand}') or ortho.get(cand)
            if og:
                n_orthogrouped += 1
            og = og or 'unclustered'

            seed_rxns = set()
            if gap_rxn in rxn_cov:
                seed_rxns.add(gap_rxn)
            if gap_ec:
                seed_rxns |= ec_to_rxns.get(gap_ec, set())
            seed_rxns &= rxn_cov

            n_sig_sum = 0
            n_sig_tot_max = 0
            n_base_sum = 0
            n_base_tot = 0
            best_lr = 0.0
            n_hits = 0
            for sr in seed_rxns:
                entry = catalog.get((og, sr))
                if not entry:
                    continue
                n_hits += 1
                a, b, c, d, lr = entry
                n_sig_sum += a
                n_sig_tot_max = max(n_sig_tot_max, b)
                n_base_sum = max(n_base_sum, c)
                n_base_tot = max(n_base_tot, d)
                if lr > best_lr:
                    best_lr = lr
            if n_hits:
                n_match += 1

            parts.extend([
                og,
                ','.join(sorted(seed_rxns)) if seed_rxns else '',
                str(n_sig_sum) if n_hits else '',
                str(n_sig_tot_max) if n_hits else '',
                str(n_base_sum) if n_hits else '',
                str(n_base_tot) if n_hits else '',
                f'{best_lr:.3f}' if n_hits else '',
                str(n_hits),
            ])
            fout.write('\t'.join(parts) + '\n')
    return n_in, n_orthogrouped, n_match


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--manifest', required=True)
    ap.add_argument('--phase4-root', required=True)
    ap.add_argument('--catalog', required=True)
    ap.add_argument('--ortho-map', required=True)
    ap.add_argument('--ec-aliases', required=True)
    ap.add_argument('--min-log-lr', type=float, default=0.3)
    args = ap.parse_args()

    print(f'[1/3] loading ortho map...', file=sys.stderr)
    ortho = load_ortho_map(args.ortho_map)
    print(f'  {len(ortho)} entries', file=sys.stderr)

    print(f'[2/3] loading ec aliases...', file=sys.stderr)
    ec_to_rxns = load_ec_aliases(args.ec_aliases)
    print(f'  {len(ec_to_rxns)} ECs', file=sys.stderr)

    print(f'[3/3] loading catalog (log_lr >= {args.min_log_lr})...',
          file=sys.stderr)
    catalog = load_catalog_filtered(args.catalog, args.min_log_lr)
    print(f'  {len(catalog)} filtered entries', file=sys.stderr)

    root = Path(args.phase4_root)
    with open(args.manifest) as f:
        f.readline()
        tags = [line.split('\t')[0].strip() for line in f if line.strip()]
    print(f'[augment] {len(tags)} genomes', file=sys.stderr)

    totals = [0, 0, 0]
    for tag in tags:
        pred = root / tag / 'dark_matter.tsv'
        out = root / tag / 'dark_matter_augmented.tsv'
        if not pred.exists():
            print(f'  [skip] {tag}: no dark_matter.tsv', file=sys.stderr)
            continue
        n_in, n_og, n_match = augment_one(str(pred), str(out), tag,
                                           ortho, catalog, ec_to_rxns)
        totals[0] += n_in
        totals[1] += n_og
        totals[2] += n_match

    print(f'[done] rows: {totals[0]}; orthogrouped: {totals[1]} '
          f'({100*totals[1]/max(1,totals[0]):.1f}%); catalog-match: '
          f'{totals[2]} ({100*totals[2]/max(1,totals[0]):.1f}%)',
          file=sys.stderr)


if __name__ == '__main__':
    main()
