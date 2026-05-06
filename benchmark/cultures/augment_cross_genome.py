#!/usr/bin/env python3
"""Add cross-genome non-anchor evidence to culture predictions.

For each (culture, candidate, gap_rxn) row:
  1. Map candidate → best panel protein hit via MMseqs2 (≥50% id, ≥80%
     qcov). Map that panel protein → orthogroup via panel
     orthogroup_map_50.tsv. Missing hit ⇒ orthogroup = 'unclustered'
     (i.e. no cross-genome anchor possible by this clustering).
  2. Resolve gap_rxn → SEED-reaction IDs via the EC of the gap.
  3. Look up (orthogroup, SEED_rxn) in nonanchor_catalog_mg1655.tsv,
     summing contributions across the SEED equivalents of the gap.
  4. Append columns:
       orthogroup, n_sig_nonanc, n_sig_total, n_base_with, n_base_total,
       log_lr  (log10 of [(n_sig_nonanc/n_sig_total + eps) /
                          (n_base_with/n_base_total + eps)])

n_sig_total is the number of panel genomes (out of 29) that had an
R-signature peak for this reaction at all. n_sig_nonanc is how many of
those genomes had a non-anchor member of the orthogroup inside that
peak. High log_lr (>0) ⇒ the orthogroup appears in R-signature contexts
more often than its baseline panel prevalence, so cross-genome evidence
*supports* this candidate for this gap. log_lr near 0 ⇒ no signal.
Negative ⇒ the orthogroup is *under*-represented in gap contexts.
"""
import argparse
import collections
import csv
import math
import re
import sys
from pathlib import Path


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


def load_panel_ortho_map(path):
    """Return {panel_protein -> orthogroup_id}. Strip any 'tag:' prefix
    from the orthogroup rep."""
    out = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            member = parts[0]
            rep = parts[1]
            if ':' in rep:
                rep = rep.split(':', 1)[1]
            out[member] = rep
    return out


def load_mmseqs_hits(path):
    """Best hit per query from MMseqs m8 output. Query has 'tag:pid'
    prefix (added at FASTA concat time); strip off the tag."""
    out = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            q = parts[0]
            if ':' in q:
                tag, pid = q.split(':', 1)
            else:
                tag, pid = '', q
            target = parts[1]
            # target format: 'tag:accession' from all_panel_proteins.faa
            if ':' in target:
                target = target.split(':', 1)[1]
            try:
                pident = float(parts[2])
            except (IndexError, ValueError):
                pident = 0.0
            key = (tag, pid)
            if key not in out or pident > out[key][1]:
                out[key] = (target, pident)
    return out


def load_catalog(path):
    """{(orthogroup, reaction_id) -> (n_sig, n_sig_total, n_base,
                                      n_base_total)}"""
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
            og, rxn = parts[0], parts[1]
            n_sig = int(parts[2])
            n_sig_tot = int(parts[3])
            n_base = int(parts[4])
            n_base_tot = int(parts[5])
            out[(og, rxn)] = (n_sig, n_sig_tot, n_base, n_base_tot)
    return out


def catalog_rxn_coverage(catalog):
    """Set of reaction IDs that have at least one catalog entry."""
    return {r for (_, r) in catalog.keys()}


def log_lr(n_sig, n_sig_tot, n_base, n_base_tot, eps=1e-3):
    if n_sig_tot == 0 or n_base_tot == 0:
        return None
    p_sig = n_sig / n_sig_tot
    p_base = n_base / n_base_tot
    return math.log10((p_sig + eps) / (p_base + eps))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--predictions', required=True)
    ap.add_argument('--mmseqs-m8', required=True)
    ap.add_argument('--panel-ortho-map', required=True)
    ap.add_argument('--catalog', required=True)
    ap.add_argument('--ec-aliases', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    hits = load_mmseqs_hits(args.mmseqs_m8)
    ortho = load_panel_ortho_map(args.panel_ortho_map)
    catalog = load_catalog(args.catalog)
    rxn_cov = catalog_rxn_coverage(catalog)
    ec_to_rxns = load_ec_aliases(args.ec_aliases)
    print(f'[info] mmseqs hits: {len(hits)}; ortho map: {len(ortho)}; '
          f'catalog rows: {len(catalog)} over {len(rxn_cov)} reactions',
          file=sys.stderr)

    n_in = 0
    n_orthogrouped = 0
    n_catalog_match = 0
    n_lr_nonzero = 0
    with open(args.predictions) as fin, open(args.out, 'w') as fout:
        h = fin.readline().rstrip('\n').split('\t')
        new_h = h + ['orthogroup', 'seed_rxn_lookup',
                     'n_sig_nonanc', 'n_sig_total',
                     'n_base_with', 'n_base_total', 'log_lr',
                     'n_catalog_lookups']
        fout.write('\t'.join(new_h) + '\n')
        culture_i = h.index('culture')
        cand_i = h.index('candidate')
        rxn_i = h.index('gap_rxn')
        ec_i = h.index('gap_ec')
        for line in fin:
            parts = line.rstrip('\n').split('\t')
            n_in += 1
            tag = parts[culture_i]
            cand = parts[cand_i]
            gap_rxn = parts[rxn_i]
            gap_ec = parts[ec_i]

            # Step 1: orthogroup of candidate
            panel_hit = hits.get((tag, cand))
            if panel_hit is None:
                og = 'unclustered'
            else:
                og = ortho.get(panel_hit[0], 'unclustered')
            if og != 'unclustered':
                n_orthogrouped += 1

            # Step 2: SEED equivalents of gap reaction
            seed_rxns = set()
            if gap_rxn in rxn_cov:
                seed_rxns.add(gap_rxn)
            if gap_ec:
                seed_rxns |= ec_to_rxns.get(gap_ec, set())
            seed_rxns &= rxn_cov

            # Step 3: catalog lookup (sum over SEED equivalents)
            n_sig_sum = 0
            n_sig_tot_max = 0
            n_base_sum = 0
            n_base_tot = 0
            n_hits = 0
            for sr in seed_rxns:
                entry = catalog.get((og, sr))
                if not entry:
                    continue
                n_hits += 1
                a, b, c, d = entry
                n_sig_sum += a
                if b > n_sig_tot_max:
                    n_sig_tot_max = b
                if c > n_base_sum:
                    n_base_sum = c
                if d > n_base_tot:
                    n_base_tot = d
            lr = log_lr(n_sig_sum, n_sig_tot_max, n_base_sum, n_base_tot) \
                if n_hits > 0 else None
            if n_hits > 0:
                n_catalog_match += 1
            if lr is not None and abs(lr) > 0.1:
                n_lr_nonzero += 1

            parts.extend([
                og,
                ','.join(sorted(seed_rxns)) if seed_rxns else '',
                str(n_sig_sum) if n_hits else '',
                str(n_sig_tot_max) if n_hits else '',
                str(n_base_sum) if n_hits else '',
                str(n_base_tot) if n_hits else '',
                f'{lr:.3f}' if lr is not None else '',
                str(n_hits),
            ])
            fout.write('\t'.join(parts) + '\n')

    print(f'\n[summary] {n_in} predictions in', file=sys.stderr)
    print(f'  orthogrouped (candidate has >=50% id hit in panel): '
          f'{n_orthogrouped} ({100 * n_orthogrouped / n_in:.1f}%)',
          file=sys.stderr)
    print(f'  catalog match (orthogroup x SEED rxn in catalog): '
          f'{n_catalog_match} ({100 * n_catalog_match / n_in:.1f}%)',
          file=sys.stderr)
    print(f'  meaningful LR (|log_lr| > 0.1): '
          f'{n_lr_nonzero} ({100 * n_lr_nonzero / n_in:.1f}%)',
          file=sys.stderr)


if __name__ == '__main__':
    main()
