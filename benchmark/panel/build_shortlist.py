#!/usr/bin/env python3
"""Phase 5: build the master dark-matter shortlist across the 97-
genome panel.

Inputs
------
- Phase 4 per-genome augmented predictions
  (phase4/<tag>/dark_matter_augmented.tsv)
- Per-genome DIAMOND blastp results against Swiss-Prot
  (phase2/<tag>/preds/diamond_results.tsv), used to compute max pident
  per candidate so we can flag dark-matter (<30%).
- genome_manifest.tsv for phylum assignment.
- seed_reactions.tsv + ec2go for reaction/EC names.
- Per-genome gapsmith Reactions.tbl for culture-specific names.

Filter
------
rank <= 3 AND log_lr >= 0.3 AND max_pident < 30

Output
------
One TSV per phylum + one combined master, sorted by log_lr desc:
  culture phylum candidate rank log_lr n_hits orthogroup
  gap_rxn reaction_name gap_ec ec_name density n_anchors
  max_pident gene_context
"""
import argparse
import collections
import csv
import re
import sys
from pathlib import Path


def load_manifest_phylum(path):
    out = {}
    with open(path) as f:
        rdr = csv.DictReader(f, delimiter='\t')
        for r in rdr:
            out[r['genome_id']] = r.get('phylum') or 'Unknown'
    return out


def load_ec2go_names(path):
    """Parse ec2go: 'EC:X.Y.Z.W > GO:<name> ; GO:<id>' → {EC: name}."""
    out = {}
    with open(path) as f:
        for line in f:
            m = re.match(r'^EC:(\S+)\s*>\s*GO:([^;]+)\s*;', line)
            if m:
                out[m.group(1)] = m.group(2).strip()
                out[m.group(1).replace('EC:', '')] = m.group(2).strip()
    return out


def load_seed_reaction_names(path):
    """seed_reactions.tsv: first column id, 3rd column name."""
    out = {}
    with open(path) as f:
        f.readline()
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 3:
                out[parts[0]] = parts[2]
    return out


def load_culture_rxn_names(tbl):
    """Per-culture gapsmith Reactions.tbl column 3 has reaction names."""
    out = {}
    if not Path(tbl).exists():
        return out
    with open(tbl) as f:
        rdr = csv.DictReader(f, delimiter='\t')
        for r in rdr:
            rxn = r.get('rxn')
            name = r.get('name')
            if rxn and name and rxn not in out:
                out[rxn] = name
    return out


def load_max_pident(diamond_path):
    """{query_id: max_pident} from a DIAMOND m6 with pident in col 3."""
    out = {}
    if not Path(diamond_path).exists():
        return out
    with open(diamond_path) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            q = parts[0]
            try:
                p = float(parts[2])
            except ValueError:
                continue
            if q not in out or p > out[q]:
                out[q] = p
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--manifest', required=True)
    ap.add_argument('--genome-manifest', required=True)
    ap.add_argument('--phase4-root', required=True)
    ap.add_argument('--phase2-root', required=True)
    ap.add_argument('--ec2go', required=True)
    ap.add_argument('--seed-reactions', required=True)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--max-rank', type=int, default=3)
    ap.add_argument('--min-log-lr', type=float, default=0.3)
    ap.add_argument('--max-pident', type=float, default=30.0)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    phylum_of = load_manifest_phylum(args.genome_manifest)
    ec_names = load_ec2go_names(args.ec2go)
    seed_names = load_seed_reaction_names(args.seed_reactions)
    print(f'[info] {len(phylum_of)} genomes in manifest; '
          f'{len(ec_names)} EC names; {len(seed_names)} SEED names',
          file=sys.stderr)

    with open(args.manifest) as f:
        f.readline()
        tags = [line.split('\t')[0].strip() for line in f if line.strip()]
    print(f'[info] processing {len(tags)} tags', file=sys.stderr)

    header = ['culture', 'phylum', 'candidate', 'rank', 'log_lr',
              'n_hits', 'orthogroup', 'gap_rxn', 'reaction_name',
              'gap_ec', 'ec_name', 'density', 'n_anchors',
              'max_pident', 'gene_context']

    rows_by_phylum = collections.defaultdict(list)
    rows_all = []

    for tag in tags:
        pred = Path(args.phase4_root) / tag / 'dark_matter_augmented.tsv'
        if not pred.exists():
            continue

        # Per-culture reaction names + per-candidate max pident
        culture_names = load_culture_rxn_names(
            Path(args.phase2_root) / tag / 'gapsmith' /
            f'{tag}-all-Reactions.tbl')
        diamond = Path(args.phase2_root) / tag / 'preds' / \
            'diamond_results.tsv'
        pident_of = load_max_pident(str(diamond))
        phylum = phylum_of.get(tag, 'Unknown')

        with open(pred) as fin:
            hdr = fin.readline().rstrip('\n').split('\t')
            idx = {c: i for i, c in enumerate(hdr)}
            cand_i = idx.get('candidate')
            rank_i = idx.get('rank')
            log_lr_i = idx.get('log_lr')
            n_hits_i = idx.get('n_hits')
            og_i = idx.get('orthogroup')
            gap_rxn_i = idx.get('gap_rxn')
            gap_ec_i = idx.get('gap_ec')
            density_i = idx.get('density')
            n_anc_i = idx.get('n_anchors')
            gctx_i = idx.get('gene_context')

            for line in fin:
                parts = line.rstrip('\n').split('\t')
                if len(parts) <= log_lr_i:
                    continue
                try:
                    rank = int(parts[rank_i])
                except ValueError:
                    continue
                if rank > args.max_rank:
                    continue
                lr_s = parts[log_lr_i]
                if not lr_s:
                    continue
                try:
                    lr = float(lr_s)
                except ValueError:
                    continue
                if lr < args.min_log_lr:
                    continue

                cand = parts[cand_i]
                max_p = pident_of.get(cand, 0.0)
                if max_p >= args.max_pident:
                    continue  # not dark matter

                gap_rxn = parts[gap_rxn_i]
                gap_ec = parts[gap_ec_i]
                rname = culture_names.get(gap_rxn) or \
                    seed_names.get(gap_rxn) or ''
                ename = ec_names.get(gap_ec) or ec_names.get(
                    f'EC:{gap_ec}') or ''

                row = [
                    tag, phylum, cand, str(rank), f'{lr:.3f}',
                    parts[n_hits_i] if n_hits_i is not None else '',
                    parts[og_i] if og_i is not None else '',
                    gap_rxn, rname, gap_ec, ename,
                    parts[density_i] if density_i is not None else '',
                    parts[n_anc_i] if n_anc_i is not None else '',
                    f'{max_p:.1f}',
                    parts[gctx_i] if gctx_i is not None else '',
                ]
                # Dedupe key: one row per (tag, candidate, gap_rxn,
                # gap_ec). Keep the highest log_lr.
                rows_by_phylum[phylum].append(row)
                rows_all.append(row)

    # Dedupe across rows_all (and per-phylum) keeping the row with
    # the highest log_lr for each (tag, candidate, gap_rxn, gap_ec).
    def dedupe(rows):
        best = {}
        for r in rows:
            key = (r[0], r[2], r[7], r[9])
            try:
                lr = float(r[4])
            except ValueError:
                lr = 0.0
            if key not in best or lr > best[key][0]:
                best[key] = (lr, r)
        return [v[1] for v in best.values()]

    for phylum in list(rows_by_phylum):
        rows_by_phylum[phylum] = dedupe(rows_by_phylum[phylum])
    rows_all = dedupe(rows_all)

    # Sort each output by log_lr desc
    def sort_key(r):
        try:
            return -float(r[4])
        except ValueError:
            return 0.0

    for phylum, rows in rows_by_phylum.items():
        rows.sort(key=sort_key)
        safe = phylum.replace(' ', '_').replace('/', '_')
        out = out_dir / f'shortlist_{safe}.tsv'
        with open(out, 'w') as f:
            f.write('\t'.join(header) + '\n')
            for r in rows:
                f.write('\t'.join(r) + '\n')
        print(f'  {phylum}: {len(rows)} rows → {out}', file=sys.stderr)

    rows_all.sort(key=sort_key)
    out = out_dir / 'shortlist_master.tsv'
    with open(out, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for r in rows_all:
            f.write('\t'.join(r) + '\n')
    print(f'[done] master: {len(rows_all)} rows → {out}', file=sys.stderr)


if __name__ == '__main__':
    main()
