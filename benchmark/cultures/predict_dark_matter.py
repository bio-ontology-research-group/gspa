#!/usr/bin/env python3
"""Prospective dark-matter predictor.

For a culture genome, enumerate reactions gapsmith failed to fill
(bad_blast / no_blast in Reactions.tbl, with EC + pathway context) and,
for each, rank candidate proteins by kernel density of neighbour-
reaction anchors, with the anchors themselves + any gapsmith-assigned
protein removed from the candidate pool.

Output TSV (one row per (gap, candidate)):
  culture gap_rxn gap_pathway gap_ec rank candidate density n_anchors
  n_nbr_gos gene_context
"""
import argparse
import collections
import csv
import math
import re
import sys
from pathlib import Path


# --- Reuse graph/layout/integrated helpers from diagnose_density_nonanchor ---


def load_ec2go(path):
    out = {}
    with open(path) as f:
        for line in f:
            if line.startswith('!') or not line.strip():
                continue
            m = re.match(r'^EC:(\S+)\s*>\s*GO:[^;]+;\s*(GO:\d+)', line)
            if m:
                out[m.group(1)] = m.group(2)
                out[m.group(1)[3:] if m.group(1).startswith('EC:')
                     else m.group(1)] = m.group(2)
    return out


def load_reaction_graph(reactions_tsv, diffusion_tsv, ec_aliases_tsv,
                          currency_percentile=99.0):
    currency = set()
    if Path(diffusion_tsv).exists():
        with open(diffusion_tsv) as f:
            f.readline()
            for line in f:
                parts = line.split('\t')
                if parts and parts[0].startswith('cpd'):
                    currency.add(parts[0].strip())
    rxn_to_mets = {}
    raw_deg = collections.Counter()
    with open(reactions_tsv) as f:
        header = f.readline().rstrip('\n').split('\t')
        id_i = header.index('id')
        try:
            stoich_i = header.index('stoichiometry')
        except ValueError:
            stoich_i = header.index('equation')
        for line in f:
            parts = line.split('\t')
            if len(parts) <= stoich_i:
                continue
            rxn = parts[id_i].strip()
            mets = set()
            for entry in parts[stoich_i].split(';'):
                entry = entry.strip()
                if not entry:
                    continue
                fields = entry.split(':', 4)
                if len(fields) < 2:
                    continue
                try:
                    float(fields[0])
                except ValueError:
                    continue
                cpd = fields[1].strip()
                if cpd:
                    mets.add(cpd)
            rxn_to_mets[rxn] = mets
            for m in mets:
                raw_deg[m] += 1
    if currency_percentile < 100.0 and raw_deg:
        sorted_deg = sorted(raw_deg.values())
        idx = max(0, min(len(sorted_deg) - 1,
                         int(math.ceil(len(sorted_deg) *
                                        currency_percentile / 100.0)) - 1))
        threshold = sorted_deg[idx]
        for mm, d in raw_deg.items():
            if d >= threshold:
                currency.add(mm)
    metabolite_to_reactions = collections.defaultdict(set)
    for rxn, mets in rxn_to_mets.items():
        for mm in mets:
            if mm not in currency:
                metabolite_to_reactions[mm].add(rxn)
    rxn_to_ecs = collections.defaultdict(set)
    ec_to_rxns = collections.defaultdict(set)
    if Path(ec_aliases_tsv).exists():
        with open(ec_aliases_tsv) as f:
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
                    if rxn and rxn in rxn_to_mets:
                        rxn_to_ecs[rxn].add(ec)
                        ec_to_rxns[ec].add(rxn)
    return (rxn_to_mets, currency, dict(metabolite_to_reactions),
            dict(rxn_to_ecs), dict(ec_to_rxns))


def bfs_neighbors(rxn_id, rxn_to_mets, metabolite_to_reactions,
                   currency, max_k, alpha):
    out = {}
    seen = {rxn_id}
    frontier = [rxn_id]
    weight = 1.0
    for d in range(1, max_k + 1):
        weight *= alpha
        nxt = []
        for f in frontier:
            for mm in rxn_to_mets.get(f, set()):
                if mm in currency:
                    continue
                for n in metabolite_to_reactions.get(mm, set()):
                    if n not in seen:
                        seen.add(n)
                        out[n] = weight
                        nxt.append(n)
        if not nxt:
            break
        frontier = nxt
    return out


def load_layout(path):
    out = {}
    by_contig = collections.defaultdict(list)
    with open(path) as f:
        f.readline()
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 5:
                continue
            pid, contig, start, end, strand = (parts[0], parts[1],
                                                 int(parts[2]),
                                                 int(parts[3]), parts[4])
            mid = (start + end) // 2
            out[pid] = (contig, start, end, strand, mid)
            by_contig[contig].append((mid, pid))
    for c in by_contig:
        by_contig[c].sort()
    return out, dict(by_contig)


def load_integrated(path):
    out = collections.defaultdict(dict)
    with open(path) as f:
        f.readline()
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 5 or parts[1] != 'GO':
                continue
            try:
                p = float(parts[4])
            except ValueError:
                continue
            out[parts[0]][parts[2]] = p
    return dict(out)


def gauss(dx, h):
    z = dx / h
    return math.exp(-0.5 * z * z)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tag', required=True)
    ap.add_argument('--reactions-tbl', required=True,
                    help='gapsmith Reactions.tbl')
    ap.add_argument('--layout', required=True)
    ap.add_argument('--integrated', required=True)
    ap.add_argument('--reactions-tsv', required=True)
    ap.add_argument('--diffusion-tsv', required=True)
    ap.add_argument('--ec-aliases-tsv', required=True)
    ap.add_argument('--ec2go', required=True)
    ap.add_argument('--bandwidth', type=float, default=5000.0)
    ap.add_argument('--radius-k', type=int, default=2)
    ap.add_argument('--alpha', type=float, default=0.5)
    ap.add_argument('--tau', type=float, default=0.3)
    ap.add_argument('--top-k', type=int, default=5)
    ap.add_argument('--min-anchors', type=int, default=3,
                    help='Skip gaps with fewer anchors than this.')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    ec2go = load_ec2go(args.ec2go)
    rxn_to_mets, currency, met2rxn, rxn_to_ecs, ec_to_rxns = \
        load_reaction_graph(args.reactions_tsv, args.diffusion_tsv,
                             args.ec_aliases_tsv)
    layout, by_contig = load_layout(args.layout)
    posteriors = load_integrated(args.integrated)
    print(f'[info] {len(layout)} proteins, {len(posteriors)} with posteriors',
          file=sys.stderr)

    # Parse Reactions.tbl. Identify: gaps (bad_blast/no_blast) with EC;
    # good_blast (assigned proteins to exclude globally).
    gap_rows = []  # (pathway, rxn, ec)
    gapsmith_assigned = set()
    gapsmith_for_rxn = collections.defaultdict(set)
    pathway_reactions = collections.defaultdict(dict)  # pathway -> rxn -> ec
    with open(args.reactions_tbl) as f:
        rdr = csv.DictReader(f, delimiter='\t')
        for row in rdr:
            pathway = row['pathway']
            rxn = row['rxn']
            ec = (row.get('ec') or '').strip()
            status = row.get('status') or ''
            if ec:
                pathway_reactions[pathway][rxn] = ec
            if status == 'good_blast':
                prot = (row.get('stitle') or '').strip()
                if prot:
                    gapsmith_assigned.add(prot)
                    gapsmith_for_rxn[rxn].add(prot)
            elif status in ('bad_blast', 'no_blast') and ec:
                gap_rows.append((pathway, rxn, ec))

    print(f'[info] {len(gap_rows)} gap rows, '
          f'{len(gapsmith_assigned)} good_blast assignments',
          file=sys.stderr)

    h = args.bandwidth
    halfwidth = 4 * h
    out_rows = []
    for pathway, rxn, ec in gap_rows:
        seed_rxns = set()
        if rxn in rxn_to_mets:
            seed_rxns.add(rxn)
        seed_rxns.update(ec_to_rxns.get(ec, set()))
        if not seed_rxns:
            continue

        nbrs = {}
        for sr in seed_rxns:
            sub = bfs_neighbors(sr, rxn_to_mets, met2rxn, currency,
                                  args.radius_k, args.alpha)
            for r, w in sub.items():
                if r in seed_rxns:
                    continue
                cur = nbrs.get(r, 0.0)
                if w > cur:
                    nbrs[r] = w
        if not nbrs:
            continue

        nbr_go_weights = {}
        for n_rxn, w in nbrs.items():
            for ec_n in rxn_to_ecs.get(n_rxn, []):
                go = ec2go.get(ec_n) or ec2go.get(f'EC:{ec_n}')
                if not go:
                    continue
                cur = nbr_go_weights.get(go)
                if cur is None or w > cur:
                    nbr_go_weights[go] = w

        anchor_set = set()
        anchor_weights = collections.defaultdict(float)
        for prot, gos in posteriors.items():
            best = 0.0
            is_anchor = False
            for go, w_r in nbr_go_weights.items():
                p = gos.get(go, 0.0)
                if p >= args.tau:
                    is_anchor = True
                    contrib = w_r * p
                    if contrib > best:
                        best = contrib
            if is_anchor:
                anchor_set.add(prot)
                anchor_weights[prot] = best
        if len(anchor_set) < args.min_anchors:
            continue

        # Find all contigs where anchors sit; compute density there.
        contigs_with_anchors = set()
        for aprot in anchor_set:
            aloc = layout.get(aprot)
            if aloc:
                contigs_with_anchors.add(aloc[0])

        # For each anchor-bearing contig, rank genes by density.
        exclude = (anchor_set | gapsmith_assigned) \
                    - gapsmith_for_rxn.get(rxn, set())
        for contig in contigs_with_anchors:
            anchor_positions = []
            for aprot in anchor_set:
                aloc = layout.get(aprot)
                if not aloc or aloc[0] != contig:
                    continue
                anchor_positions.append((aloc[4], anchor_weights[aprot]))
            anchor_positions.sort()
            gene_density = {}
            for gmid, gpid in by_contig.get(contig, []):
                d = 0.0
                for (amid, aw) in anchor_positions:
                    if amid - gmid > halfwidth:
                        break
                    if gmid - amid > halfwidth:
                        continue
                    d += aw * gauss(abs(amid - gmid), h)
                gene_density[gpid] = d
            candidates = [(g, d) for g, d in gene_density.items()
                          if g not in exclude and d > 0]
            candidates.sort(key=lambda x: -x[1])
            for rank, (cand, den) in enumerate(candidates[:args.top_k], 1):
                cloc = layout.get(cand, ('', 0, 0, '', 0))
                # Context: list nearest 3 anchors and their weights
                near = sorted(anchor_positions,
                               key=lambda ap:
                                 abs(ap[0] - cloc[4]))[:3]
                near_desc = ';'.join(f'anc@{ap[0]}w={ap[1]:.3f}'
                                      for ap in near)
                out_rows.append({
                    'culture': args.tag, 'gap_pathway': pathway,
                    'gap_rxn': rxn, 'gap_ec': ec, 'rank': rank,
                    'candidate': cand, 'density': f'{den:.4f}',
                    'n_anchors': len(anchor_positions),
                    'n_nbr_gos': len(nbr_go_weights),
                    'contig': cloc[0], 'pos': cloc[4],
                    'near_anchors': near_desc,
                })

    print(f'[info] {len(out_rows)} candidate predictions across gaps',
          file=sys.stderr)
    cols = ['culture', 'gap_pathway', 'gap_rxn', 'gap_ec', 'rank',
             'candidate', 'density', 'n_anchors', 'n_nbr_gos',
             'contig', 'pos', 'near_anchors']
    with open(args.out, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in out_rows:
            f.write('\t'.join(str(r[c]) for c in cols) + '\n')
    print(f'[info] wrote {args.out}', file=sys.stderr)


if __name__ == '__main__':
    main()
