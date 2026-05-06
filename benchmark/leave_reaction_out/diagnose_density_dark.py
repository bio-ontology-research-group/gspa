#!/usr/bin/env python3
"""Diagnostic (D): for each dark LRO case, compute kernel density
D(pos(g), R) at *every* gene position in mg1655 and report the rank
of the true dark target among all genes.

If the dark target's density-rank is consistently good (top 20), the
signal exists and the RLC suggester's current candidate filter is
throwing it away; candidate generation should include these positions.
If the rank is consistently poor, no amount of ranker engineering on
top of density will recover dark matter — the signal isn't there.
"""
import argparse
import collections
import json
import math
import re
import sys
from pathlib import Path


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


def load_integrated(path, tau=0.3):
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
            if p <= tau:
                continue
            out[parts[0]][parts[2]] = p
    return dict(out)


def gauss(dx, h):
    z = dx / h
    return math.exp(-0.5 * z * z)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--reactions-tsv', required=True)
    ap.add_argument('--diffusion-tsv', required=True)
    ap.add_argument('--ec-aliases-tsv', required=True)
    ap.add_argument('--ec2go', required=True)
    ap.add_argument('--integrated', required=True,
                    help='mg1655_integrated.tsv — anchor posteriors')
    ap.add_argument('--layout', required=True)
    ap.add_argument('--cases', required=True,
                    help='Dark cases TSV')
    ap.add_argument('--bandwidth', type=float, default=5000.0)
    ap.add_argument('--radius-k', type=int, default=2)
    ap.add_argument('--alpha', type=float, default=0.5)
    ap.add_argument('--tau', type=float, default=0.5,
                    help='anchor posterior threshold')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    ec2go = load_ec2go(args.ec2go)
    rxn_to_mets, currency, met2rxn, rxn_to_ecs, ec_to_rxns = \
        load_reaction_graph(args.reactions_tsv, args.diffusion_tsv,
                             args.ec_aliases_tsv)
    print(f'[info] graph: {len(rxn_to_mets)} rxns, '
          f'{len(currency)} currency, {len(rxn_to_ecs)} EC-bearing',
          file=sys.stderr)

    layout, by_contig = load_layout(args.layout)
    print(f'[info] layout: {len(layout)} proteins', file=sys.stderr)

    posteriors = load_integrated(args.integrated, tau=0.0)
    print(f'[info] integrated posteriors for '
          f'{sum(len(v) for v in posteriors.values())} (prot, GO) '
          f'across {len(posteriors)} proteins', file=sys.stderr)

    # Read dark cases
    cases = []
    with open(args.cases) as f:
        header = f.readline().rstrip('\n').split('\t')
        for line in f:
            parts = line.rstrip('\n').split('\t')
            cases.append(dict(zip(header, parts)))
    print(f'[info] {len(cases)} dark cases', file=sys.stderr)

    h = args.bandwidth
    halfwidth = 4 * h
    with open(args.out, 'w') as fout:
        fout.write('\t'.join(['case_idx', 'protein_id', 'reaction_id',
                              'go_term', 'max_pident',
                              'n_nbr_rxn', 'n_nbr_gos', 'n_anchors',
                              'target_density', 'target_rank',
                              'total_genes', 'top1_gene',
                              'top1_density']) + '\n')
        for i, c in enumerate(cases):
            target = c['protein_id']
            target_rxn = c['reaction_id']
            target_go = c['go_term']
            max_pid = c.get('max_pident_annot_homolog', '')
            target_loc = layout.get(target)
            if target_loc is None:
                continue
            target_contig, tstart, tend, tstrand, target_mid = target_loc

            # Resolve target reaction to SEED reaction IDs via EC bridge
            # (target_rxn is usually MetaCyc; SEED graph uses rxnNNNNN).
            seed_rxns = set()
            if target_rxn in rxn_to_mets:
                seed_rxns.add(target_rxn)
            target_ec = c.get('ec', '').strip()
            if target_ec:
                seed_rxns.update(ec_to_rxns.get(target_ec, set()))

            # BFS from each SEED equivalent; union neighbour weights.
            nbrs = {}
            for sr in seed_rxns:
                sub = bfs_neighbors(sr, rxn_to_mets, met2rxn,
                                      currency, args.radius_k, args.alpha)
                for r, w in sub.items():
                    if r in seed_rxns:
                        continue
                    cur = nbrs.get(r, 0.0)
                    if w > cur:
                        nbrs[r] = w

            # Map to neighbour GOs with weights
            nbr_go_weights = {}
            for n_rxn, w in nbrs.items():
                for ec in rxn_to_ecs.get(n_rxn, []):
                    go = ec2go.get(ec) or ec2go.get(f'EC:{ec}')
                    if not go:
                        continue
                    cur = nbr_go_weights.get(go)
                    if cur is None or w > cur:
                        nbr_go_weights[go] = w

            # Build anchor list (protein_id -> summed weight)
            # Anchor = any protein (NOT the target) with posterior on any
            # neighbour GO above tau.
            anchor_weights = collections.defaultdict(float)
            for prot, gos in posteriors.items():
                if prot == target:
                    continue
                best = 0.0
                for go, w_r in nbr_go_weights.items():
                    p = gos.get(go, 0.0)
                    if p < args.tau:
                        continue
                    contrib = w_r * p
                    if contrib > best:
                        best = contrib
                if best > 0:
                    anchor_weights[prot] = best

            # Compute density at every gene position on the target's
            # contig (cross-contig anchors contribute 0).
            # Use sorted anchor list for efficiency.
            anchor_positions = []
            for aprot, aw in anchor_weights.items():
                aloc = layout.get(aprot)
                if not aloc or aloc[0] != target_contig:
                    continue
                amid = aloc[4]
                anchor_positions.append((amid, aw))
            anchor_positions.sort()

            gene_density = {}
            genes_on_contig = by_contig.get(target_contig, [])
            for gmid, gpid in genes_on_contig:
                d = 0.0
                for (amid, aw) in anchor_positions:
                    if amid - gmid > halfwidth:
                        break
                    if gmid - amid > halfwidth:
                        continue
                    d += aw * gauss(abs(amid - gmid), h)
                gene_density[gpid] = d

            # Rank target
            sorted_g = sorted(gene_density.items(), key=lambda x: -x[1])
            rank = 0
            tden = gene_density.get(target, 0.0)
            for r, (g, d) in enumerate(sorted_g, 1):
                if g == target:
                    rank = r
                    break
            top1_gene = sorted_g[0][0] if sorted_g else ''
            top1_d = sorted_g[0][1] if sorted_g else 0.0

            fout.write('\t'.join([
                str(i), target, target_rxn, target_go, str(max_pid),
                str(len(nbrs)), str(len(nbr_go_weights)),
                str(len(anchor_positions)),
                f'{tden:.4f}', str(rank),
                str(len(gene_density)),
                top1_gene, f'{top1_d:.4f}',
            ]) + '\n')
            if (i + 1) % 5 == 0:
                print(f'  [{i+1}/{len(cases)}]', file=sys.stderr)

    print(f'[info] wrote {args.out}', file=sys.stderr)


if __name__ == '__main__':
    main()
