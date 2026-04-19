#!/usr/bin/env python3
"""Diagnostic (1 / refined): for each dark case, rank mg1655 genes by
density D(pos(g), R), but *exclude* anchors from the candidate pool.

Anchor = any mg1655 protein (not the target) whose integrated posterior
on any neighbour-reaction GO exceeds tau. Candidates = all other genes
on the target's contig, ranked by D. Reports rank of the true target.

If anchor-exclusion materially improves dark rank vs the unfiltered
diagnostic, the "non-anchor within an anchor cluster" framing is the
right signal for dark-matter prediction.
"""
import argparse
import collections
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


def load_integrated(path, tau=0.0):
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
    ap.add_argument('--integrated', required=True)
    ap.add_argument('--layout', required=True)
    ap.add_argument('--cases', required=True)
    ap.add_argument('--bandwidth', type=float, default=5000.0)
    ap.add_argument('--radius-k', type=int, default=2)
    ap.add_argument('--alpha', type=float, default=0.5)
    ap.add_argument('--tau', type=float, default=0.3,
                    help='anchor posterior threshold')
    ap.add_argument('--reactions-tbl', default=None,
                    help='gapsmith Reactions.tbl — excludes any good_blast '
                         'protein from candidate pool (option D).')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    ec2go = load_ec2go(args.ec2go)
    rxn_to_mets, currency, met2rxn, rxn_to_ecs, ec_to_rxns = \
        load_reaction_graph(args.reactions_tsv, args.diffusion_tsv,
                             args.ec_aliases_tsv)
    layout, by_contig = load_layout(args.layout)
    posteriors = load_integrated(args.integrated, tau=0.0)
    print(f'[info] graph: {len(rxn_to_mets)} rxns; layout: {len(layout)} '
          f'prots; posteriors for {len(posteriors)} prots',
          file=sys.stderr)

    cases = []
    with open(args.cases) as f:
        header = f.readline().rstrip('\n').split('\t')
        for line in f:
            parts = line.rstrip('\n').split('\t')
            cases.append(dict(zip(header, parts)))
    print(f'[info] {len(cases)} cases', file=sys.stderr)

    # Gapsmith good_blast proteins (globally assigned = "already
    # occupied" and therefore not dark-candidate material).
    gapsmith_assigned = set()
    # Also: per-reaction, the proteins gapsmith assigned *to that reaction*.
    gapsmith_for_rxn = collections.defaultdict(set)
    if args.reactions_tbl:
        import csv as _csv
        with open(args.reactions_tbl) as f:
            rdr = _csv.DictReader(f, delimiter='\t')
            for row in rdr:
                if row.get('status') != 'good_blast':
                    continue
                stitle = (row.get('stitle') or '').strip()
                if not stitle:
                    continue
                gapsmith_assigned.add(stitle)
                gapsmith_for_rxn[row['rxn']].add(stitle)
        print(f'[info] gapsmith good_blast proteins: '
              f'{len(gapsmith_assigned)}', file=sys.stderr)

    h = args.bandwidth
    halfwidth = 4 * h
    with open(args.out, 'w') as fout:
        fout.write('\t'.join(['case_idx', 'protein_id', 'reaction_id',
                              'go_term', 'max_pident',
                              'n_nbr_rxn', 'n_nbr_gos', 'n_anchors',
                              'target_density', 'rank_all',
                              'rank_nonanchor', 'n_nonanchor',
                              'rank_unassigned', 'n_unassigned',
                              'top1_unassigned',
                              'target_is_gapsmith_assigned']) + '\n')
        for i, c in enumerate(cases):
            target = c['protein_id']
            target_rxn = c['reaction_id']
            target_go = c['go_term']
            max_pid = c.get('max_pident_annot_homolog', '')
            target_loc = layout.get(target)
            if target_loc is None:
                continue
            target_contig = target_loc[0]

            # Resolve target reaction to SEED rxns
            seed_rxns = set()
            if target_rxn in rxn_to_mets:
                seed_rxns.add(target_rxn)
            target_ec = c.get('ec', '').strip()
            if target_ec:
                seed_rxns.update(ec_to_rxns.get(target_ec, set()))

            # BFS neighbours
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

            # Map neighbours to GOs with weights
            nbr_go_weights = {}
            for n_rxn, w in nbrs.items():
                for ec in rxn_to_ecs.get(n_rxn, []):
                    go = ec2go.get(ec) or ec2go.get(f'EC:{ec}')
                    if not go:
                        continue
                    cur = nbr_go_weights.get(go)
                    if cur is None or w > cur:
                        nbr_go_weights[go] = w

            # Build anchor set: any protein (NOT target) with posterior
            # on any nbr GO >= tau. Also accumulate anchor weights for
            # density computation (weight = max_go w_R * posterior).
            anchor_set = set()
            anchor_weights = collections.defaultdict(float)
            for prot, gos in posteriors.items():
                if prot == target:
                    continue
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

            # Density at every gene position (same as before).
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

            # Rank A: all genes.
            sorted_all = sorted(gene_density.items(),
                                 key=lambda x: -x[1])
            rank_all = 0
            for r, (g, d) in enumerate(sorted_all, 1):
                if g == target:
                    rank_all = r
                    break
            tden = gene_density.get(target, 0.0)

            # Rank B: non-anchor only.
            non_anchor = [(g, d) for g, d in gene_density.items()
                           if g not in anchor_set]
            non_anchor.sort(key=lambda x: -x[1])
            rank_nonanc = 0
            for r, (g, d) in enumerate(non_anchor, 1):
                if g == target:
                    rank_nonanc = r
                    break
            top1_nonanc_gene = non_anchor[0][0] if non_anchor else ''

            # Rank C: exclude gapsmith-assigned globally (anywhere).
            # Do NOT exclude proteins gapsmith assigned *to this reaction
            # R* — those ARE the correct answer (if gapsmith put them
            # there) and we want density to rank them.
            exclude = (anchor_set | gapsmith_assigned) \
                        - gapsmith_for_rxn.get(target_rxn, set())
            # But never exclude the target itself:
            exclude.discard(target)
            unassigned = [(g, d) for g, d in gene_density.items()
                           if g not in exclude]
            unassigned.sort(key=lambda x: -x[1])
            rank_unassigned = 0
            for r, (g, d) in enumerate(unassigned, 1):
                if g == target:
                    rank_unassigned = r
                    break
            top1_unassigned_gene = unassigned[0][0] if unassigned else ''
            target_in_gapsmith = int(target in gapsmith_assigned)

            fout.write('\t'.join([
                str(i), target, target_rxn, target_go, str(max_pid),
                str(len(nbrs)), str(len(nbr_go_weights)),
                str(len(anchor_positions)),
                f'{tden:.4f}', str(rank_all), str(rank_nonanc),
                str(len(non_anchor)),
                str(rank_unassigned), str(len(unassigned)),
                top1_unassigned_gene, str(target_in_gapsmith),
            ]) + '\n')
            if (i + 1) % 20 == 0:
                print(f'  [{i+1}/{len(cases)}]', file=sys.stderr)
    print(f'[info] wrote {args.out}', file=sys.stderr)


if __name__ == '__main__':
    main()
