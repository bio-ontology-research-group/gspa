#!/usr/bin/env python3
"""Part 2: cross-genome non-anchor co-occurrence catalog.

For each (orthogroup C, reaction R) pair, count across panel genomes G':
  n_sig_with_nonanchor_C = #genomes where
      (a) R has anchor-density peaks in G'
      (b) some C-member in G' sits IN a density peak
      (c) that C-member is NOT itself an anchor for R's neighbours

  n_sig_total = #genomes with property (a)
  n_base_with_C = #genomes where C has any member
  n_base_total = panel size

LR(C, R) = (n_sig_nonanchor / n_sig_total) / (n_base_with_C / n_base_total)

This differs from M2's catalog (build_catalog.py) by restricting the
"sig_with" count to NON-ANCHOR C-members — aligning with the user's
"[dark] in an operon of A/B/C/D anchors" formulation.
"""
import argparse
import collections
import json
import math
import os
import re
import sys
from pathlib import Path


def read_orthogroups(path):
    out = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            member, rep = parts[0], parts[1]
            if ':' in rep:
                rep = rep.split(':', 1)[1]
            out[member] = rep
    return out


def load_ec2go(path):
    out = {}
    with open(path) as f:
        for line in f:
            if line.startswith('!') or not line.strip():
                continue
            m = re.match(r'^EC:(\S+)\s*>\s*GO:[^;]+;\s*(GO:\d+)', line)
            if m:
                out[m.group(1)] = m.group(2)
                out[m.group(1)[3:]] = m.group(2)
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


def load_layout(path):
    out = {}
    by_contig = collections.defaultdict(list)
    with open(path) as f:
        f.readline()
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 5:
                continue
            pid, contig, start, end, _ = (parts[0], parts[1],
                                            int(parts[2]), int(parts[3]),
                                            parts[4])
            mid = (start + end) // 2
            out[pid] = (contig, mid)
            by_contig[contig].append((mid, pid))
    for c in by_contig:
        by_contig[c].sort()
    return out, dict(by_contig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--manifest', required=True)
    ap.add_argument('--root', required=True)
    ap.add_argument('--orthogroup-map', required=True)
    ap.add_argument('--reactions-tsv', required=True)
    ap.add_argument('--diffusion-tsv', required=True)
    ap.add_argument('--ec-aliases-tsv', required=True)
    ap.add_argument('--ec2go', required=True)
    ap.add_argument('--exclude-tag', default=None)
    ap.add_argument('--tau', type=float, default=0.3)
    ap.add_argument('--bandwidth-bp', type=int, default=5000)
    ap.add_argument('--radius-k', type=int, default=2)
    ap.add_argument('--alpha', type=float, default=0.5)
    ap.add_argument('--restrict-to-gaps-file', default=None)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    exclude_tags = set((args.exclude_tag or '').split(',')) \
        if args.exclude_tag else set()
    tags = []
    with open(args.manifest) as f:
        f.readline()
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if parts and parts[0] and parts[0] not in exclude_tags:
                tags.append(parts[0])
    print(f'[info] panel: {len(tags)} genomes', file=sys.stderr)

    ortho = read_orthogroups(args.orthogroup_map)
    ec2go = load_ec2go(args.ec2go)
    rxn_to_mets, currency, met2rxn, rxn_to_ecs, ec_to_rxns = \
        load_reaction_graph(args.reactions_tsv, args.diffusion_tsv,
                             args.ec_aliases_tsv)

    # Pre-load per-genome (posteriors, layout) + baseline orthogroup
    # membership.
    per_genome = {}
    base_presence = collections.defaultdict(set)  # og -> {tags}
    # ortho keys can be plain (<pid>) or tag-prefixed (<tag>:<pid>) —
    # support both so the catalog builds regardless of which FASTA
    # prefixing convention the clustering used.
    for tag in tags:
        intg = Path(args.root) / 'integrated' / f'{tag}_integrated.tsv'
        lay = Path(args.root) / 'layout' / f'{tag}_layout.tsv'
        if not intg.exists() or not lay.exists():
            continue
        posts = load_integrated(str(intg))
        layout, by_contig = load_layout(str(lay))
        per_genome[tag] = (posts, layout, by_contig)
        for pid in layout:
            og = ortho.get(pid) or ortho.get(f'{tag}:{pid}')
            if og:
                base_presence[og].add(tag)
    panel_size = len(per_genome)
    print(f'[info] per-genome data loaded for {panel_size} tags',
          file=sys.stderr)

    # Reaction universe
    target_rxns = None
    if args.restrict_to_gaps_file and os.path.exists(args.restrict_to_gaps_file):
        target_rxns = set()
        with open(args.restrict_to_gaps_file) as f:
            for line in f:
                try:
                    obj = json.loads(line)
                except json.JSONDecodeError:
                    continue
                ec = obj.get('ec_number')
                if ec:
                    for r in ec_to_rxns.get(ec, []):
                        target_rxns.add(r)
                rxn = obj.get('reaction_id')
                if rxn and rxn in rxn_to_mets:
                    target_rxns.add(rxn)
        print(f'[info] restricted to {len(target_rxns)} rxns',
              file=sys.stderr)
    rxn_iter = target_rxns if target_rxns else list(rxn_to_mets.keys())

    bandwidth = args.bandwidth_bp
    halfwidth = 4 * bandwidth

    # Per-(og, rxn) counts
    sig_with_nonanc = collections.Counter()
    sig_total = collections.Counter()

    done = 0
    for rxn_id in rxn_iter:
        done += 1
        if done % 200 == 0:
            print(f'  rxn {done}/{len(rxn_iter)}', file=sys.stderr)
        nbrs = bfs_neighbors(rxn_id, rxn_to_mets, met2rxn, currency,
                              args.radius_k, args.alpha)
        if not nbrs:
            continue
        nbr_go_weights = {}
        for n_rxn, w in nbrs.items():
            for ec in rxn_to_ecs.get(n_rxn, []):
                go = ec2go.get(f'EC:{ec}') or ec2go.get(ec)
                if not go:
                    continue
                if nbr_go_weights.get(go, 0) < w:
                    nbr_go_weights[go] = w
        if not nbr_go_weights:
            continue

        for tag, (posts, layout, by_contig) in per_genome.items():
            # Anchors in this genome
            anchor_set = set()
            anchor_positions_by_contig = collections.defaultdict(list)
            for prot, gos in posts.items():
                best = 0.0
                is_anc = False
                for go, w_r in nbr_go_weights.items():
                    p = gos.get(go, 0.0)
                    if p >= args.tau:
                        is_anc = True
                        contrib = w_r * p
                        if contrib > best:
                            best = contrib
                if is_anc:
                    anchor_set.add(prot)
                    loc = layout.get(prot)
                    if loc:
                        anchor_positions_by_contig[loc[0]].append(
                            (loc[1], best))
            if not anchor_set:
                continue
            sig_total[rxn_id] += 1

            # For each non-anchor protein in this genome with a layout,
            # check if its position falls within halfwidth of any anchor.
            # Collect which orthogroups have such members.
            ogs_with_nonanc = set()
            for contig, anchors in anchor_positions_by_contig.items():
                anchors.sort()
                for gmid, gpid in by_contig.get(contig, []):
                    if gpid in anchor_set:
                        continue
                    og = ortho.get(gpid) or ortho.get(f'{tag}:{gpid}')
                    if not og:
                        continue
                    # check any anchor within halfwidth of gmid
                    for amid, _aw in anchors:
                        if amid - gmid > halfwidth:
                            break
                        if gmid - amid > halfwidth:
                            continue
                        ogs_with_nonanc.add(og)
                        break
            for og in ogs_with_nonanc:
                sig_with_nonanc[(og, rxn_id)] += 1

    # Write TSV
    print(f'[info] writing {args.out}', file=sys.stderr)
    n_rows = 0
    with open(args.out, 'w') as f:
        f.write('# NonAnchorCatalog\n')
        f.write(f'# panel_size={panel_size}\n')
        f.write('orthogroup_id\treaction_id\tn_sig_nonanc\tn_sig_total\t'
                'n_base_with\tn_base_total\n')
        for (og, rxn), n in sig_with_nonanc.items():
            tot = sig_total.get(rxn, 0)
            if tot == 0:
                continue
            n_base = len(base_presence.get(og, set()))
            f.write(f'{og}\t{rxn}\t{n}\t{tot}\t{n_base}\t{panel_size}\n')
            n_rows += 1
    print(f'[info] wrote {n_rows} rows to {args.out}', file=sys.stderr)


if __name__ == '__main__':
    main()
