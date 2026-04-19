#!/usr/bin/env python3
"""Phase 12 M2: build a ReactionLocusCatalog TSV from a panel of
genomes + an orthogroup map.

For each panel genome g and each reaction R, we:
  1. Compute anchor proteins (proteins with posterior > τ on GO terms
     catalysing any r' ∈ N_k(R)).
  2. Compute R-signature-window membership: proteins whose locus falls
     in any window around an anchor (halfWidth = 4 * bandwidth).
  3. If any anchor exists for R in g, g is "R-signature-present".

Then for each (orthogroup C, reaction R) we aggregate:
  n_sig_with = #genomes where R-signature is present AND any C-member
               locus lies in any R-signature window
  n_sig_total = #genomes where R-signature is present
  n_base_with = #genomes where C has any member
  n_base_total = total panel size

Output TSV (same header as ReactionLocusCatalog.writeTo):
  # panel_size=N
  orthogroup_id  reaction_id  n_sig_with  n_sig_total  n_base_with  n_base_total
"""
import argparse
import collections
import gzip
import json
import math
import os
import re
import sys


# ------------------------------------------------------------------
# IO helpers
# ------------------------------------------------------------------

def read_tsv(path, want_header=True):
    with open(path) as f:
        if want_header:
            header = f.readline().rstrip('\n').split('\t')
        else:
            header = None
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            yield header, parts


def read_layout(path):
    """Return {proteinId: (contig, start, end, strand)}."""
    out = {}
    for header, parts in read_tsv(path):
        if len(parts) < 5:
            continue
        out[parts[0]] = (parts[1], int(parts[2]), int(parts[3]), parts[4])
    return out


def read_integrated(path, tau=0.3):
    """Return {goTerm: {proteinId: prob}} filtered to prob > tau.

    integrated.tsv from gspa-cli has columns:
      protein_id type function_id go_aspect posterior_prob ...
    """
    out = collections.defaultdict(dict)
    for header, parts in read_tsv(path):
        if len(parts) < 5 or parts[1] != 'GO':
            continue
        try:
            prob = float(parts[4])
        except ValueError:
            continue
        if prob <= tau:
            continue
        out[parts[2]][parts[0]] = prob
    return dict(out)


def read_claims_as_posteriors(path, tau=0.3):
    """Fallback when integrated.tsv is unavailable: use raw claim
    scores as an approximation. Returns {goTerm: {proteinId: prob}}."""
    out = collections.defaultdict(dict)
    with open(path) as f:
        for line in f:
            try:
                obj = json.loads(line)
            except json.JSONDecodeError:
                continue
            if obj.get('function_type') != 'GO':
                continue
            pid = obj.get('protein_id')
            go = obj.get('function_id')
            score = float(obj.get('raw_score', 0.0))
            if score < tau:
                continue
            cur = out[go].get(pid, 0.0)
            if score > cur:
                out[go][pid] = score
    return dict(out)


def read_orthogroup_map(path):
    """Return {proteinId: orthogroupId}.

    Format varies:
      col1=member col2=rep       (bench_gtdb30's orthogroup_map_50.tsv)
      col1=member col2=ns:rep    (bench_ecoli's orthogroup_map.tsv)
    We take col2 as the orthogroup ID, stripping any 'prefix:' namespace.
    """
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


# ------------------------------------------------------------------
# Reaction graph (minimal port of the Groovy loader)
# ------------------------------------------------------------------

def load_reaction_graph(reactions_tsv, diffusion_tsv, ec_aliases_tsv,
                       currency_percentile=99.0):
    """Return (reactions, currency_mets, ec_to_rxns, rxn_to_ecs,
               metabolite_to_reactions)."""
    currency = set()
    if diffusion_tsv and os.path.exists(diffusion_tsv):
        with open(diffusion_tsv) as f:
            f.readline()
            for line in f:
                parts = line.split('\t')
                if parts and parts[0].startswith('cpd'):
                    currency.add(parts[0].strip())

    rxn_to_mets = {}  # rxnId -> set of cpd ids (all touched)
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

    # Degree percentile
    if currency_percentile < 100.0 and raw_deg:
        sorted_deg = sorted(raw_deg.values())
        idx = max(0, min(len(sorted_deg) - 1,
                         int(math.ceil(len(sorted_deg) * currency_percentile / 100.0)) - 1))
        threshold = sorted_deg[idx]
        for m, d in raw_deg.items():
            if d >= threshold:
                currency.add(m)

    # Metabolite -> reactions (post-currency)
    metabolite_to_reactions = collections.defaultdict(set)
    for rxn, mets in rxn_to_mets.items():
        for m in mets:
            if m not in currency:
                metabolite_to_reactions[m].add(rxn)

    # EC aliases
    ec_to_rxns = collections.defaultdict(set)
    rxn_to_ecs = collections.defaultdict(set)
    if ec_aliases_tsv and os.path.exists(ec_aliases_tsv):
        with open(ec_aliases_tsv) as f:
            header = f.readline()
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
                        ec_to_rxns[ec].add(rxn)
                        rxn_to_ecs[rxn].add(ec)

    return (rxn_to_mets, currency, dict(ec_to_rxns), dict(rxn_to_ecs),
            dict(metabolite_to_reactions))


def bfs_neighbors(rxn_id, rxn_to_mets, metabolite_to_reactions,
                  currency, max_k, alpha):
    """Return {rxnId: weight} for 1..max_k hops."""
    out = {}
    seen = {rxn_id}
    frontier = [rxn_id]
    weight = 1.0
    for d in range(1, max_k + 1):
        weight *= alpha
        nxt = []
        for f in frontier:
            for m in rxn_to_mets.get(f, set()):
                if m in currency:
                    continue
                for n in metabolite_to_reactions.get(m, set()):
                    if n not in seen:
                        seen.add(n)
                        out[n] = weight
                        nxt.append(n)
        if not nxt:
            break
        frontier = nxt
    return out


# ------------------------------------------------------------------
# ec2go
# ------------------------------------------------------------------

def load_ec2go(path):
    out = {}
    with open(path) as f:
        for line in f:
            if line.startswith('!') or not line.strip():
                continue
            m = re.match(r'^(EC:\S+)\s*>\s*GO:[^;]+;\s*(GO:\d+)', line)
            if m:
                out[m.group(1)] = m.group(2)
                # Also store bare form
                out[m.group(1)[3:]] = m.group(2)
    return out


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------

def anchor_set(genome_posteriors, nbr_go_weights, tau):
    """Return {proteinId: weight}."""
    out = {}
    for go, w in nbr_go_weights.items():
        per_prot = genome_posteriors.get(go)
        if not per_prot:
            continue
        for pid, prob in per_prot.items():
            if prob < tau:
                continue
            contrib = w * prob
            cur = out.get(pid)
            if cur is None or contrib > cur:
                out[pid] = contrib
    return out


def window_proteins(anchors, layout, half_width):
    """Return set of protein IDs in any window around anchor proteins.
    A "window" is ±half_width bp on the same contig."""
    # Group layout by contig, sorted by midpoint
    by_contig = collections.defaultdict(list)
    for pid, (contig, start, end, strand) in layout.items():
        mid = (start + end) // 2
        by_contig[contig].append((mid, pid))
    for c in by_contig:
        by_contig[c].sort()

    out = set()
    for aid in anchors:
        loc = layout.get(aid)
        if not loc:
            continue
        contig, start, end, strand = loc
        center = (start + end) // 2
        lo = center - half_width
        hi = center + half_width
        for mid, pid in by_contig.get(contig, []):
            if mid < lo:
                continue
            if mid > hi:
                break
            out.add(pid)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--manifest', required=True,
                    help='TSV with header + tag column (first col).')
    ap.add_argument('--root', required=True,
                    help='Panel root (contains integrated/, layout/, gaps/)')
    ap.add_argument('--orthogroup-map', required=True)
    ap.add_argument('--reactions-tsv', required=True)
    ap.add_argument('--diffusion-tsv', required=True)
    ap.add_argument('--ec-aliases-tsv', required=True)
    ap.add_argument('--ec2go', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--radius-k', type=int, default=2)
    ap.add_argument('--alpha', type=float, default=0.5)
    ap.add_argument('--tau', type=float, default=0.3,
                    help='min posterior to be an anchor')
    ap.add_argument('--bandwidth-bp', type=int, default=5000,
                    help='kernel bandwidth (4σ halfwidth used as window)')
    ap.add_argument('--exclude-tag', default=None,
                    help='Tag(s) to hold out of the catalog (comma-sep).')
    ap.add_argument('--use-claims-as-posteriors', action='store_true',
                    help='Fallback when integrated.tsv is unavailable.')
    ap.add_argument('--reactions-per-gap', type=int, default=0,
                    help='Limit the reaction universe to those in any '
                         'panel genome gaps.jsonl (0 = all graph reactions).')
    ap.add_argument('--restrict-to-gaps-file', default=None,
                    help='gaps.jsonl file — restrict reaction universe '
                         'to the EC-equivalents of reactions listed here.')
    args = ap.parse_args()

    exclude = set((args.exclude_tag or '').split(',')) if args.exclude_tag else set()

    # Load panel manifest
    tags = []
    with open(args.manifest) as f:
        f.readline()
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if parts and parts[0] and parts[0] not in exclude:
                tags.append(parts[0])
    print(f'[info] panel size: {len(tags)} (excluded {len(exclude)})', file=sys.stderr)

    ortho = read_orthogroup_map(args.orthogroup_map)
    print(f'[info] orthogroup map: {len(ortho)} members', file=sys.stderr)

    ec2go = load_ec2go(args.ec2go)
    print(f'[info] ec2go: {len(ec2go)} entries', file=sys.stderr)

    (rxn_to_mets, currency, ec_to_rxns, rxn_to_ecs,
     metabolite_to_reactions) = load_reaction_graph(
        args.reactions_tsv, args.diffusion_tsv, args.ec_aliases_tsv)
    print(f'[info] reaction graph: {len(rxn_to_mets)} rxns, '
          f'{len(currency)} currency, {len(ec_to_rxns)} ECs',
          file=sys.stderr)

    # Optional: constrain to gap-reaction-equivalent set
    target_reactions = None
    if args.restrict_to_gaps_file and os.path.exists(args.restrict_to_gaps_file):
        target_reactions = set()
        with open(args.restrict_to_gaps_file) as f:
            for line in f:
                try:
                    obj = json.loads(line)
                except json.JSONDecodeError:
                    continue
                ec = obj.get('ec_number')
                if ec:
                    for r in ec_to_rxns.get(ec, []):
                        target_reactions.add(r)
                rxn = obj.get('reaction_id')
                if rxn and rxn in rxn_to_mets:
                    target_reactions.add(rxn)
        print(f'[info] gap-file-restricted reaction universe: {len(target_reactions)}',
              file=sys.stderr)
    elif args.reactions_per_gap > 0:
        target_reactions = set()
        for tag in tags:
            gaps_path = os.path.join(args.root, 'gaps', f'{tag}_gaps.jsonl')
            if not os.path.exists(gaps_path):
                continue
            with open(gaps_path) as f:
                for line in f:
                    try:
                        obj = json.loads(line)
                    except json.JSONDecodeError:
                        continue
                    ec = obj.get('ec_number')
                    if ec:
                        for r in ec_to_rxns.get(ec, []):
                            target_reactions.add(r)
        print(f'[info] gap-derived reaction universe: {len(target_reactions)}',
              file=sys.stderr)

    # Stage 1: per-genome posteriors + layout + per-genome orthogroup membership
    per_genome = {}
    for tag in tags:
        layout_path = os.path.join(args.root, 'layout', f'{tag}_layout.tsv')
        if not os.path.exists(layout_path):
            print(f'[warn] missing layout for {tag}; skipping', file=sys.stderr)
            continue
        if args.use_claims_as_posteriors:
            claims_path = os.path.join(args.root, 'claims', f'{tag}_dp_claims.jsonl')
            if not os.path.exists(claims_path):
                print(f'[warn] missing claims for {tag}; skipping', file=sys.stderr)
                continue
            posts = read_claims_as_posteriors(claims_path, args.tau)
        else:
            intg_path = os.path.join(args.root, 'integrated', f'{tag}_integrated.tsv')
            if not os.path.exists(intg_path):
                print(f'[warn] missing integrated.tsv for {tag}; skipping',
                      file=sys.stderr)
                continue
            posts = read_integrated(intg_path, args.tau)
        layout = read_layout(layout_path)
        per_genome[tag] = (posts, layout)
        print(f'[info] {tag}: {sum(len(v) for v in posts.values())} posteriors, '
              f'{len(layout)} loci', file=sys.stderr)

    panel_size = len(per_genome)

    # Stage 2: baseline presence per orthogroup per genome
    # C is "present in g" if any member of C has a layout entry in g.
    base_presence = collections.defaultdict(set)  # og -> {tags}
    for tag, (posts, layout) in per_genome.items():
        gs = set()
        for pid in layout:
            og = ortho.get(pid)
            if og:
                base_presence[og].add(tag)

    # Stage 3: per-reaction anchor windows.
    # Iterate reactions in target_reactions (or all graph reactions).
    reactions_iter = target_reactions if target_reactions else rxn_to_mets.keys()
    print(f'[info] catalog over {len(reactions_iter)} reactions × '
          f'{panel_size} genomes', file=sys.stderr)

    half_width = 4 * args.bandwidth_bp
    # (og, rxn) -> n_sig_with
    sig_with = collections.Counter()
    # rxn -> #genomes where R-signature present
    sig_total = collections.Counter()

    n_progress = 0
    for rxn_id in reactions_iter:
        n_progress += 1
        if n_progress % 200 == 0:
            print(f'  ... {n_progress}/{len(reactions_iter)} reactions',
                  file=sys.stderr)
        # Compute neighbour GO weights for this reaction.
        nbr = bfs_neighbors(rxn_id, rxn_to_mets, metabolite_to_reactions,
                           currency, args.radius_k, args.alpha)
        if not nbr:
            continue
        nbr_go_weights = {}
        for n_rxn, w in nbr.items():
            for ec in rxn_to_ecs.get(n_rxn, []):
                go = ec2go.get(f'EC:{ec}') or ec2go.get(ec)
                if not go:
                    continue
                cur = nbr_go_weights.get(go)
                if cur is None or w > cur:
                    nbr_go_weights[go] = w
        if not nbr_go_weights:
            continue
        for tag, (posts, layout) in per_genome.items():
            anchors = anchor_set(posts, nbr_go_weights, args.tau)
            if not anchors:
                continue
            sig_total[rxn_id] += 1
            window_pids = window_proteins(anchors, layout, half_width)
            # Per-og membership in window
            ogs_in_window = set()
            for pid in window_pids:
                og = ortho.get(pid)
                if og:
                    ogs_in_window.add(og)
            for og in ogs_in_window:
                sig_with[(og, rxn_id)] += 1

    # Write TSV
    print(f'[info] writing {args.out}', file=sys.stderr)
    n_rows = 0
    with open(args.out, 'w') as f:
        f.write('# ReactionLocusCatalog\n')
        f.write(f'# panel_size={panel_size}\n')
        f.write('orthogroup_id\treaction_id\tn_sig_with\tn_sig_total\t'
                'n_base_with\tn_base_total\n')
        for (og, rxn), n_sig_w in sig_with.items():
            n_sig_t = sig_total.get(rxn, 0)
            if n_sig_t == 0:
                continue
            n_base_w = len(base_presence.get(og, set()))
            n_base_t = panel_size
            f.write(f'{og}\t{rxn}\t{n_sig_w}\t{n_sig_t}\t{n_base_w}\t{n_base_t}\n')
            n_rows += 1
    print(f'[info] wrote {n_rows} rows to {args.out}', file=sys.stderr)


if __name__ == '__main__':
    main()
