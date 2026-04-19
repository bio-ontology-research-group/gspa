#!/usr/bin/env python3
"""Leave-reaction-out test-set builder.

Reads gapsmith output, ground truth, operons, orthogroup map, and GFF,
emits a TSV of (protein, reaction, ec, go, pathway, ...) tuples where
the reaction is a good candidate for LRO evaluation:

  - gapsmith status=good_blast for this (rxn, protein)
  - reaction has EC -> GO mapping
  - reaction has >=2 neighbor reactions in the same gapsmith pathway
    catalyzed by other proteins (good_blast), AND those neighbor
    proteins lie within a 20 kb genomic window of the target protein
  - target protein's orthogroup has <=3 members (inside this genome)
  - ground truth contains (protein, GO) -- i.e. the ablation target
    is corroborated by GOA

Output columns:
  protein_id  reaction_id  pathway_id  ec  go_term  contig  start  end
  strand  operon_idx  n_neighbors_in_pathway  n_neighbors_local
  neighbor_proteins  orthogroup_size  truth_ok
"""
import argparse
import collections
import csv
import json
import os
import re
import sys


def parse_ec2go(path):
    """Return {ec_4digit: [go_term, ...]}. Only keep 4-digit ECs."""
    ec2go = collections.defaultdict(list)
    with open(path) as f:
        for line in f:
            if line.startswith('!') or not line.strip():
                continue
            # Format: "EC:1.1.1.1 > GO:alcohol dehydrogenase (NAD+) activity ; GO:0004022"
            m = re.match(r'^EC:(\S+)\s+>\s+GO:.*;\s*(GO:\d+)\s*$', line.strip())
            if not m:
                continue
            ec, go = m.group(1), m.group(2)
            if '-' in ec:
                continue
            ec2go[ec].append(go)
    return dict(ec2go)


def parse_reactions_tbl(path):
    """Yield dict-rows from gapsmith *-all-Reactions.tbl."""
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            yield row


def parse_gff_cds(gff_path, refseq_to_uniprot):
    """Return {uniprot_acc: {contig, start, end, strand}} from GFF CDS rows."""
    coords = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9 or parts[2] != 'CDS':
                continue
            contig, start, end, strand = parts[0], int(parts[3]), int(parts[4]), parts[6]
            attrs = parts[8]
            m = re.search(r'protein_id=([^;]+)', attrs)
            if not m:
                continue
            refseq = m.group(1)
            uniprot = refseq_to_uniprot.get(refseq)
            if not uniprot:
                continue
            # Some proteins have multiple CDS entries (multi-exon -- rare in
            # prokaryotes but possible). Keep the outermost span.
            if uniprot in coords:
                c = coords[uniprot]
                if c['contig'] != contig:
                    continue
                c['start'] = min(c['start'], start)
                c['end'] = max(c['end'], end)
            else:
                coords[uniprot] = {'contig': contig, 'start': start, 'end': end, 'strand': strand}
    return coords


def parse_map_tsv(path):
    """mg1655_map.tsv: refseq_id <TAB> uniprot_acc."""
    out = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                out[parts[0]] = parts[1]
    return out


def parse_operons(path):
    """Return list-of-lists and a {protein: idx} index."""
    operons = []
    idx = {}
    with open(path) as f:
        for i, line in enumerate(f):
            members = line.strip().split('\t')
            members = [m for m in members if m]
            if not members:
                continue
            operons.append(members)
            for m in members:
                idx[m] = len(operons) - 1
    return operons, idx


def parse_orthogroups(path, genome_members=None):
    """Count orthogroup size restricted to genome_members (for paralog count).

    orthogroup_map.tsv rows: member <TAB> rep_namespace:rep_protein
    We want, per protein, the size of its orthogroup restricted to the
    target-genome's own members (within-genome paralog count), not the
    cross-strain cluster size.
    """
    protein_to_rep = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            member, rep = parts[0], parts[1]
            if ':' in rep:
                rep_genome, rep_acc = rep.split(':', 1)
            else:
                rep_genome, rep_acc = 'unknown', rep
            protein_to_rep[member] = (rep_genome, rep_acc)
    rep_to_local = collections.defaultdict(set)
    for member, rep in protein_to_rep.items():
        if genome_members is not None and member not in genome_members:
            continue
        rep_to_local[rep].add(member)
    size_by_protein = {}
    for member, rep in protein_to_rep.items():
        if genome_members is not None and member not in genome_members:
            continue
        size_by_protein[member] = len(rep_to_local[rep])
    return size_by_protein


def parse_truth(path):
    """Return {(protein, go): aspect}."""
    out = {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            acc, aspect, go = parts[0], parts[1], parts[2]
            out[(acc, go)] = aspect
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--reactions-tbl', required=True)
    ap.add_argument('--gff', required=True)
    ap.add_argument('--map-tsv', required=True)
    ap.add_argument('--operons', required=True)
    ap.add_argument('--ortho-map', required=True)
    ap.add_argument('--truth', required=True)
    ap.add_argument('--ec2go', required=True)
    ap.add_argument('--genome-tag', default='mg1655')
    ap.add_argument('--window-kb', type=int, default=20)
    ap.add_argument('--max-ortho-size', type=int, default=3)
    ap.add_argument('--min-neighbors-local', type=int, default=2)
    ap.add_argument('--out', required=True)
    ap.add_argument('--verbose', action='store_true')
    args = ap.parse_args()

    print(f'[info] loading EC2GO from {args.ec2go}', file=sys.stderr)
    ec2go = parse_ec2go(args.ec2go)

    print(f'[info] loading RefSeq->UniProt map from {args.map_tsv}', file=sys.stderr)
    refseq2uniprot = parse_map_tsv(args.map_tsv)
    print(f'         {len(refseq2uniprot)} entries', file=sys.stderr)

    print(f'[info] loading GFF coords from {args.gff}', file=sys.stderr)
    coords = parse_gff_cds(args.gff, refseq2uniprot)
    print(f'         {len(coords)} proteins with coords', file=sys.stderr)

    print(f'[info] loading operons from {args.operons}', file=sys.stderr)
    operons, protein_to_operon = parse_operons(args.operons)
    print(f'         {len(operons)} operons, {len(protein_to_operon)} proteins', file=sys.stderr)

    print(f'[info] loading orthogroups from {args.ortho_map}', file=sys.stderr)
    # Restrict to members present in the mg1655 map (as UniProt acc)
    genome_members = set(refseq2uniprot.values())
    ortho_size = parse_orthogroups(args.ortho_map, genome_members=genome_members)
    print(f'         {len(ortho_size)} proteins with orthogroup', file=sys.stderr)

    print(f'[info] loading truth from {args.truth}', file=sys.stderr)
    truth = parse_truth(args.truth)
    print(f'         {len(truth)} (protein, GO) pairs', file=sys.stderr)

    # Step 1: parse Reactions.tbl, index good_blast hits by pathway.
    # We store (pathway -> [(rxn, protein, ec, stitle, pident), ...]).
    # A given rxn can have multiple good_blast hits (multi-subunit or paralogs);
    # we keep the top-scoring by pident.
    print(f'[info] parsing gapsmith Reactions.tbl', file=sys.stderr)
    pathway_good = collections.defaultdict(list)
    best_by_rxn = {}  # (pathway, rxn) -> (protein, ec, pident)
    for row in parse_reactions_tbl(args.reactions_tbl):
        if row['status'] != 'good_blast':
            continue
        pathway = row['pathway']
        rxn = row['rxn']
        ec = row['ec'].strip()
        protein = row['stitle'].strip()
        try:
            pident = float(row['pident'])
        except (TypeError, ValueError):
            pident = 0.0
        if not protein:
            continue
        key = (pathway, rxn)
        # keep best (highest pident) if multiple candidates
        cur = best_by_rxn.get(key)
        if cur is None or pident > cur[2]:
            best_by_rxn[key] = (protein, ec, pident)

    # Now index by pathway for neighbor lookups.
    pathway_rxns = collections.defaultdict(dict)  # pathway -> {rxn: (protein, ec)}
    for (pathway, rxn), (protein, ec, _) in best_by_rxn.items():
        pathway_rxns[pathway][rxn] = (protein, ec)

    # Step 2: iterate over candidate targets and apply filters.
    window_bp = args.window_kb * 1000
    out_rows = []
    skipped = collections.Counter()
    for (pathway, rxn), (protein, ec, pident) in best_by_rxn.items():
        # filter 1: has EC
        if not ec:
            skipped['no_ec'] += 1
            continue
        # filter 2: has EC -> GO mapping (prefer exact 4-digit)
        if ec not in ec2go:
            skipped['no_ec2go'] += 1
            continue
        gos = ec2go[ec]
        if not gos:
            skipped['no_ec2go'] += 1
            continue
        # filter 3: has coords
        if protein not in coords:
            skipped['no_coords'] += 1
            continue
        pcoord = coords[protein]
        # filter 4: orthogroup size
        osize = ortho_size.get(protein, 1)
        if osize > args.max_ortho_size:
            skipped['ortho_too_large'] += 1
            continue
        # filter 5: neighbor reactions in the same gapsmith pathway catalyzed
        # by OTHER proteins, with coords within window
        neighbor_proteins = []
        n_in_pathway = 0
        for nrxn, (nprot, nec) in pathway_rxns[pathway].items():
            if nrxn == rxn:
                continue
            if nprot == protein:
                continue
            n_in_pathway += 1
            nc = coords.get(nprot)
            if not nc:
                continue
            if nc['contig'] != pcoord['contig']:
                continue
            if abs(nc['start'] - pcoord['start']) > window_bp:
                continue
            neighbor_proteins.append(nprot)
        if len(neighbor_proteins) < args.min_neighbors_local:
            skipped['too_few_local_neighbors'] += 1
            continue
        # filter 6: truth_all has at least one of protein's target GO terms
        truth_ok = any((protein, g) in truth for g in gos)
        if not truth_ok:
            skipped['not_in_truth'] += 1
            continue
        # Pick the single best GO: the one in truth, else the most specific.
        go = next((g for g in gos if (protein, g) in truth), gos[0])
        operon_idx = protein_to_operon.get(protein, -1)
        out_rows.append({
            'protein_id': protein,
            'reaction_id': rxn,
            'pathway_id': f'MetaCyc:{pathway}',
            'ec': ec,
            'go_term': go,
            'contig': pcoord['contig'],
            'start': pcoord['start'],
            'end': pcoord['end'],
            'strand': pcoord['strand'],
            'operon_idx': operon_idx,
            'n_neighbors_in_pathway': n_in_pathway,
            'n_neighbors_local': len(neighbor_proteins),
            'neighbor_proteins': ','.join(neighbor_proteins),
            'orthogroup_size': osize,
            'truth_ok': int(truth_ok),
        })
    print(f'[info] {len(out_rows)} test cases after filtering', file=sys.stderr)
    print(f'[info] rejection counts:', file=sys.stderr)
    for k, v in skipped.most_common():
        print(f'         {k}: {v}', file=sys.stderr)

    with open(args.out, 'w') as f:
        if not out_rows:
            f.write('')
            return
        cols = list(out_rows[0].keys())
        f.write('\t'.join(cols) + '\n')
        for r in out_rows:
            f.write('\t'.join(str(r[c]) for c in cols) + '\n')
    print(f'[info] wrote {args.out}', file=sys.stderr)


if __name__ == '__main__':
    main()
