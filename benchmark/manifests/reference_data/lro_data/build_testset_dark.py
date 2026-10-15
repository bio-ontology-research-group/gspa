#!/usr/bin/env python3
"""Dark-matter LRO test set builder (v2).

Starts from gapsmith Reactions.tbl (pathway, rxn, EC) entries — any
status — and, for each, finds mg1655 proteins whose GOA truth asserts
the EC's GO term. This lets the target be a *bad_blast* / *no_blast*
protein (dark) as well as *good_blast* (bright), giving us the full
homology-depth spectrum.

Context check: the candidate case needs ≥2 gapsmith good_blast anchors
in the same pathway within 20 kb of the target — these are the RLC
anchors that let genomic context have a chance.

Dark matter: max pident (from claims JSONL) to any Swiss-Prot entry
asserting the same GO on this target ≤ --dark-pident.
"""
import argparse
import collections
import csv
import json
import re
import sys


def parse_ec2go(path):
    ec2go = collections.defaultdict(list)
    with open(path) as f:
        for line in f:
            if line.startswith('!') or not line.strip():
                continue
            m = re.match(r'^EC:(\S+)\s+>\s+GO:.*;\s*(GO:\d+)\s*$',
                          line.strip())
            if not m:
                continue
            ec, go = m.group(1), m.group(2)
            if '-' in ec:
                continue
            ec2go[ec].append(go)
    return dict(ec2go)


def parse_reactions_tbl(path):
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            yield row


def parse_gff_cds(gff_path, refseq_to_uniprot):
    coords = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9 or parts[2] != 'CDS':
                continue
            contig, start, end, strand = (parts[0], int(parts[3]),
                                            int(parts[4]), parts[6])
            m = re.search(r'protein_id=([^;]+)', parts[8])
            if not m:
                continue
            refseq = m.group(1)
            uniprot = refseq_to_uniprot.get(refseq)
            if not uniprot:
                continue
            if uniprot in coords:
                c = coords[uniprot]
                if c['contig'] != contig:
                    continue
                c['start'] = min(c['start'], start)
                c['end'] = max(c['end'], end)
            else:
                coords[uniprot] = {'contig': contig, 'start': start,
                                     'end': end, 'strand': strand}
    return coords


def parse_map_tsv(path):
    out = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                out[parts[0]] = parts[1]
    return out


def parse_truth(path):
    """Return {go: set(protein)} for MF aspect."""
    go_to_prots = collections.defaultdict(set)
    prot_to_gos = collections.defaultdict(set)
    with open(path) as f:
        f.readline()
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3 or parts[1] != 'MF':
                continue
            go_to_prots[parts[2]].add(parts[0])
            prot_to_gos[parts[0]].add(parts[2])
    return dict(go_to_prots), dict(prot_to_gos)


def parse_claims_pidents(path):
    """Max pident per (protein, GO) across Swiss-Prot DIAMOND hits
    asserting that GO. If (protein, GO) not in claims -> 0."""
    out = {}
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
            meta = obj.get('metadata', {}) or {}
            try:
                ident = float(meta.get('pident', 0.0))
            except (TypeError, ValueError):
                ident = 0.0
            key = (pid, go)
            cur = out.get(key, 0.0)
            if ident > cur:
                out[key] = ident
    return out


def parse_orthogroups(path, genome_members):
    protein_to_rep = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            member, rep = parts[0], parts[1]
            if ':' in rep:
                rep = rep.split(':', 1)[1]
            protein_to_rep[member] = rep
    rep_to_local = collections.defaultdict(set)
    for m, r in protein_to_rep.items():
        if m in genome_members:
            rep_to_local[r].add(m)
    return {m: len(rep_to_local[r])
             for m, r in protein_to_rep.items() if m in genome_members}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--reactions-tbl', required=True)
    ap.add_argument('--gff', required=True)
    ap.add_argument('--map-tsv', required=True)
    ap.add_argument('--ortho-map', required=True)
    ap.add_argument('--truth', required=True)
    ap.add_argument('--claims', required=True)
    ap.add_argument('--ec2go', required=True)
    ap.add_argument('--genome-tag', default='mg1655')
    ap.add_argument('--window-kb', type=int, default=20)
    ap.add_argument('--min-neighbors-local', type=int, default=2)
    ap.add_argument('--dark-pident', type=float, default=30.0)
    ap.add_argument('--out', required=True)
    ap.add_argument('--out-dark', default=None)
    args = ap.parse_args()

    ec2go = parse_ec2go(args.ec2go)
    refseq2uniprot = parse_map_tsv(args.map_tsv)
    coords = parse_gff_cds(args.gff, refseq2uniprot)
    genome_members = set(refseq2uniprot.values())
    ortho_size = parse_orthogroups(args.ortho_map, genome_members)
    go_to_prots, prot_to_gos = parse_truth(args.truth)
    pident_lookup = parse_claims_pidents(args.claims)

    print(f'[info] ECs in ec2go: {len(ec2go)}; truth MF proteins: '
          f'{len(prot_to_gos)}; claims pidents: {len(pident_lookup)}',
          file=sys.stderr)

    # Parse Reactions.tbl. Collect, per pathway, the good_blast anchors
    # and per (pathway, rxn) the EC.
    pathway_anchors = collections.defaultdict(dict)  # pathway -> {rxn: prot}
    pathway_rxn_ec = {}  # (pathway, rxn) -> ec
    rxn_ec = {}  # rxn -> ec (pick any)
    for row in parse_reactions_tbl(args.reactions_tbl):
        pathway = row['pathway']
        rxn = row['rxn']
        ec = (row.get('ec') or '').strip()
        if ec:
            pathway_rxn_ec[(pathway, rxn)] = ec
            rxn_ec[rxn] = ec
        if row.get('status') == 'good_blast':
            prot = (row.get('stitle') or '').strip()
            if prot:
                pathway_anchors[pathway][rxn] = prot

    print(f'[info] {len(pathway_anchors)} pathways have anchors; '
          f'{len(pathway_rxn_ec)} (pathway, rxn) have EC',
          file=sys.stderr)

    # For each (pathway, rxn) with EC, map EC -> GO(s), GO -> truth-protein(s)
    window_bp = args.window_kb * 1000
    out_rows = []
    skipped = collections.Counter()
    seen = set()
    for (pathway, rxn), ec in pathway_rxn_ec.items():
        gos = ec2go.get(ec, [])
        if not gos:
            skipped['no_ec2go'] += 1
            continue
        # candidate proteins (mg1655-truth-annotated with these GOs)
        cand_prots = set()
        for go in gos:
            cand_prots.update(go_to_prots.get(go, set()))
        cand_prots &= set(coords.keys())
        if not cand_prots:
            skipped['no_truth_protein'] += 1
            continue
        anchors = pathway_anchors.get(pathway, {})
        for protein in cand_prots:
            pcoord = coords[protein]
            # Find in-genome neighbours within window
            nbrs = []
            for nrxn, nprot in anchors.items():
                if nrxn == rxn:
                    continue
                if nprot == protein:
                    continue
                nc = coords.get(nprot)
                if not nc or nc['contig'] != pcoord['contig']:
                    continue
                if abs(nc['start'] - pcoord['start']) > window_bp:
                    continue
                nbrs.append(nprot)
            if len(nbrs) < args.min_neighbors_local:
                skipped['too_few_neighbors'] += 1
                continue
            key = (protein, rxn, pathway)
            if key in seen:
                continue
            seen.add(key)
            # Pick the best matching GO (one that truth has AND ec2go has)
            matched_go = None
            for g in gos:
                if protein in go_to_prots.get(g, set()):
                    matched_go = g
                    break
            if not matched_go:
                continue
            max_pid = pident_lookup.get((protein, matched_go), 0.0)
            out_rows.append({
                'protein_id': protein,
                'reaction_id': rxn,
                'pathway_id': f'MetaCyc:{pathway}',
                'ec': ec,
                'go_term': matched_go,
                'contig': pcoord['contig'],
                'start': pcoord['start'],
                'end': pcoord['end'],
                'strand': pcoord['strand'],
                'operon_idx': -1,
                'n_neighbors_in_pathway': len(anchors) - 1,
                'n_neighbors_local': len(nbrs),
                'neighbor_proteins': ','.join(nbrs),
                'orthogroup_size': ortho_size.get(protein, 1),
                'truth_ok': 1,
                'max_pident_annot_homolog': f'{max_pid:.2f}',
            })

    print(f'[info] {len(out_rows)} candidate cases; skips: {dict(skipped)}',
          file=sys.stderr)

    cols = ['protein_id', 'reaction_id', 'pathway_id', 'ec', 'go_term',
             'contig', 'start', 'end', 'strand', 'operon_idx',
             'n_neighbors_in_pathway', 'n_neighbors_local',
             'neighbor_proteins', 'orthogroup_size', 'truth_ok',
             'max_pident_annot_homolog']
    with open(args.out, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in out_rows:
            f.write('\t'.join(str(r[c]) for c in cols) + '\n')
    print(f'[info] wrote {args.out}', file=sys.stderr)

    # Dark subset + histogram
    bins = [(0, 0), (0.01, 20), (20, 30), (30, 40), (40, 60),
             (60, 80), (80, 101)]
    for lo, hi in bins:
        n = sum(1 for r in out_rows
                 if lo <= float(r['max_pident_annot_homolog']) <= hi)
        print(f'  pident {lo:5.1f}-{hi:5.1f}: {n}', file=sys.stderr)
    dark = [r for r in out_rows
             if float(r['max_pident_annot_homolog']) <= args.dark_pident]
    print(f'[info] dark subset (pident ≤ {args.dark_pident}): '
          f'{len(dark)} / {len(out_rows)}', file=sys.stderr)
    if args.out_dark:
        with open(args.out_dark, 'w') as f:
            f.write('\t'.join(cols) + '\n')
            for r in dark:
                f.write('\t'.join(str(r[c]) for c in cols) + '\n')
        print(f'[info] wrote {args.out_dark}', file=sys.stderr)


if __name__ == '__main__':
    main()
