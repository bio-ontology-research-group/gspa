#!/usr/bin/env python3
"""
Build a GSPA pathway TSV from KEGG's ec_pathway link + pathway names +
ec2go mapping.

Output columns (per PathwayLoader):
  pathway_id  pathway_name  go_term  reaction_id  ec_number  depends_on

We collapse KEGG's (ec, pathway) edges into one row per (pathway, ec), drop
the 'map' (reference) pathway_id prefix in favour of the 'ec' prefix, fill
the GO term by looking up ec2go, use EC as the reaction_id, leave the
depends_on column blank (we don't have step ordering from KEGG REST).
"""
import argparse
import re


def load_ec2go(path):
    """Parse ec2go lines like 'EC:1.1.1.1 > GO:alcohol dehydrogenase activity ; GO:0004022'."""
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith('!') or line.startswith('#'):
                continue
            m = re.match(r'^(EC:\S+)\s*>\s*GO:[^;]+;\s*(GO:\d{7})', line)
            if m:
                out[m.group(1).replace('EC:', '')] = m.group(2)
    return out


def load_pathway_names(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t', 1)
            if len(parts) == 2:
                pid = parts[0].replace('path:', '').replace('map', '')
                out[pid] = parts[1]
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--ec-pathway', required=True, help='KEGG ec_pathway.list')
    p.add_argument('--pathway-names', required=True, help='KEGG list/pathway output')
    p.add_argument('--ec2go', required=True, help='ec2go file')
    p.add_argument('--out', required=True)
    p.add_argument('--skip-global', action='store_true', default=True,
                   help='Drop the global "metabolic pathways" umbrellas 01100-01130')
    args = p.parse_args()

    ec2go = load_ec2go(args.ec2go)
    print(f'EC->GO mappings: {len(ec2go):,}')
    names = load_pathway_names(args.pathway_names)
    print(f'KEGG pathway names: {len(names):,}')

    # Drop umbrellas that cover everything
    drop = {'01100', '01110', '01120', '01200', '01210', '01230', '01240', '01250'}

    rows = []
    pathway_reactions = {}
    with open(args.ec_pathway) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) != 2:
                continue
            ec = parts[0].replace('ec:', '')
            pid_raw = parts[1].replace('path:', '')
            # Prefer the 'ec' prefix so (ec00010) and (map00010) do not duplicate
            if pid_raw.startswith('ec'):
                pid_numeric = pid_raw[2:]
            elif pid_raw.startswith('map'):
                continue  # handled via the ec-prefixed duplicate
            else:
                pid_numeric = pid_raw
            if args.skip_global and pid_numeric in drop:
                continue
            go_term = ec2go.get(ec, '')
            # Reaction ID = pathway-local index (we don't have KEGG reactions linked here)
            pathway_reactions.setdefault(pid_numeric, []).append((ec, go_term))

    with open(args.out, 'w') as out:
        out.write('pathway_id\tpathway_name\tgo_term\treaction_id\tec_number\tdepends_on\n')
        total_rows = 0
        for pid_numeric, ecs in pathway_reactions.items():
            name = names.get(pid_numeric, f'pathway {pid_numeric}')
            pathway_id = f'KEGG:{pid_numeric}'
            for idx, (ec, go) in enumerate(ecs):
                rxn = f'RXN-{pid_numeric}-{idx+1}'
                out.write(f'{pathway_id}\t{name}\t{go}\t{rxn}\tEC:{ec}\t\n')
                total_rows += 1
    print(f'wrote {total_rows:,} rows, {len(pathway_reactions):,} pathways -> {args.out}')
    # Count how many rows have GO terms vs blank
    with_go = sum(1 for row in pathway_reactions.values() for (_, go) in row if go)
    print(f'rows with GO term: {with_go:,} / {total_rows:,}')


if __name__ == '__main__':
    main()
