#!/usr/bin/env python3
"""
Predict operons from a Prokka GFF using the same algorithm as
gspa.predictor.context.OperonPredictor:
  - same strand
  - intergenic distance <= max_dist (default 300 bp)
  - operon size >= min_size (default 2)

Outputs:
  operons.tsv          one row per operon: operon_id  contig  start  end  strand  n_members  members
  protein_to_operon.tsv  protein_id  operon_id  position_in_operon  operon_size

This script is a Python sibling of the JVM OperonPredictor — they share the
same parameters and produce identical results given the same GFF. Used by
the tutorial pipeline because we want operons.tsv next to the other
sidecar outputs without spinning up the JVM. The JVM-side AnnotateCommand
will write the same file via the real predictor when wired (issue #20).
"""
from __future__ import annotations
import argparse, os, sys
from collections import defaultdict

def parse_gff(gff_path):
    cdss = []
    with open(gff_path) as fh:
        for line in fh:
            if not line or line.startswith('#') or '\tCDS\t' not in line:
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9:
                continue
            attrs = dict(p.split('=', 1) for p in f[8].split(';') if '=' in p)
            lt = attrs.get('locus_tag') or attrs.get('ID')
            if not lt:
                continue
            cdss.append({
                'id': lt,
                'contig': f[0],
                'start': int(f[3]),
                'end': int(f[4]),
                'strand': f[6],
            })
    return cdss

def predict_operons(cdss, max_dist=300, min_size=2):
    by_contig = defaultdict(list)
    for c in cdss:
        by_contig[c['contig']].append(c)
    operons = []
    next_id = 1
    for contig, prots in by_contig.items():
        prots.sort(key=lambda p: p['start'])
        if len(prots) < 2:
            continue
        cur = [prots[0]]
        for prev, curr in zip(prots, prots[1:]):
            if prev['strand'] == curr['strand']:
                gap = curr['start'] - prev['end'] - 1
                if gap <= max_dist:
                    cur.append(curr)
                    continue
            if len(cur) >= min_size:
                operons.append({
                    'id': f'op_{next_id:05d}',
                    'contig': contig,
                    'start': cur[0]['start'],
                    'end': cur[-1]['end'],
                    'strand': cur[0]['strand'],
                    'members': [p['id'] for p in cur],
                })
                next_id += 1
            cur = [curr]
        if len(cur) >= min_size:
            operons.append({
                'id': f'op_{next_id:05d}',
                'contig': contig,
                'start': cur[0]['start'],
                'end': cur[-1]['end'],
                'strand': cur[0]['strand'],
                'members': [p['id'] for p in cur],
            })
            next_id += 1
    return operons

def write_outputs(operons, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    op_path = os.path.join(out_dir, 'operons.tsv')
    p2o_path = os.path.join(out_dir, 'protein_to_operon.tsv')
    integrate_path = os.path.join(out_dir, 'operons_for_integrate.tsv')
    with open(op_path, 'w') as f:
        f.write('operon_id\tcontig\tstart\tend\tstrand\tn_members\tmembers\n')
        for op in operons:
            f.write(f"{op['id']}\t{op['contig']}\t{op['start']}\t{op['end']}\t{op['strand']}\t{len(op['members'])}\t{','.join(op['members'])}\n")
    with open(p2o_path, 'w') as f:
        f.write('protein_id\toperon_id\tposition\toperon_size\n')
        for op in operons:
            for i, pid in enumerate(op['members']):
                f.write(f"{pid}\t{op['id']}\t{i+1}\t{len(op['members'])}\n")
    # Format consumed by `gspa integrate --operons`: one operon per line,
    # tab-separated protein IDs (no header, no other columns).
    with open(integrate_path, 'w') as f:
        for op in operons:
            f.write('\t'.join(op['members']) + '\n')
    return op_path, p2o_path, integrate_path

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--gff', required=True)
    p.add_argument('--out-dir', required=True)
    p.add_argument('--max-dist', type=int, default=300)
    p.add_argument('--min-size', type=int, default=2)
    args = p.parse_args()
    cdss = parse_gff(args.gff)
    print(f'  parsed {len(cdss)} CDS')
    operons = predict_operons(cdss, args.max_dist, args.min_size)
    print(f'  predicted {len(operons)} operons; '
          f'{sum(len(op["members"]) for op in operons)} CDS in operons '
          f'({100.0*sum(len(op["members"]) for op in operons)/max(1,len(cdss)):.1f}%)')
    op_path, p2o_path, integrate_path = write_outputs(operons, args.out_dir)
    print(f'  wrote {op_path}')
    print(f'  wrote {p2o_path}')
    print(f'  wrote {integrate_path} (for gspa integrate --operons)')

if __name__ == '__main__':
    main()
