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

RELIABILITY = {'distance': 0.85, 'strict': 0.95, 'functional': 0.70}
PAIR_THRESHOLD = 0.50

# GO BP terms too broad to count as functional evidence.
BP_DENYLIST = {
    'GO:0008150', 'GO:0008152', 'GO:0009987', 'GO:0050896', 'GO:0065007',
    'GO:0050789', 'GO:0050794', 'GO:0044238', 'GO:0044237', 'GO:0071704',
    'GO:0006807', 'GO:0019222',
}

def _intergenic(a, b):
    if a['end'] < b['start']: return b['start'] - a['end'] - 1
    if b['end'] < a['start']: return a['start'] - b['end'] - 1
    return -(min(a['end'], b['end']) - max(a['start'], b['start']) + 1)

def _call_pair(a, b, bp_terms_per_protein, distance_max=300, strict_max=50, functional_max=1000):
    """Return (support_set, posterior, gap)."""
    if a['strand'] != b['strand']:
        return set(), 0.0, _intergenic(a, b)
    gap = _intergenic(a, b)
    support = set()
    if gap <= distance_max: support.add('distance')
    if gap <= strict_max:   support.add('strict')
    if gap <= functional_max:
        ta = bp_terms_per_protein.get(a['id'], set())
        tb = bp_terms_per_protein.get(b['id'], set())
        if ta and tb and (ta & tb):
            support.add('functional')
    one_minus = 1.0
    for k in support: one_minus *= (1.0 - RELIABILITY[k])
    return support, 1.0 - one_minus, gap

def predict_operons(cdss, max_dist=300, min_size=2, ensemble=False, bp_terms_per_protein=None):
    """Detect operons.

    With ensemble=False (default) this is the original distance + same-strand
    rule used by gspa.predictor.context.OperonPredictor — still emitted so the
    output stays backwards compatible with `gspa integrate --operons`.

    With ensemble=True a small collection of predictors votes on each adjacent
    gene-pair and per-pair Noisy-OR posteriors decide membership. Each operon
    additionally carries `support_set`, `min_pair_posterior`, `mean_pair_posterior`.
    """
    by_contig = defaultdict(list)
    for c in cdss:
        by_contig[c['contig']].append(c)
    bp_terms_per_protein = bp_terms_per_protein or {}
    operons = []
    next_id = 1
    for contig, prots in by_contig.items():
        prots.sort(key=lambda p: p['start'])
        if len(prots) < 2:
            continue
        cur = [prots[0]]
        cur_edges = []  # list of (support, posterior) for ensemble mode
        cur_support = set()
        for prev, curr in zip(prots, prots[1:]):
            if ensemble:
                support, posterior, gap = _call_pair(prev, curr, bp_terms_per_protein,
                                                     distance_max=max_dist)
                add_to_op = posterior >= PAIR_THRESHOLD
            else:
                add_to_op = (prev['strand'] == curr['strand']
                             and (curr['start'] - prev['end'] - 1) <= max_dist)
                support, posterior = set(), 0.0
            if add_to_op:
                cur.append(curr)
                if ensemble:
                    cur_edges.append((support, posterior))
                    cur_support |= support
                continue
            if len(cur) >= min_size:
                operons.append(_make_op(next_id, contig, cur, cur_edges, cur_support, ensemble))
                next_id += 1
            cur = [curr]
            cur_edges = []
            cur_support = set()
        if len(cur) >= min_size:
            operons.append(_make_op(next_id, contig, cur, cur_edges, cur_support, ensemble))
            next_id += 1
    return operons

def _make_op(op_idx, contig, members, edges, support, ensemble):
    op = {
        'id': f'op_{op_idx:05d}',
        'contig': contig,
        'start': members[0]['start'],
        'end': members[-1]['end'],
        'strand': members[0]['strand'],
        'members': [p['id'] for p in members],
    }
    if ensemble:
        if edges:
            posts = [e[1] for e in edges]
            op['support_set'] = sorted(support)
            op['min_pair_posterior'] = round(min(posts), 4)
            op['mean_pair_posterior'] = round(sum(posts)/len(posts), 4)
        else:
            op['support_set'] = []
            op['min_pair_posterior'] = 0.0
            op['mean_pair_posterior'] = 0.0
    return op

def write_outputs(operons, out_dir, ensemble=False):
    os.makedirs(out_dir, exist_ok=True)
    op_path = os.path.join(out_dir, 'operons.tsv')
    p2o_path = os.path.join(out_dir, 'protein_to_operon.tsv')
    integrate_path = os.path.join(out_dir, 'operons_for_integrate.tsv')
    with open(op_path, 'w') as f:
        if ensemble:
            f.write('operon_id\tcontig\tstart\tend\tstrand\tn_members\tsupport_set\tmin_pair_posterior\tmean_pair_posterior\tmembers\n')
            for op in operons:
                f.write(f"{op['id']}\t{op['contig']}\t{op['start']}\t{op['end']}\t{op['strand']}\t{len(op['members'])}\t{','.join(op.get('support_set', []))}\t{op.get('min_pair_posterior', 0):.4f}\t{op.get('mean_pair_posterior', 0):.4f}\t{','.join(op['members'])}\n")
        else:
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

def _load_bp_terms_from_gaf(gaf_path, posterior_threshold=0.5):
    """Load BP GO terms per protein from a GAF — used for the functional
    sub-predictor in ensemble mode. Tolerates missing/malformed files."""
    bp = defaultdict(set)
    if not gaf_path or not os.path.exists(gaf_path):
        return bp
    with open(gaf_path) as fh:
        for line in fh:
            if line.startswith('!') or '\t' not in line: continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9: continue
            pid, term, aspect = f[1], f[4], f[8]
            if aspect != 'P': continue
            if term in BP_DENYLIST: continue
            bp[pid].add(term)
    return bp

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--gff', required=True)
    p.add_argument('--out-dir', required=True)
    p.add_argument('--max-dist', type=int, default=300)
    p.add_argument('--min-size', type=int, default=2)
    p.add_argument('--ensemble', action='store_true',
                   help='Use the 3-predictor ensemble (distance + strict + functional) '
                        'with Noisy-OR per-pair posteriors. Operons gain support_set + '
                        'min/mean pair posterior columns.')
    p.add_argument('--gaf', help='GAF file for the functional sub-predictor (--ensemble only). '
                                 'BP-aspect annotations only. Optional.')
    args = p.parse_args()
    cdss = parse_gff(args.gff)
    print(f'  parsed {len(cdss)} CDS')
    bp_terms = _load_bp_terms_from_gaf(args.gaf) if args.ensemble else {}
    if args.ensemble and bp_terms:
        print(f'  loaded BP terms for {len(bp_terms)} proteins from {args.gaf}')
    operons = predict_operons(cdss, args.max_dist, args.min_size,
                              ensemble=args.ensemble, bp_terms_per_protein=bp_terms)
    print(f'  predicted {len(operons)} operons '
          f'({"ensemble" if args.ensemble else "distance-only"}); '
          f'{sum(len(op["members"]) for op in operons)} CDS in operons '
          f'({100.0*sum(len(op["members"]) for op in operons)/max(1,len(cdss)):.1f}%)')
    op_path, p2o_path, integrate_path = write_outputs(operons, args.out_dir, ensemble=args.ensemble)
    print(f'  wrote {op_path}')
    print(f'  wrote {p2o_path}')
    print(f'  wrote {integrate_path} (for gspa integrate --operons)')

if __name__ == '__main__':
    main()
