#!/usr/bin/env python3
"""Analyze leave-reaction-out results.

Reads results.tsv from run_ablation.py (and optionally the test-cases
TSV from build_testset.py for stratification columns), computes:

  - overall hit@1, hit@3, hit@5
  - mean reciprocal rank (excluding cases with rank=0, and a variant
    including them as 1/inf=0)
  - mean candidate-set size
  - hit@1 margin distribution
  - hit@1 stratified by:
      - n_neighbors_local (bucketed: 2, 3-4, >=5)
      - orthogroup_size (bucketed: 1, 2-3)

Outputs a plain-text report to stdout and optionally a CSV.
"""
import argparse
import collections
import statistics
import sys


def read_tsv(path):
    rows = []
    with open(path) as f:
        header = f.readline().rstrip('\n').split('\t')
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < len(header):
                parts += [''] * (len(header) - len(parts))
            rows.append(dict(zip(header, parts)))
    return rows


def compute(rows, label=None):
    n = len(rows)
    if n == 0:
        return {}
    hit1 = sum(1 for r in rows if r.get('rank_of_p') == '1')
    hit3 = sum(1 for r in rows if r.get('rank_of_p') not in ('', '0') and int(r['rank_of_p']) <= 3)
    hit5 = sum(1 for r in rows if r.get('rank_of_p') not in ('', '0') and int(r['rank_of_p']) <= 5)
    hit_any = sum(1 for r in rows if r.get('rank_of_p') not in ('', '0'))
    mrr_all = 0.0
    for r in rows:
        rk = r.get('rank_of_p', '0')
        try:
            rki = int(rk)
        except ValueError:
            rki = 0
        mrr_all += (1.0 / rki) if rki > 0 else 0.0
    mrr_all /= n
    # Mean |cand|
    ncs = [int(r['n_candidates']) for r in rows if r.get('n_candidates', '').isdigit()]
    mean_nc = statistics.mean(ncs) if ncs else 0
    # Margins for hit@1
    margins = [float(r['margin']) for r in rows
               if r.get('rank_of_p') == '1' and r.get('margin', '') not in ('', 'nan')]
    mean_margin = statistics.mean(margins) if margins else 0
    return {
        'label': label or 'all',
        'n': n,
        'hit@1': hit1 / n,
        'hit@3': hit3 / n,
        'hit@5': hit5 / n,
        'recall (any rank)': hit_any / n,
        'MRR': mrr_all,
        'mean |cand|': mean_nc,
        'mean margin (hit@1)': mean_margin,
    }


def fmt(m):
    return (
        f"  n={m['n']:4d}  "
        f"hit@1={m['hit@1']:.3f}  hit@3={m['hit@3']:.3f}  hit@5={m['hit@5']:.3f}  "
        f"recall={m['recall (any rank)']:.3f}  MRR={m['MRR']:.3f}  "
        f"|cand|={m['mean |cand|']:.1f}  margin={m['mean margin (hit@1)']:.3f}"
    )


def bucket_neighbors(v):
    try:
        n = int(v)
    except (TypeError, ValueError):
        return 'unk'
    if n <= 2:
        return '2'
    if n <= 4:
        return '3-4'
    return '>=5'


def bucket_ortho(v):
    try:
        n = int(v)
    except (TypeError, ValueError):
        return 'unk'
    if n == 1:
        return '1'
    return '2-3'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--results', required=True)
    ap.add_argument('--group-by', default=None,
                    help='Column name to stratify on (e.g. n_neighbors_local)')
    args = ap.parse_args()

    rows = read_tsv(args.results)
    print(f'read {len(rows)} rows from {args.results}')
    if not rows:
        sys.exit(0)

    print('\n== Overall ==')
    print(fmt(compute(rows)))

    print('\n== By ablation mode ==')
    by_mode = collections.defaultdict(list)
    for r in rows:
        by_mode[r.get('mode', '?')].append(r)
    for mode, rs in sorted(by_mode.items()):
        print(f'[mode={mode}]')
        print(fmt(compute(rs, label=mode)))

    print('\n== By n_neighbors_local ==')
    by_nb = collections.defaultdict(list)
    for r in rows:
        by_nb[bucket_neighbors(r.get('n_neighbors_local'))].append(r)
    for k in ['2', '3-4', '>=5']:
        if not by_nb.get(k):
            continue
        print(f'[neighbors={k}]')
        print(fmt(compute(by_nb[k], label=k)))

    print('\n== By orthogroup_size ==')
    by_o = collections.defaultdict(list)
    for r in rows:
        by_o[bucket_ortho(r.get('orthogroup_size'))].append(r)
    for k in ['1', '2-3']:
        if not by_o.get(k):
            continue
        print(f'[ortho={k}]')
        print(fmt(compute(by_o[k], label=k)))

    # Report a small sample of successes and failures
    print('\n== 10 sample successes (hit@1) ==')
    succ = [r for r in rows if r.get('rank_of_p') == '1'][:10]
    for r in succ:
        print(f"  {r['protein_id']:15s}  {r['reaction_id']:20s}  nbrs={r.get('n_neighbors_local','')}  ortho={r.get('orthogroup_size','')}  margin={r.get('margin','')}")

    print('\n== 10 sample misses (rank=0 or not in cand list) ==')
    miss = [r for r in rows if r.get('rank_of_p') in ('', '0')][:10]
    for r in miss:
        print(f"  {r['protein_id']:15s}  {r['reaction_id']:20s}  nbrs={r.get('n_neighbors_local','')}  ortho={r.get('orthogroup_size','')}  top={r.get('top_candidate','')}")


if __name__ == '__main__':
    main()
