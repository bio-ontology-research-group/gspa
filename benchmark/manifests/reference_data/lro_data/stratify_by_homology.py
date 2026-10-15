#!/usr/bin/env python3
"""Stratify LRO result files by panel-homology depth.

For each case (target_protein, target_go):
  max_pident = max over panel-proteins-annotated-with-target_go: pident(target, panel)
(If no such panel protein has that GO in truth, max_pident = 0 → fully dark.)

Buckets:
  dark      < 30
  twilight  30-50
  bright    >= 50

Reports hit@1 per bucket per result file.
"""
import argparse
import collections
import sys
from pathlib import Path


def load_truth_dir(truth_dir, panel_tags):
    """Return {(tag, protein): set(GO)}."""
    out = collections.defaultdict(set)
    for tag in panel_tags:
        p = Path(truth_dir) / f'{tag}_truth_all.tsv'
        if not p.exists():
            continue
        with open(p) as f:
            f.readline()
            for line in f:
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 3:
                    continue
                out[(tag, parts[0])].add(parts[2])
    return dict(out)


def load_diamond(path):
    """Return {query: [(subject, pident, tag, sub_prot)]}."""
    out = collections.defaultdict(list)
    with open(path) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            q, s, pid = parts[0], parts[1], float(parts[2])
            # panel subject format: "tag:accession"
            if ':' in s:
                tag, sub = s.split(':', 1)
            else:
                tag, sub = '', s
            out[q].append((s, pid, tag, sub))
    return dict(out)


def load_panel_tags(manifest, exclude=None):
    exclude = set(exclude or [])
    tags = []
    with open(manifest) as f:
        f.readline()
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if parts and parts[0] and parts[0] not in exclude:
                tags.append(parts[0])
    return tags


def load_cases(path):
    rows = []
    with open(path) as f:
        header = f.readline().rstrip('\n').split('\t')
        for line in f:
            parts = line.rstrip('\n').split('\t')
            rows.append(dict(zip(header, parts)))
    return rows


def load_results(path):
    rows = []
    with open(path) as f:
        header = f.readline().rstrip('\n').split('\t')
        for line in f:
            parts = line.rstrip('\n').split('\t')
            rows.append(dict(zip(header, parts)))
    return rows


def bucket_of(pid):
    if pid <= 0.0:
        return 'no_annot_homolog'
    if pid < 30.0:
        return 'dark_<30'
    if pid < 50.0:
        return 'twilight_30_50'
    if pid < 70.0:
        return 'near_50_70'
    if pid < 90.0:
        return 'close_70_90'
    return 'very_close_90+'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--manifest', required=True)
    ap.add_argument('--truth-dir', required=True)
    ap.add_argument('--diamond-tsv', required=True,
                    help='DIAMOND blastp output: q s pident')
    ap.add_argument('--cases', required=True)
    ap.add_argument('--results', nargs='+', required=True)
    ap.add_argument('--exclude-tag', default='mg1655')
    ap.add_argument('--out', default=None)
    args = ap.parse_args()

    tags = load_panel_tags(args.manifest,
                            exclude=args.exclude_tag.split(',')
                            if args.exclude_tag else None)
    truth = load_truth_dir(args.truth_dir, tags)
    print(f'[info] truth entries: {sum(len(v) for v in truth.values())} '
          f'across {len(truth)} (tag, protein) pairs', file=sys.stderr)

    diamond = load_diamond(args.diamond_tsv)
    print(f'[info] diamond queries: {len(diamond)}', file=sys.stderr)

    cases = load_cases(args.cases)
    # Compute per-case max_pident over panel proteins annotated with target GO
    # + max_pident over ANY panel protein (regardless of annotation).
    case_depth = {}
    case_anyhit = {}
    for case in cases:
        prot = case['protein_id']
        go = case['go_term']
        hits = diamond.get(prot, [])
        best_annot = 0.0
        best_any = 0.0
        for (subj, pid, tag, sub_prot) in hits:
            if not tag or tag == args.exclude_tag:
                continue
            if pid > best_any:
                best_any = pid
            if go in truth.get((tag, sub_prot), set()):
                if pid > best_annot:
                    best_annot = pid
        case_depth[prot] = best_annot
        case_anyhit[prot] = best_any
    depths = list(case_depth.values())
    anyhits = list(case_anyhit.values())
    depths.sort()
    anyhits.sort()
    print(f'[info] annot-depth pctl (min/25/50/75/max): '
          f'{depths[0]:.1f} / {depths[len(depths)//4]:.1f} / '
          f'{depths[len(depths)//2]:.1f} / {depths[3*len(depths)//4]:.1f} / '
          f'{depths[-1]:.1f}', file=sys.stderr)
    print(f'[info] any-hit    pctl (min/25/50/75/max): '
          f'{anyhits[0]:.1f} / {anyhits[len(anyhits)//4]:.1f} / '
          f'{anyhits[len(anyhits)//2]:.1f} / '
          f'{anyhits[3*len(anyhits)//4]:.1f} / {anyhits[-1]:.1f}',
          file=sys.stderr)
    # Cross-tabulate: annot_depth bucket × any_hit bucket
    crosstab = collections.Counter()
    for p in case_depth:
        ab = bucket_of(case_depth[p])
        hb = bucket_of(case_anyhit[p])
        crosstab[(ab, hb)] += 1
    print('\n[info] crosstab: annotated-homolog pident × any-hit pident',
          file=sys.stderr)
    for (ab, hb), n in sorted(crosstab.items()):
        print(f'  annot={ab:20s}  any={hb:20s}  n={n}', file=sys.stderr)

    # Compute per-bucket hit@k for each results file
    out_lines = []
    header = (f'{"variant":40s}  {"bucket":10s}  {"n":>4}  '
              f'{"hit@1":>6}  {"hit@3":>6}  {"hit@5":>6}  {"MRR":>6}')
    print(header)
    out_lines.append(header)
    for rp in args.results:
        rows = load_results(rp)
        bucket_rows = collections.defaultdict(list)
        for r in rows:
            prot = r.get('protein_id', '')
            d = case_depth.get(prot)
            if d is None:
                continue
            bucket = bucket_of(d)
            bucket_rows[bucket].append(r)
            bucket_rows['ALL'].append(r)
        for bucket in ['ALL', 'no_annot_homolog', 'dark_<30',
                       'twilight_30_50', 'near_50_70',
                       'close_70_90', 'very_close_90+']:
            rs = bucket_rows.get(bucket, [])
            if not rs:
                continue
            n = len(rs)
            ranks = [int(r.get('rank_of_p', '0') or 0) for r in rs]
            hit1 = sum(1 for x in ranks if x == 1) / n
            hit3 = sum(1 for x in ranks if 1 <= x <= 3) / n
            hit5 = sum(1 for x in ranks if 1 <= x <= 5) / n
            mrr = sum(1.0 / x for x in ranks if x > 0) / n
            name = Path(rp).stem.replace('results_', '')
            line = (f'{name:40s}  {bucket:10s}  {n:4d}  '
                    f'{hit1:6.3f}  {hit3:6.3f}  {hit5:6.3f}  {mrr:6.3f}')
            print(line)
            out_lines.append(line)
        print('')
        out_lines.append('')

    if args.out:
        with open(args.out, 'w') as f:
            f.write('\n'.join(out_lines) + '\n')
        print(f'[info] saved to {args.out}', file=sys.stderr)


if __name__ == '__main__':
    main()
