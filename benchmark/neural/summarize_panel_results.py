#!/usr/bin/env python3
"""Consolidate per-(genome, predictor, truth) JSON eval rows into a human
summary. Reads eval_all.jsonl (concatenated indent=2 JSON objects) and
emits:

1. A wide TSV with one row per (genome, predictor, truth).
2. A Markdown summary table aggregating per predictor (mean ± SD across
   the 21 genomes) for fmax_overall / fmax_cafa_overall / smin_overall,
   plus per-aspect numbers when present.
"""
from __future__ import annotations

import argparse
import json
import statistics
from collections import defaultdict
from pathlib import Path


def read_concat_json(path: Path) -> list[dict]:
    txt = path.read_text()
    decoder = json.JSONDecoder()
    out = []
    i = 0
    n = len(txt)
    while i < n:
        while i < n and txt[i] in ' \n\t\r':
            i += 1
        if i >= n:
            break
        obj, end = decoder.raw_decode(txt[i:])
        out.append(obj)
        i += end
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--in', dest='inp', required=True)
    ap.add_argument('--tsv', required=True)
    ap.add_argument('--md', required=True)
    args = ap.parse_args()

    rows = read_concat_json(Path(args.inp))
    if not rows:
        print('no rows'); return

    # Collect all column names preserving order of first appearance.
    cols = []
    seen = set()
    for r in rows:
        for k in r:
            if k not in seen:
                cols.append(k); seen.add(k)

    def fmt_cell(v):
        if isinstance(v, float):
            return f'{v:.4f}'
        if isinstance(v, list):
            return '[' + ','.join(f'{x:.4f}' if isinstance(x, float) else str(x) for x in v) + ']'
        return '' if v is None else str(v)

    with Path(args.tsv).open('w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for r in rows:
            fh.write('\t'.join(fmt_cell(r.get(c)) for c in cols) + '\n')

    # Group by (predictor, truth, annotation_type)
    grp = defaultdict(list)
    for r in rows:
        key = (r.get('method') or r.get('predictor', '?'),
               r.get('truth', '?'),
               r.get('annotation_type', '?'))
        grp[key].append(r)

    def agg(vals):
        vals = [v for v in vals if isinstance(v, (int, float))]
        if not vals:
            return '—', '—', 0
        m = statistics.mean(vals)
        sd = statistics.stdev(vals) if len(vals) > 1 else 0.0
        return f'{m:.3f}', f'{sd:.3f}', len(vals)

    MD_METRICS = ['fmax_overall', 'fmax_cafa_overall', 'smin_overall',
                  'fmax_MF', 'fmax_BP', 'fmax_CC',
                  'ec_level1_fmax', 'ec_level2_fmax',
                  'ec_level3_fmax', 'ec_level4_fmax']

    with Path(args.md).open('w') as fh:
        fh.write('# Panel evaluation summary\n\n')
        for (pred, truth, ann), items in sorted(grp.items()):
            fh.write(f'## predictor={pred}  truth={truth}  type={ann}\n\n')
            cols_md = ['metric', 'mean', 'sd', 'n']
            fh.write('| ' + ' | '.join(cols_md) + ' |\n')
            fh.write('|' + '|'.join(['---'] * len(cols_md)) + '|\n')
            for m in MD_METRICS:
                vals = [i.get(m) for i in items if m in i]
                mean, sd, n = agg(vals)
                if n > 0:
                    fh.write(f'| {m} | {mean} | {sd} | {n} |\n')
            fh.write('\n')

    print('wrote', args.tsv, args.md)


if __name__ == '__main__':
    main()
