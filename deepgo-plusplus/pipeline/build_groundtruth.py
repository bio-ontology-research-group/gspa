#!/usr/bin/env python3
"""Build a faithful CAFA6 ground-truth from a dated GOA filter (Phase-0).

Input: `testsuperset_exp_annots.tsv` produced on IBEX by filtering the current
goa_uniprot_all.gaf.gz to test-superset proteins + experimental evidence codes:

    accession \\t GO_id \\t aspect(F|P|C) \\t evidence \\t date(YYYYMMDD)

Applies the official CAFA6 protocol:
  * t0 = 2026-02-02 (final submission deadline).
  * A protein is a TEST target for aspect a if it accumulated >=1 experimental
    annotation in a with date > t0.
  * Knowledge class is per-protein from the pre-t0 state across the 3 aspects:
      no-knowledge      = no exp annotation in ANY aspect on/before t0
      limited-knowledge = exp annotation in 1-2 aspects on/before t0
      partial-knowledge = exp annotation in all 3 aspects on/before t0
  * "Full evaluation mode": the ground truth for a test (protein, aspect) is the
    protein's COMPLETE experimental annotation set in that aspect (any date) --
    leaf terms; CAFA-evaluator propagates ancestors itself.

Writes CAFA-evaluator ground-truth files (protein \\t GO):
    gt_all.tsv, gt_no.tsv, gt_limited.tsv, gt_partial.tsv
and a `test_proteins.txt` (all test targets) for pre-filtering predictions.
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict

ASP = {'F': 'MF', 'P': 'BP', 'C': 'CC'}
T0 = 20260202


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--annots', required=True,
                    help='testsuperset_exp_annots.tsv (acc, GO, F/P/C, evidence, date)')
    ap.add_argument('--t0', type=int, default=T0)
    ap.add_argument('--outdir', default='.')
    args = ap.parse_args()

    # per (protein, aspect): all terms; and earliest/any pre-t0 presence
    terms = defaultdict(set)              # (prot, aspect) -> {GO}
    pre = defaultdict(set)                # prot -> {aspects with date<=t0}
    post = defaultdict(set)              # prot -> {aspects with date>t0}
    n = 0
    with open(args.annots) as fh:
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) < 5:
                continue
            acc, go, asp, _ev, date = p[0], p[1], ASP.get(p[2]), p[3], p[4]
            if asp is None:
                continue
            try:
                d = int(date)
            except ValueError:
                continue
            n += 1
            terms[(acc, asp)].add(go)
            (pre if d <= args.t0 else post)[acc].add(asp)

    # knowledge class per protein (from pre-t0 aspects)
    def kclass(prot):
        k = len(pre.get(prot, ()))
        return 'no' if k == 0 else ('partial' if k == 3 else 'limited')

    # test targets: (prot, aspect) where aspect in post[prot]
    gt = {'all': [], 'no': [], 'limited': [], 'partial': []}
    test_proteins = set()
    n_targets = defaultdict(int)
    for prot, asps in post.items():
        kc = kclass(prot)
        for asp in asps:
            rows = [(prot, go) for go in terms[(prot, asp)]]
            gt['all'].extend(rows)
            gt[kc].extend(rows)
            test_proteins.add(prot)
            n_targets[(kc, asp)] += 1

    import os
    os.makedirs(args.outdir, exist_ok=True)
    for k, rows in gt.items():
        with open(os.path.join(args.outdir, f'gt_{k}.tsv'), 'w') as out:
            for prot, go in rows:
                out.write(f'{prot}\t{go}\n')
    with open(os.path.join(args.outdir, 'test_proteins.txt'), 'w') as out:
        for p in sorted(test_proteins):
            out.write(p + '\n')

    print(f'annotations read           : {n:,}')
    print(f'test target proteins (uniq): {len(test_proteins):,}')
    print('test targets by class/aspect:')
    for kc in ('no', 'limited', 'partial'):
        row = '  '.join(f'{a}={n_targets.get((kc, a), 0)}' for a in ('MF', 'BP', 'CC'))
        print(f'  {kc:<8} {row}')
    print('wrote gt_all.tsv gt_no.tsv gt_limited.tsv gt_partial.tsv test_proteins.txt '
          f'-> {args.outdir}', file=sys.stderr)


if __name__ == '__main__':
    main()
