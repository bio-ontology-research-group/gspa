#!/usr/bin/env python3
"""Net-KNN component for the DeepGO-PlusPlus: guilt-by-association GO
transfer over the STRING PPI network (GOAlpha's Net-KNN). For each query
protein, vote its STRING neighbours' pre-t0 GO labels, weighted by the
STRING combined-confidence score.

Temporal integrity: edges are STRING v12.0 (2023) and neighbour labels are
the pre-t0 train_terms -- both before t0 = 2026-02-02.

Inputs:
  --index        text_string_index.tsv  (accession\\ttaxon\\tstring_id\\t...)
  --train-terms  pre-t0 labels (EntryID\\tterm\\taspect)
  --queries      query protein IDs (one per line)
  --string-dir   directory of per-species STRING files named
                 <taxid>.protein.links.v12.0.txt.gz  (space-separated:
                 protein1 protein2 combined_score, with header)
Output: protein\\tterm\\tscore  (neighbour-vote, normalised per protein).
"""
from __future__ import annotations

import argparse
import gzip
import os
import sys
import time
from collections import defaultdict

import numpy as np


def log(m):
    print(f'[{time.strftime("%H:%M:%S")}] {m}', file=sys.stderr, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--index', required=True)
    ap.add_argument('--train-terms', required=True)
    ap.add_argument('--queries', required=True)
    ap.add_argument('--string-dir', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--min-conf', type=int, default=400,
                    help='STRING combined-score cutoff (default 400 = medium)')
    ap.add_argument('--topk', type=int, default=50, help='max neighbours to vote')
    ap.add_argument('--min-score', type=float, default=0.01)
    args = ap.parse_args()

    queries = {l.strip() for l in open(args.queries) if l.strip()}
    log(f'queries: {len(queries):,}')

    # accession -> (taxon, string_id) ; string_id -> accession (first wins)
    ac_taxon, ac_sid, sid_ac = {}, {}, {}
    with open(args.index) as fh:
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) < 3 or not p[2]:
                continue
            ac, tax, sid = p[0], p[1], p[2]
            ac_taxon[ac] = tax
            ac_sid[ac] = sid
            sid_ac.setdefault(sid, ac)

    # string_id -> set(GO) from pre-t0 train labels (the transferable neighbour labels)
    sid_labels = defaultdict(set)
    with open(args.train_terms) as fh:
        next(fh, None)
        for line in fh:
            row = line.rstrip('\n').split('\t')
            p, t = row[0], row[1]
            sid = ac_sid.get(p)
            if sid:
                sid_labels[sid].add(t)
    log(f'string nodes with pre-t0 labels: {len(sid_labels):,}')

    # group queries by taxon so we read each species' STRING file once
    q_by_taxon = defaultdict(list)
    for q in queries:
        tax = ac_taxon.get(q)
        sid = ac_sid.get(q)
        if tax and sid:
            q_by_taxon[tax].append((q, sid))
    log(f'taxa with query proteins: {len(q_by_taxon)}')

    n_out = 0
    with open(args.out, 'w') as out:
        for tax, qlist in sorted(q_by_taxon.items(), key=lambda kv: -len(kv[1])):
            path = os.path.join(args.string_dir, f'{tax}.protein.links.v12.0.txt.gz')
            if not os.path.exists(path):
                continue
            qsids = {sid for _q, sid in qlist}
            # adjacency for query nodes only: sid -> list[(neighbour_sid, conf)]
            adj = defaultdict(list)
            opener = gzip.open if path.endswith('.gz') else open
            try:
                with opener(path, 'rt') as fh:
                    next(fh, None)  # header
                    for line in fh:
                        a, b, sc = line.split()
                        c = int(sc)
                        if c < args.min_conf:
                            continue
                        if a in qsids:
                            adj[a].append((b, c))
                        if b in qsids:
                            adj[b].append((a, c))
            except (OSError, EOFError, ValueError) as e:
                log(f'  taxon {tax}: SKIPPED (unreadable STRING file: {e})')
                continue
            for q, sid in qlist:
                nbrs = adj.get(sid)
                if not nbrs:
                    continue
                nbrs.sort(key=lambda x: -x[1])
                vote = defaultdict(float)
                for nb_sid, conf in nbrs[:args.topk]:
                    w = conf / 1000.0
                    for t in sid_labels.get(nb_sid, ()):  # absent -> no labels
                        vote[t] += w
                if not vote:
                    continue
                mx = max(vote.values())
                for t, v in vote.items():
                    s = v / mx
                    if s >= args.min_score:
                        out.write(f'{q}\t{t}\t{s:.4f}\n')
                        n_out += 1
            log(f'  taxon {tax}: {len(qlist)} queries done ({n_out:,} rows so far)')
    log(f'done: {n_out:,} rows -> {args.out}')


if __name__ == '__main__':
    main()
