#!/usr/bin/env python3
"""Homology-bridged Net-KNN — extend STRING guilt-by-association to proteins
that are NOT in STRING (the realistic "novel protein" case).

Plain `net` (build_net_component.py) only fires for a query that is itself a
STRING node. A genuinely novel protein usually is not. This bridges the gap with
DIAMOND: query q -> best homolog(s) h that ARE STRING nodes -> h's STRING
neighbours N(h) -> transfer their pre-t0 GO labels to q, weighted by the
homology strength (q->h) and the STRING confidence (h-neighbour).

CAFA temporal integrity (the whole point):
  * DIAMOND targets are the pre-t0 train proteins (train.fasta) -> homologs h are
    pre-t0; and STRING v12.0 (2023) + train_terms labels are pre-t0.
  * No-knowledge queries have NO pre-t0 labels, so they are absent from the train
    DB -> DIAMOND can never self-match (no leak). h != q and neighbours != q are
    also enforced.
  * q's own (post-t0) annotations are never read.

Note the species hop: neighbours of h live in *h's* species STRING file, not q's.

Inputs:
  --index        text_string_index.tsv  (accession\\ttaxon\\tstring_id\\t...)
  --train-terms  pre-t0 labels (EntryID\\tterm\\taspect)
  --diamond      DIAMOND m8 of queries vs the train DB, columns:
                 qseqid sseqid bitscore pident  (outfmt 6)
  --string-dir   per-species <taxid>.protein.links.v12.0.txt.gz
Output: protein\\tterm\\tscore  (homology-bridged neighbour vote, per-protein max-normalised)
"""
from __future__ import annotations

import argparse
import gzip
import os
import sys
import time
from collections import defaultdict


def log(m):
    print(f'[{time.strftime("%H:%M:%S")}] {m}', file=sys.stderr, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--index', required=True)
    ap.add_argument('--train-terms', required=True)
    ap.add_argument('--diamond', required=True, help='m8: qseqid sseqid bitscore pident')
    ap.add_argument('--string-dir', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--min-conf', type=int, default=400)
    ap.add_argument('--topk-homologs', type=int, default=5,
                    help='STRING-member homologs to bridge through per query')
    ap.add_argument('--topk-neighbours', type=int, default=50)
    ap.add_argument('--min-score', type=float, default=0.01)
    args = ap.parse_args()

    # accession -> (taxon, string_id)
    ac_taxon, ac_sid = {}, {}
    with open(args.index) as fh:
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) < 3 or not p[2]:
                continue
            ac_taxon[p[0]], ac_sid[p[0]] = p[1], p[2]
    log(f'string-id accessions: {len(ac_sid):,}')

    # string_id -> set(GO) from pre-t0 train labels
    sid_labels = defaultdict(set)
    with open(args.train_terms) as fh:
        next(fh, None)
        for line in fh:
            row = line.rstrip('\n').split('\t')
            sid = ac_sid.get(row[0])
            if sid:
                sid_labels[sid].add(row[1])
    log(f'string nodes with pre-t0 labels: {len(sid_labels):,}')

    # diamond hits: q -> top-k homologs h that are STRING members, with per-query
    # normalised bitscore as the homology weight.
    q_homologs = defaultdict(list)
    with open(args.diamond) as fh:
        for line in fh:
            c = line.rstrip('\n').split('\t')
            if len(c) < 3:
                continue
            q, h = c[0], c[1]
            if h == q or h not in ac_sid:
                continue
            try:
                bits = float(c[2])
            except ValueError:
                continue
            q_homologs[q].append((h, bits))
    for q in q_homologs:
        hs = sorted(q_homologs[q], key=lambda x: -x[1])[:args.topk_homologs]
        mx = hs[0][1] if hs else 1.0
        q_homologs[q] = [(h, b / mx) for h, b in hs]   # normalise to best=1
    log(f'queries with >=1 STRING-member homolog: {len(q_homologs):,}')

    # group the homolog STRING nodes we must look up, by their species file
    need_by_taxon = defaultdict(set)        # taxon -> {homolog sid}
    for q, hs in q_homologs.items():
        for h, _w in hs:
            need_by_taxon[ac_taxon[h]].add(ac_sid[h])

    # adjacency for the needed homolog nodes only: sid_h -> [(neighbour_sid, conf)]
    adj = defaultdict(list)
    for tax, want in sorted(need_by_taxon.items(), key=lambda kv: -len(kv[1])):
        path = os.path.join(args.string_dir, f'{tax}.protein.links.v12.0.txt.gz')
        if not os.path.exists(path):
            continue
        try:
            with gzip.open(path, 'rt') as fh:
                next(fh, None)
                for line in fh:
                    a, b, sc = line.split()
                    c = int(sc)
                    if c < args.min_conf:
                        continue
                    if a in want:
                        adj[a].append((b, c))
                    if b in want:
                        adj[b].append((a, c))
        except (OSError, EOFError, ValueError) as e:
            log(f'  taxon {tax}: SKIPPED ({e})')
    for s in adj:
        adj[s].sort(key=lambda x: -x[1])
    log(f'homolog nodes with neighbours: {len(adj):,}')

    n_out = 0
    with open(args.out, 'w') as out:
        for q, hs in q_homologs.items():
            qsid = ac_sid.get(q)            # may be None (set A); used only to skip self
            vote = defaultdict(float)
            for h, w_qh in hs:
                hsid = ac_sid[h]
                for nb_sid, conf in adj.get(hsid, ())[:args.topk_neighbours]:
                    if nb_sid == qsid:
                        continue
                    w = w_qh * (conf / 1000.0)
                    for t in sid_labels.get(nb_sid, ()):
                        vote[t] += w
            if not vote:
                continue
            mx = max(vote.values())
            for t, v in vote.items():
                s = v / mx
                if s >= args.min_score:
                    out.write(f'{q}\t{t}\t{s:.4f}\n')
                    n_out += 1
    log(f'done: {n_out:,} rows -> {args.out}')


if __name__ == '__main__':
    main()
