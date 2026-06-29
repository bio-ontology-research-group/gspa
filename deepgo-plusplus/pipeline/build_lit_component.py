#!/usr/bin/env python3
"""Literature component for the DeepGO-PlusPlus: a BM25 text-kNN that
transfers GO terms from textually-similar *training* proteins to query
proteins (the spirit of GOAlpha's GORetriever: retrieve-and-transfer over
text, minus the cross-encoder rerank).

Temporal integrity (the key choice -- see build_text_string_index.py):
  * CORPUS  = training proteins (pre-t0 train_terms labels), document text =
              name + characterization (full), since their labels are pre-t0.
  * QUERY   = test proteins, query text = *name fields only* (identification,
              temporally safe), NEVER the post-t0 CC FUNCTION characterization,
              so a novel protein cannot retrieve its own later-discovered
              function. Query proteins are also excluded from the corpus.

Output: protein\\tterm\\tscore (same component format the LTR integrator
loads), score in (0,1] = BM25-weighted neighbour vote, normalised per protein.
"""
from __future__ import annotations

import argparse
import re
import sys
import time
from collections import defaultdict

import numpy as np
from rank_bm25 import BM25Okapi

_TOK = re.compile(r'[a-z0-9]+')
_STOP = {
    'the', 'a', 'an', 'of', 'and', 'or', 'to', 'in', 'is', 'are', 'for', 'by',
    'with', 'as', 'at', 'be', 'this', 'that', 'from', 'protein', 'putative',
    'probable', 'uncharacterized', 'fragment', 'isoform', 'family', 'member',
}


def log(m):
    print(f'[{time.strftime("%H:%M:%S")}] {m}', file=sys.stderr, flush=True)


def tok(s: str) -> list[str]:
    return [t for t in _TOK.findall(s.lower()) if len(t) > 2 and t not in _STOP]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--index', required=True, help='text_string_index.tsv')
    ap.add_argument('--train-terms', required=True)
    ap.add_argument('--queries', required=True, help='file of query protein IDs (one per line)')
    ap.add_argument('--out', required=True)
    ap.add_argument('--topk', type=int, default=30, help='neighbours to vote (default 30)')
    ap.add_argument('--min-score', type=float, default=0.01)
    ap.add_argument('--shard', default='0/1',
                    help='i/N: process only query shard i of N (for parallel runs)')
    args = ap.parse_args()

    si, sn = (int(x) for x in args.shard.split('/'))
    queries = {l.strip() for l in open(args.queries) if l.strip()}
    log(f'queries: {len(queries):,} (shard {si}/{sn})')

    # corpus labels (pre-t0); exclude query proteins from the corpus
    labels = defaultdict(set)
    with open(args.train_terms) as fh:
        next(fh, None)
        for line in fh:
            p, t, _asp = line.rstrip('\n').split('\t')[:3]
            if p not in queries:
                labels[p].add(t)
    log(f'corpus-labelled proteins (train, minus queries): {len(labels):,}')

    # text from the index
    name_text, char_text = {}, {}
    with open(args.index) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 5:
                continue
            ac, _tax, _sid, nm, ch = parts[:5]
            name_text[ac] = nm
            char_text[ac] = ch

    # build corpus docs (name + characterization) for labelled training proteins with text
    corpus_ids, corpus_docs = [], []
    for p in labels:
        nm = name_text.get(p, '')
        ch = char_text.get(p, '')
        toks = tok(nm + ' ' + ch)
        if toks:
            corpus_ids.append(p)
            corpus_docs.append(toks)
    log(f'corpus docs with text: {len(corpus_docs):,}')
    bm25 = BM25Okapi(corpus_docs)
    corpus_labels = [labels[p] for p in corpus_ids]

    topk = args.topk
    n_out = 0
    shard_queries = [q for i, q in enumerate(sorted(queries)) if i % sn == si]
    log(f'this shard processes {len(shard_queries):,} queries')
    with open(args.out, 'w') as out:
        for qi, q in enumerate(shard_queries):
            qtoks = tok(name_text.get(q, ''))  # NAME ONLY -- temporally safe
            if not qtoks:
                continue
            scores = bm25.get_scores(qtoks)
            if not scores.any():
                continue
            nbr = np.argpartition(scores, -topk)[-topk:] if len(scores) > topk \
                else np.arange(len(scores))
            vote = defaultdict(float)
            for j in nbr:
                s = float(scores[j])
                if s <= 0:
                    continue
                for t in corpus_labels[j]:
                    vote[t] += s
            if not vote:
                continue
            mx = max(vote.values())
            for t, v in vote.items():
                sc = v / mx
                if sc >= args.min_score:
                    out.write(f'{q}\t{t}\t{sc:.4f}\n')
                    n_out += 1
            if (qi + 1) % 5000 == 0:
                log(f'  {qi+1:,} queries, {n_out:,} rows')
    log(f'done: {n_out:,} rows -> {args.out}')


if __name__ == '__main__':
    main()
