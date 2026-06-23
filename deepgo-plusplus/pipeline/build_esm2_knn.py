#!/usr/bin/env python3
"""goPredSim-style ESM2 embedding-kNN GO component (the orphan-tier workhorse).

Transfers GO labels from the nearest pre-t0 train proteins by cosine similarity in
ESM2-35M embedding space. CPU-only. Two modes:

  * LOO  (--query == --ref, --loo): leave-one-out over the train store -> the
    leak-free TRAINING feature for Integrator-B (synthetic orphans). Excludes self;
    also drops near-identical neighbours (sim >= --max-sim) so the feature matches
    the orphan inference regime (a novel protein has no train twin).
  * transfer (--query test store, --ref train store): the INFERENCE feature — test
    orphans vs the train reference.

Out: protein<TAB>GO<TAB>score   (sim-weighted vote, max-normalized per protein).
"""
import argparse, sys, time
from collections import defaultdict
import numpy as np


def log(m): print(f'[{time.strftime("%H:%M:%S")}] {m}', file=sys.stderr, flush=True)


def load_store(path):
    d = np.load(path, allow_pickle=True)
    emb = d['emb'].astype(np.float32)
    emb /= (np.linalg.norm(emb, axis=1, keepdims=True) + 1e-8)
    return list(d['ids']), emb


def load_labels(path):
    y = defaultdict(set)
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) >= 2:
                y[p[0]].add(p[1])
    return y


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--query', required=True, help='query .npz store')
    ap.add_argument('--ref', required=True, help='reference (train) .npz store')
    ap.add_argument('--ref-labels', required=True, help='train_terms (EntryID,term,aspect)')
    ap.add_argument('--out', required=True)
    ap.add_argument('--topk', type=int, default=10)
    ap.add_argument('--loo', action='store_true', help='leave-one-out (query==ref)')
    ap.add_argument('--max-sim', type=float, default=0.999,
                    help='drop neighbours with sim >= this (near-identical) in LOO mode')
    ap.add_argument('--min-score', type=float, default=0.01)
    ap.add_argument('--chunk', type=int, default=1000)
    args = ap.parse_args()

    qids, Q = load_store(args.query)
    rids, R = load_store(args.ref)
    ylab = load_labels(args.ref_labels)
    log(f'query {Q.shape}  ref {R.shape}  labelled-refs {len(ylab):,}')
    rid_index = {r: i for i, r in enumerate(rids)} if args.loo else None

    fout = open(args.out, 'w')
    RT = R.T
    for s in range(0, len(qids), args.chunk):
        chunk = qids[s:s + args.chunk]
        sims = Q[s:s + args.chunk] @ RT               # (c, nref)
        for r, qid in enumerate(chunk):
            row = sims[r]
            if args.loo:
                j = rid_index.get(qid)
                if j is not None:
                    row[j] = -2.0                     # exclude self
                row[row >= args.max_sim] = -2.0       # exclude near-identical twins
            nn = np.argpartition(-row, args.topk)[:args.topk]
            vote, wsum = defaultdict(float), 0.0
            for k in nn:
                w = float(row[k])
                if w <= 0:
                    continue
                wsum += w
                for t in ylab.get(rids[k], ()):
                    vote[t] += w
            if wsum <= 0:
                continue
            vmax = max(vote.values())
            for t, v in vote.items():
                sc = v / vmax
                if sc >= args.min_score:
                    fout.write(f'{qid}\t{t}\t{sc:.4f}\n')
        log(f'  {min(s+args.chunk,len(qids)):,}/{len(qids):,}')
    fout.close()
    log(f'wrote {args.out}')


if __name__ == '__main__':
    main()
