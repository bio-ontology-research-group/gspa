#!/usr/bin/env python3
"""Build the STRING Net-KNN components from the FULL STRING aggregate, keyed on
string_id (NOT UniProt OX taxon) so the per-species taxon mismatch is gone.

Background: STRING names everything by its OWN taxon (the prefix of the string_id,
e.g. `4932.YAL001C`), which differs from the UniProt OX taxon for many model
organisms (yeast OX 559292 -> STRING 4932; E. coli K-12 OX 83333 -> STRING 511145).
The old per-species builders grouped/downloaded by OX taxon, so the PPI networks of
the most-studied organisms were silently absent. This builder subsets the full
STRING aggregate (`protein.links.v12.0.txt.gz`) by the exact anchor string_ids, so
every edge is matched by STRING ID and the mismatch can't occur.

Two-stage:
  1. (awk, outside) subset the aggregate to edges touching an anchor string_id
     (anchor = query string-node UNION homolog string-node), score >= min-conf,
     -> edges.tsv  (p1 p2 score).
  2. (this script) build `net` (direct), `net_bridge` (DIAMOND homology hop) and
     `net_union` (max merge) for the query proteins.

Scoring matches the original build_net_component.py / build_net_bridge.py exactly
(confidence-weighted neighbour vote, per-protein max-normalised) so the eval
isolates the DATA effect (full STRING + mismatch fix) from any scoring change.
A `--ic-weight` flag optionally switches to the 5th-place CAFA6 scheme
(IC-weighted vote, normalised by sum-of-weights) for a separate ablation.

Inputs:
  --index        net_index.tsv  (accession\\ttaxon\\tstring_id\\t...)
  --train-terms  pre-t0 labels (EntryID\\tterm\\taspect)
  --queries      query FASTA (or one id per line)
  --edges        edges.tsv from the awk subset (p1\\tp2\\tscore, whitespace ok)
  --diamond      m8 queries vs train DB: qseqid sseqid bitscore [pident]
  --out-prefix   writes <prefix>_net.tsv, <prefix>_net_bridge.tsv, <prefix>_net_union.tsv
"""
from __future__ import annotations

import argparse
import math
import sys
import time
from collections import defaultdict


def log(m):
    print(f'[{time.strftime("%H:%M:%S")}] {m}', file=sys.stderr, flush=True)


def load_queries(path):
    ids = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            ids.append(line[1:].split()[0] if line.startswith('>') else line.split()[0])
    return list(dict.fromkeys(ids))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--index', required=True)
    ap.add_argument('--train-terms', required=True)
    ap.add_argument('--queries', required=True)
    ap.add_argument('--edges', required=True)
    ap.add_argument('--diamond', required=True)
    ap.add_argument('--out-prefix', required=True)
    ap.add_argument('--min-conf', type=int, default=400)
    ap.add_argument('--topk', type=int, default=50)
    ap.add_argument('--topk-homologs', type=int, default=5)
    ap.add_argument('--min-score', type=float, default=0.01)
    ap.add_argument('--ic-weight', action='store_true',
                    help='5th-place scheme: multiply votes by term IC, normalise by sum-of-weights')
    args = ap.parse_args()

    # index: accession -> string_id
    ac_sid = {}
    with open(args.index) as fh:
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) >= 3 and p[2]:
                ac_sid[p[0]] = p[2]
    log(f'accessions with string_id: {len(ac_sid):,}')

    # string_id -> set(GO) from pre-t0 train labels; also term frequency for IC
    sid_labels = defaultdict(set)
    term_freq = defaultdict(int)
    with open(args.train_terms) as fh:
        next(fh, None)
        for line in fh:
            row = line.rstrip('\n').split('\t')
            sid = ac_sid.get(row[0])
            if sid:
                sid_labels[sid].add(row[1])
    n_train_prot = 0
    if args.ic_weight:
        # IC from training-set occurrence probability (per the 5th-place formula)
        prot_terms = defaultdict(set)
        with open(args.train_terms) as fh:
            next(fh, None)
            for line in fh:
                row = line.rstrip('\n').split('\t')
                prot_terms[row[0]].add(row[1])
        n_train_prot = len(prot_terms)
        for terms in prot_terms.values():
            for t in terms:
                term_freq[t] += 1
    log(f'string nodes with pre-t0 labels: {len(sid_labels):,}')

    # queries and their string nodes
    queries = load_queries(args.queries)
    q_sid = {q: ac_sid[q] for q in queries if q in ac_sid}
    log(f'queries: {len(queries):,}; with string node: {len(q_sid):,}')

    # adjacency for anchor nodes from the subset edges
    adj = defaultdict(list)   # anchor_sid -> [(nbr_sid, conf)]
    anchors = set(q_sid.values())
    # homolog anchors from diamond
    qset = set(queries)
    homolog_sid = {}          # homolog accession -> string_id
    with open(args.diamond) as fh:
        for line in fh:
            c = line.rstrip('\n').split('\t')
            if len(c) < 3 or c[0] == c[1] or c[0] not in qset:
                continue
            h = c[1]
            if h in ac_sid:
                homolog_sid[h] = ac_sid[h]
    anchors |= set(homolog_sid.values())
    log(f'anchor string nodes: {len(anchors):,}')

    n_edge = 0
    with open(args.edges) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 3:
                continue
            a, b, sc = parts[0], parts[1], parts[2]
            try:
                c = int(sc)
            except ValueError:
                continue
            if c < args.min_conf:
                continue
            if a in anchors:
                adj[a].append((b, c))
            if b in anchors:
                adj[b].append((a, c))
            n_edge += 1
    log(f'edges kept: {n_edge:,}; anchors with neighbours: {len(adj):,}')

    max_ic = 1.0
    if args.ic_weight and term_freq:
        max_ic = max(-math.log2(f / n_train_prot) for f in term_freq.values())

    def vote(nbrs):
        """neighbour list -> {term: score} per-protein normalised."""
        nbrs = sorted(nbrs, key=lambda x: -x[1])[:args.topk]
        v = defaultdict(float)
        wsum = 0.0
        for nb_sid, conf in nbrs:
            w = conf / 1000.0
            wsum += w
            for t in sid_labels.get(nb_sid, ()):
                if args.ic_weight:
                    ic = (-math.log2(term_freq[t] / n_train_prot) / max_ic) if term_freq.get(t) else 0.0
                    v[t] += w * ic
                else:
                    v[t] += w
        if not v:
            return {}
        if args.ic_weight and wsum > 0:
            for t in v:
                v[t] /= wsum
            mx = max(v.values()) or 1.0
            return {t: s / mx for t, s in v.items()}   # keep [0,1] for downstream merge
        mx = max(v.values())
        return {t: s / mx for t, s in v.items()}

    # ---- direct net ----
    net = {}
    for q, sid in q_sid.items():
        if sid in adj:
            vt = vote(adj[sid])
            if vt:
                net[q] = vt
    # ---- bridge ----
    # precompute net(h) for homolog nodes
    hsid_vote = {}
    for h, sid in homolog_sid.items():
        if sid in adj and sid not in hsid_vote:
            hsid_vote[sid] = vote(adj[sid])
    # diamond homologs per query
    q_hom = defaultdict(list)
    with open(args.diamond) as fh:
        for line in fh:
            c = line.rstrip('\n').split('\t')
            if len(c) < 3 or c[0] == c[1] or c[0] not in qset:
                continue
            h = c[1]
            if h not in homolog_sid:
                continue
            try:
                q_hom[c[0]].append((h, float(c[2])))
            except ValueError:
                pass
    bridge = {}
    for q, hs in q_hom.items():
        hs = sorted(hs, key=lambda x: -x[1])[:args.topk_homologs]
        if not hs:
            continue
        mxb = hs[0][1] or 1.0
        v = defaultdict(float)
        qsid = q_sid.get(q)
        for h, b in hs:
            sid = homolog_sid[h]
            if sid == qsid:
                continue
            w = b / mxb
            for t, s in hsid_vote.get(sid, {}).items():
                v[t] += w * s
        if v:
            mx = max(v.values())
            bridge[q] = {t: s / mx for t, s in v.items()}

    def write(d, path):
        n = 0
        with open(path, 'w') as out:
            for q, vt in d.items():
                for t, s in vt.items():
                    if s >= args.min_score:
                        out.write(f'{q}\t{t}\t{s:.4f}\n')
                        n += 1
        log(f'  {n:,} rows -> {path}')

    write(net, f'{args.out_prefix}_net.tsv')
    write(bridge, f'{args.out_prefix}_net_bridge.tsv')
    # union = max merge
    union = {}
    for q in set(net) | set(bridge):
        m = dict(net.get(q, {}))
        for t, s in bridge.get(q, {}).items():
            if s > m.get(t, 0):
                m[t] = s
        union[q] = m
    write(union, f'{args.out_prefix}_net_union.tsv')


if __name__ == '__main__':
    main()
