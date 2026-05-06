#!/usr/bin/env python3
"""Build per-EC PLM centroids from panel catalysts.

Keyed by EC number (not reaction_id) so it aligns with training-data
candidates that have an `ec` column in MetaCyc-namespace.

For each EC number E, the centroid is the L2-normalised mean embedding
of all panel proteins whose integrated posterior on the GO term for E
exceeds tau.

Output:
  {out_dir}/ec_centroids.npy       (float32, N_ec x 640)
  {out_dir}/ec_centroids.index.tsv (ec -> row, n_catalysts)
"""
import argparse
import collections
import re
import sys
from pathlib import Path

import numpy as np


def load_ec2go(path):
    ec2go = {}
    with open(path) as f:
        for line in f:
            if line.startswith('!') or not line.strip():
                continue
            m = re.match(r'^EC:(\S+)\s*>\s*GO:[^;]+;\s*(GO:\d+)', line)
            if m:
                ec2go[m.group(1)] = m.group(2)
    return ec2go


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


def load_integrated_go_for_tag(path, tau=0.3):
    out = collections.defaultdict(set)
    with open(path) as f:
        f.readline()
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 5 or parts[1] != 'GO':
                continue
            try:
                p = float(parts[4])
            except ValueError:
                continue
            if p <= tau:
                continue
            out[parts[2]].add(parts[0])
    return dict(out)


def load_embeddings(plm_dir, tag):
    npy = Path(plm_dir) / f'{tag}_esm2t30.npy'
    idx = Path(plm_dir) / f'{tag}_esm2t30.index.tsv'
    if not npy.exists() or not idx.exists():
        return None, {}
    arr = np.load(npy, mmap_mode='r')
    pid_to_row = {}
    with open(idx) as f:
        f.readline()
        for line in f:
            pid, r = line.rstrip('\n').split('\t')
            pid_to_row[pid] = int(r)
    return arr, pid_to_row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--manifest', required=True)
    ap.add_argument('--root', required=True)
    ap.add_argument('--plm-dir', required=True)
    ap.add_argument('--ec2go', required=True)
    ap.add_argument('--exclude-tag', default='mg1655')
    ap.add_argument('--tau', type=float, default=0.3)
    ap.add_argument('--out-dir', required=True)
    args = ap.parse_args()

    tags = load_panel_tags(
        args.manifest,
        exclude=args.exclude_tag.split(',') if args.exclude_tag else None)
    print(f'[info] panel size: {len(tags)}', file=sys.stderr)

    ec2go = load_ec2go(args.ec2go)
    print(f'[info] ec2go: {len(ec2go)} entries', file=sys.stderr)

    # EC -> accumulated sum, count
    ec_sum = {}
    ec_cnt = collections.Counter()

    for tag in tags:
        intg = Path(args.root) / 'integrated' / f'{tag}_integrated.tsv'
        if not intg.exists():
            continue
        go_to_prots = load_integrated_go_for_tag(str(intg), args.tau)
        if not go_to_prots:
            continue
        arr, pid_to_row = load_embeddings(args.plm_dir, tag)
        if arr is None:
            continue
        touched = 0
        for ec, go in ec2go.items():
            if ec.startswith('EC:'):
                ec = ec[3:]
            prots = go_to_prots.get(go)
            if not prots:
                continue
            seen_rows = set()
            for pid in prots:
                row = pid_to_row.get(pid)
                if row is None:
                    continue
                seen_rows.add(row)
            for r in seen_rows:
                v = arr[r].astype(np.float32)
                n = np.linalg.norm(v) + 1e-9
                v = v / n
                if ec not in ec_sum:
                    ec_sum[ec] = np.zeros(640, dtype=np.float32)
                ec_sum[ec] += v
                ec_cnt[ec] += 1
            if seen_rows:
                touched += 1
        print(f'[info] {tag}: touched {touched} ECs', file=sys.stderr)

    if not ec_sum:
        print('[error] no EC centroids built', file=sys.stderr)
        sys.exit(1)

    ecs = sorted(ec_sum.keys())
    centroids = np.zeros((len(ecs), 640), dtype=np.float32)
    for i, e in enumerate(ecs):
        c = ec_sum[e] / ec_cnt[e]
        n = np.linalg.norm(c) + 1e-9
        centroids[i] = c / n

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    np.save(out_dir / 'ec_centroids.npy', centroids)
    with open(out_dir / 'ec_centroids.index.tsv', 'w') as f:
        f.write('ec\trow\tn_catalysts\n')
        for i, e in enumerate(ecs):
            f.write(f'{e}\t{i}\t{ec_cnt[e]}\n')
    print(f'[info] wrote {len(ecs)} EC-centroids to {out_dir}',
          file=sys.stderr)


if __name__ == '__main__':
    main()
