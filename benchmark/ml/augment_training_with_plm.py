#!/usr/bin/env python3
"""Add PLM features to per-genome training.tsv, keyed by EC.

Adds: plm_cos_centroid_EC, plm_has_emb
"""
import argparse
import sys
from pathlib import Path

import numpy as np


def load_centroids(centroids_dir):
    arr = np.load(Path(centroids_dir) / 'ec_centroids.npy')
    ec_to_row = {}
    with open(Path(centroids_dir) / 'ec_centroids.index.tsv') as f:
        f.readline()
        for line in f:
            ec, r, n = line.rstrip('\n').split('\t')
            ec_to_row[ec] = int(r)
    return arr, ec_to_row


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


def cos_sim(a, b):
    na = np.linalg.norm(a) + 1e-9
    nb = np.linalg.norm(b) + 1e-9
    return float(np.dot(a, b) / (na * nb))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--in-dir', required=True)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--plm-dir', required=True)
    ap.add_argument('--centroids-dir', required=True)
    args = ap.parse_args()

    centroids, ec_to_row = load_centroids(args.centroids_dir)
    print(f'[info] centroids: {centroids.shape}, {len(ec_to_row)} ECs',
          file=sys.stderr)

    in_root = Path(args.in_dir)
    out_root = Path(args.out_dir)
    out_root.mkdir(parents=True, exist_ok=True)

    n_ok = 0
    for tag_dir in sorted(in_root.iterdir()):
        if not tag_dir.is_dir():
            continue
        tag = tag_dir.name
        fin = tag_dir / 'training.tsv'
        if not fin.exists() or fin.stat().st_size == 0:
            continue
        arr, pid_to_row = load_embeddings(args.plm_dir, tag)
        if arr is None:
            continue
        out_tag = out_root / tag
        out_tag.mkdir(parents=True, exist_ok=True)
        fout = out_tag / 'training.tsv'
        n_rows = 0
        hit_ec = 0
        with open(fin) as f, open(fout, 'w') as g:
            header = f.readline().rstrip('\n').split('\t')
            new_header = header + ['plm_cos_centroid_EC', 'plm_has_emb']
            g.write('\t'.join(new_header) + '\n')
            p_i = header.index('protein_id')
            ec_i = header.index('ec')
            for line in f:
                parts = line.rstrip('\n').split('\t')
                if len(parts) < len(header):
                    continue
                pid = parts[p_i]
                ec = parts[ec_i].strip()
                if ec.startswith('EC:'):
                    ec = ec[3:]
                val = 0.0
                has = 0.0
                row = pid_to_row.get(pid)
                if row is not None:
                    has = 1.0
                erow = ec_to_row.get(ec)
                if row is not None and erow is not None:
                    val = cos_sim(arr[row].astype(np.float32),
                                  centroids[erow])
                    hit_ec += 1
                g.write('\t'.join(parts + [f'{val:.6f}', f'{has:.1f}']) + '\n')
                n_rows += 1
        print(f'[info] {tag}: {n_rows} rows, {hit_ec} EC-matched',
              file=sys.stderr)
        n_ok += 1
    print(f'[info] done, {n_ok} tags augmented', file=sys.stderr)


if __name__ == '__main__':
    main()
