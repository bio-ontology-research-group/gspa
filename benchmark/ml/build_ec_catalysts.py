#!/usr/bin/env python3
"""Emit per-EC catalyst contributor list + merged panel embeddings.

Per-EC centroid = mean of catalyst contributors. For strict Track A we
need to exclude contributors sequence-similar to a per-case target, so
we persist *all* contributors (not just the final centroid).

Outputs:
  ec_panel.npy               (N_panel, 640) — L2-normalised panel embs
  ec_panel.index.tsv         global_row tag protein_id
  ec_catalysts.tsv           ec global_row  (one row per (ec, catalyst))
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

    # First pass: build merged panel embedding matrix
    all_rows = []  # (tag, protein_id)
    all_embs = []
    tag_to_offset = {}
    for tag in tags:
        npy = Path(args.plm_dir) / f'{tag}_esm2t30.npy'
        idx = Path(args.plm_dir) / f'{tag}_esm2t30.index.tsv'
        if not npy.exists() or not idx.exists():
            continue
        arr = np.load(npy)
        pids = []
        with open(idx) as f:
            f.readline()
            for line in f:
                pid, r = line.rstrip('\n').split('\t')
                pids.append((int(r), pid))
        pids.sort()
        tag_to_offset[tag] = (len(all_rows), [p for _, p in pids])
        for _, pid in pids:
            all_rows.append((tag, pid))
        # L2-normalise rows
        norm = np.linalg.norm(arr, axis=1, keepdims=True) + 1e-9
        all_embs.append((arr / norm).astype(np.float32))

    panel = np.concatenate(all_embs, axis=0)
    print(f'[info] merged panel: {panel.shape}', file=sys.stderr)

    # Build pid_to_global per tag
    tag_pid_to_global = {}
    for tag, (offset, pids) in tag_to_offset.items():
        tag_pid_to_global[tag] = {p: offset + i for i, p in enumerate(pids)}

    # Second pass: per-tag integrated → per-EC contributors
    ec_catalysts = collections.defaultdict(set)  # ec -> {global_row}
    for tag in tags:
        intg = Path(args.root) / 'integrated' / f'{tag}_integrated.tsv'
        if not intg.exists():
            continue
        go_to_prots = load_integrated_go_for_tag(str(intg), args.tau)
        tag_map = tag_pid_to_global.get(tag, {})
        for ec, go in ec2go.items():
            if ec.startswith('EC:'):
                ec = ec[3:]
            prots = go_to_prots.get(go)
            if not prots:
                continue
            for pid in prots:
                g = tag_map.get(pid)
                if g is not None:
                    ec_catalysts[ec].add(g)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    np.save(out_dir / 'ec_panel.npy', panel)
    with open(out_dir / 'ec_panel.index.tsv', 'w') as f:
        f.write('global_row\ttag\tprotein_id\n')
        for i, (t, p) in enumerate(all_rows):
            f.write(f'{i}\t{t}\t{p}\n')
    with open(out_dir / 'ec_catalysts.tsv', 'w') as f:
        f.write('ec\tglobal_row\n')
        n = 0
        for ec in sorted(ec_catalysts.keys()):
            for g in sorted(ec_catalysts[ec]):
                f.write(f'{ec}\t{g}\n')
                n += 1
    print(f'[info] wrote {len(ec_catalysts)} ECs with {n} catalyst rows',
          file=sys.stderr)


if __name__ == '__main__':
    main()
