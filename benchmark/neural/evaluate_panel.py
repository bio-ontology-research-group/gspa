#!/usr/bin/env python3
"""Evaluate per-genome predictor output against SwissProt truth.

Produces a TSV with one row per (genome, predictor, truth-evidence-set)
with columns: fmax_overall, fmax_cafa, fmax_MF/BP/CC, smin_overall, and
(if EC) ec_level1..4 F-max. Optional bootstrap CIs.

Truth file format (tab-separated, header): protein_id<TAB>aspect<TAB>function_id
Predictor TSV format (header): protein_id<TAB>term<TAB>score<TAB>annotation_type

The GO track loads only rows with annotation_type=GO. The EC track loads
only rows with annotation_type=EC and uses a separate EC truth file.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

# Import F-max / Smin / EC helpers from the unified benchmark module.
BENCH_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(BENCH_DIR))
import benchmark_pgap_v2 as B  # noqa: E402


def load_predictor_tsv(path: Path, ann_type: str,
                       go_aspect_map: dict[str, str] | None = None,
                       protein_filter: set | None = None):
    """Load a predictor TSV into preds[acc][aspect]->{term:score}.
    GO: aspect is looked up per term from go_aspect_map (required if GO).
    EC: aspect is '' for all rows.
    If ``protein_filter`` is given, keep only rows whose protein_id is in the set.
    """
    out = defaultdict(lambda: defaultdict(dict))
    with path.open() as fh:
        header = next(fh).rstrip('\n').split('\t')
        idx = {c: i for i, c in enumerate(header)}
        needed = ('protein_id', 'term', 'score', 'annotation_type')
        for c in needed:
            if c not in idx:
                raise SystemExit(f"{path} missing column {c}")
        for line in fh:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 4:
                continue
            t = fields[idx['annotation_type']]
            if t != ann_type:
                continue
            pid = fields[idx['protein_id']]
            if protein_filter is not None and pid not in protein_filter:
                continue
            term = fields[idx['term']]
            try:
                score = float(fields[idx['score']])
            except ValueError:
                continue
            if ann_type == 'GO':
                asp = go_aspect_map.get(term, 'BP') if go_aspect_map else 'BP'
            else:
                asp = ''
            prev = out[pid][asp].get(term, 0.0)
            if score > prev:
                out[pid][asp][term] = score
    return out


def load_go_aspect_map(path: Path) -> dict[str, str]:
    """Load GO term → aspect map (MF/BP/CC). Two formats supported:
      - tsv with header: term<TAB>aspect
      - OBO: parse [Term] blocks for id and namespace.
    """
    m: dict[str, str] = {}
    txt = path.read_text()
    if '[Term]' in txt:
        cur_id = None
        for line in txt.split('\n'):
            line = line.strip()
            if line == '[Term]':
                cur_id = None
            elif line.startswith('id: '):
                cur_id = line[4:].strip()
            elif line.startswith('namespace: ') and cur_id:
                ns = line[11:].strip()
                asp = {'molecular_function': 'MF',
                       'biological_process': 'BP',
                       'cellular_component': 'CC'}.get(ns)
                if asp:
                    m[cur_id] = asp
    else:
        with path.open() as fh:
            header = fh.readline()
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 2:
                    m[parts[0]] = parts[1]
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--predictor-tsv', required=True,
                    help='path to <tag>.<predictor>.tsv')
    ap.add_argument('--truth', required=True,
                    help='truth TSV: protein_id\\taspect\\tfunction_id')
    ap.add_argument('--annotation-type', choices=['GO', 'EC'], default='GO')
    ap.add_argument('--go-aspect-map', default=None,
                    help='for GO: TSV(term,aspect) or go.obo for aspect lookup')
    ap.add_argument('--tag', required=True)
    ap.add_argument('--predictor', required=True)
    ap.add_argument('--truth-name', default='all')
    ap.add_argument('--n-bootstrap', type=int, default=0)
    args = ap.parse_args()

    # 1) truth
    truth = B.load_truth(args.truth)
    ic_map = B.compute_ic(args.truth)
    n_truth_prot = len(truth)
    n_truth_ann = sum(sum(len(t) for t in asps.values()) for asps in truth.values())

    # 2) predictor — filter to proteins with truth to bound memory
    aspect_map = None
    if args.annotation_type == 'GO' and args.go_aspect_map:
        aspect_map = load_go_aspect_map(Path(args.go_aspect_map))
    preds = load_predictor_tsv(Path(args.predictor_tsv), args.annotation_type,
                               aspect_map, protein_filter=set(truth.keys()))

    n_pred_prot = len(preds)
    n_pred_ann = sum(sum(len(t) for t in asps.values()) for asps in preds.values())

    row = B.summarize(args.predictor, preds, truth, ic_map, args.n_bootstrap,
                      ann_type=args.annotation_type)
    row['tag'] = args.tag
    row['truth'] = args.truth_name
    row['annotation_type'] = args.annotation_type
    row['n_truth_proteins'] = n_truth_prot
    row['n_truth_annotations'] = n_truth_ann
    row['n_pred_proteins'] = n_pred_prot
    row['n_pred_annotations'] = n_pred_ann
    print(json.dumps(row, indent=2))


if __name__ == '__main__':
    main()
