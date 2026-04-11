#!/usr/bin/env python3
"""
Compute head-to-head metrics for GSPA vs PGAP on a held-out proteome.

Inputs per genome:
  - truth.tsv                   (experimental GOA ground truth, key = UniProt acc)
  - gspa_integrated.tsv         (gspa-cli integrate output, key = UniProt acc)
  - pgap_annotations.tsv        (extract_pgap_annotations.py output)

Metrics (all computed per aspect + overall):
  - F-max (threshold sweep over predicted scores)
  - Coverage (fraction of truth proteins with any prediction)
  - Unique GO terms predicted
  - IC-weighted recall (each correctly-predicted term contributes its IC)

Information Content is the -log2 of the term's frequency in a background
set (we use GOA's Swiss-Prot frequency).
"""

import argparse
import json
import math
import os
import sys
from collections import defaultdict


def load_truth(path):
    """Parse 'accession aspect go_term' into {acc: {aspect: set(go)}}"""
    out = defaultdict(lambda: defaultdict(set))
    with open(path) as fh:
        header = next(fh, None)
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            out[parts[0]][parts[1]].add(parts[2])
    return out


def load_gspa_integrated(path, acc_mapper=None):
    """Parse gspa integrate output and optionally remap keys.

    `acc_mapper` is a dict from the key in integrate.tsv (RefSeq or UniProt)
    to the desired key (UniProt for comparison with truth).
    Returns {acc: {aspect: {go: prob}}}.
    """
    out = defaultdict(lambda: defaultdict(dict))
    with open(path) as fh:
        header_line = next(fh, '')
        headers = header_line.rstrip('\n').split('\t')
        idx_pid = headers.index('protein_id')
        idx_type = headers.index('type')
        idx_func = headers.index('function_id')
        idx_aspect = headers.index('go_aspect')
        idx_prob = headers.index('posterior_prob')
        for line in fh:
            fields = line.rstrip('\n').split('\t')
            if len(fields) <= idx_prob or fields[idx_type] != 'GO':
                continue
            try:
                prob = float(fields[idx_prob])
            except ValueError:
                continue
            pid = fields[idx_pid]
            if acc_mapper is not None:
                pid = acc_mapper.get(pid, pid)
            aspect = fields[idx_aspect] or 'BP'
            out[pid][aspect][fields[idx_func]] = prob
    return out


def load_pgap_binary(path):
    """Parse PGAP TSV as binary predictions (score=1.0)."""
    out = defaultdict(lambda: defaultdict(dict))
    with open(path) as fh:
        header = next(fh, None)
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            out[parts[0]][parts[1]][parts[2]] = 1.0
    return out


def fmax(predictions, truth, aspect=None, thresholds=None):
    """Return (fmax, threshold)."""
    if thresholds is None:
        thresholds = [i * 0.05 for i in range(21)]
    best_f1 = 0.0
    best_t = 0.0
    for t in thresholds:
        tp = fp = fn = 0
        for acc in truth:
            for asp in truth[acc]:
                if aspect and asp != aspect:
                    continue
                true_set = set(truth[acc][asp])
                pred_set = {g for g, s in predictions.get(acc, {}).get(asp, {}).items() if s >= t}
                tp += len(true_set & pred_set)
                fp += len(pred_set - true_set)
                fn += len(true_set - pred_set)
        p = tp / (tp + fp) if tp + fp else 0.0
        r = tp / (tp + fn) if tp + fn else 0.0
        f1 = 2 * p * r / (p + r) if p + r else 0.0
        if f1 > best_f1:
            best_f1 = f1
            best_t = t
    return best_f1, best_t


def coverage(predictions, truth, threshold=0.5):
    """Fraction of truth proteins with >=1 prediction >= threshold."""
    denom = 0
    num = 0
    for acc in truth:
        denom += 1
        has = False
        for asp, terms in predictions.get(acc, {}).items():
            if any(s >= threshold for s in terms.values()):
                has = True
                break
        if has:
            num += 1
    return num / denom if denom else 0.0


def unique_terms(predictions, threshold=0.5):
    seen = set()
    for acc_preds in predictions.values():
        for asp_preds in acc_preds.values():
            for g, s in asp_preds.items():
                if s >= threshold:
                    seen.add(g)
    return len(seen)


def compute_ic(truth_set_path):
    """
    IC from truth_set_path frequency: ic(term) = -log2(count(term) / total_annotated_proteins).
    Returns dict go_term -> ic.
    """
    by_term = defaultdict(set)
    total_proteins = set()
    with open(truth_set_path) as fh:
        header = next(fh, None)
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            acc, aspect, go = parts[0], parts[1], parts[2]
            by_term[go].add(acc)
            total_proteins.add(acc)
    total = max(1, len(total_proteins))
    return {g: -math.log2(len(a) / total) for g, a in by_term.items() if a}


def ic_recall(predictions, truth, ic_map, threshold=0.5, aspect=None):
    """Sum IC of correctly-predicted (per-aspect) terms, normalized by total IC of truth."""
    total_ic = 0.0
    recovered_ic = 0.0
    for acc in truth:
        for asp in truth[acc]:
            if aspect and asp != aspect:
                continue
            true_set = set(truth[acc][asp])
            pred_set = {g for g, s in predictions.get(acc, {}).get(asp, {}).items() if s >= threshold}
            for g in true_set:
                ic = ic_map.get(g, 0.0)
                total_ic += ic
                if g in pred_set:
                    recovered_ic += ic
    if total_ic == 0:
        return 0.0
    return recovered_ic / total_ic


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--truth', required=True)
    p.add_argument('--gspa', required=True, help='GSPA integrated TSV')
    p.add_argument('--pgap', default=None, help='PGAP annotations TSV (optional)')
    p.add_argument('--gspa-key-map', default=None,
                   help='TSV RefSeq\\tUniProt to remap GSPA protein IDs (optional)')
    p.add_argument('--tag', required=True)
    args = p.parse_args()

    truth = load_truth(args.truth)
    truth_total_annotations = sum(
        sum(len(terms) for terms in asps.values())
        for asps in truth.values()
    )
    truth_proteins = len(truth)

    # Optional protein ID remap for GSPA (RefSeq -> UniProt).
    acc_mapper = None
    if args.gspa_key_map:
        acc_mapper = {}
        with open(args.gspa_key_map) as fh:
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 2:
                    acc_mapper[parts[0]] = parts[1]

    gspa = load_gspa_integrated(args.gspa, acc_mapper=acc_mapper)
    pgap = load_pgap_binary(args.pgap) if args.pgap else None

    ic_map = compute_ic(args.truth)

    def summarize(name, preds):
        row = {'method': name}
        fmax_o, t_o = fmax(preds, truth)
        row['fmax_overall'] = fmax_o
        row['fmax_t'] = t_o
        for asp in ('MF', 'BP', 'CC'):
            f, _ = fmax(preds, truth, aspect=asp)
            row[f'fmax_{asp}'] = f
        row['coverage'] = coverage(preds, truth, threshold=0.5)
        row['unique_terms'] = unique_terms(preds, threshold=0.5)
        row['ic_recall'] = ic_recall(preds, truth, ic_map, threshold=0.5)
        return row

    results = [summarize('GSPA', gspa)]
    if pgap is not None:
        results.append(summarize('PGAP', pgap))

    print(json.dumps({
        'tag': args.tag,
        'truth_proteins': truth_proteins,
        'truth_annotations': truth_total_annotations,
        'results': results,
    }, indent=2))


if __name__ == '__main__':
    main()
