#!/usr/bin/env python3
"""Patched benchmark script: per-genome micro-averaged F-max + CAFA-style protein-centric F-max."""
import argparse
import json
import math
import random
from collections import defaultdict


def load_truth(path):
    out = defaultdict(lambda: defaultdict(set))
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            out[parts[0]][parts[1]].add(parts[2])
    return out


def load_gspa_integrated(path, acc_mapper=None):
    out = defaultdict(lambda: defaultdict(dict))
    with open(path) as fh:
        headers = next(fh, '').rstrip('\n').split('\t')
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
    out = defaultdict(lambda: defaultdict(dict))
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            out[parts[0]][parts[1]][parts[2]] = 1.0
    return out


# ============================================================
# Per-genome MICRO-averaged F-max (existing behaviour)
# TP/FP/FN summed across all (protein, GO-term) pairs in the genome,
# one F1 from those global sums, max over thresholds.
# ============================================================

def _count_per_acc(predictions, truth, threshold, aspect):
    per_acc = []
    for acc in truth:
        acc_tp = acc_fp = acc_fn = 0
        for asp in truth[acc]:
            if aspect and asp != aspect:
                continue
            true_set = set(truth[acc][asp])
            pred_set = {g for g, s in predictions.get(acc, {}).get(asp, {}).items() if s >= threshold}
            acc_tp += len(true_set & pred_set)
            acc_fp += len(pred_set - true_set)
            acc_fn += len(true_set - pred_set)
        per_acc.append((acc, acc_tp, acc_fp, acc_fn))
    return per_acc


def _f1_from_totals(tp, fp, fn):
    p = tp / (tp + fp) if tp + fp else 0.0
    r = tp / (tp + fn) if tp + fn else 0.0
    return 2 * p * r / (p + r) if p + r else 0.0


def fmax_with_ci(predictions, truth, aspect=None, thresholds=None,
                 n_bootstrap=200, seed=42):
    if thresholds is None:
        thresholds = [i * 0.05 for i in range(1, 21)]
    counts_by_t = {t: _count_per_acc(predictions, truth, t, aspect) for t in thresholds}
    totals = {t: (sum(c[1] for c in counts_by_t[t]),
                  sum(c[2] for c in counts_by_t[t]),
                  sum(c[3] for c in counts_by_t[t])) for t in thresholds}
    best_f = 0.0
    best_t = 0.0
    for t in thresholds:
        f = _f1_from_totals(*totals[t])
        if f > best_f:
            best_f = f
            best_t = t

    accs_at_t = counts_by_t[best_t] if best_t else None
    if not accs_at_t or n_bootstrap <= 0:
        return {'fmax': best_f, 'threshold': best_t, 'ci_low': best_f, 'ci_high': best_f}
    rng = random.Random(seed)
    n = len(accs_at_t)
    if n == 0:
        return {'fmax': best_f, 'threshold': best_t, 'ci_low': 0.0, 'ci_high': 0.0}
    boots = []
    for _ in range(n_bootstrap):
        sample_idx = [rng.randrange(n) for _ in range(n)]
        best_b = 0.0
        for t in thresholds:
            ca = counts_by_t[t]
            tp = fp = fn = 0
            for i in sample_idx:
                _, a_tp, a_fp, a_fn = ca[i]
                tp += a_tp
                fp += a_fp
                fn += a_fn
            f = _f1_from_totals(tp, fp, fn)
            if f > best_b:
                best_b = f
        boots.append(best_b)
    boots.sort()
    lo = boots[int(0.025 * n_bootstrap)]
    hi = boots[int(0.975 * n_bootstrap) - 1]
    return {'fmax': best_f, 'threshold': best_t, 'ci_low': lo, 'ci_high': hi}


# ============================================================
# CAFA-style protein-centric F-max
# For each threshold:
#   precision_p = |pred_p ∩ truth_p| / |pred_p|   (only for proteins with >=1 prediction)
#   recall_p    = |pred_p ∩ truth_p| / |truth_p|  (for all proteins with >=1 truth annotation)
#   avg_precision = mean(precision_p) over m(t) proteins with predictions
#   avg_recall    = mean(recall_p)   over n_e proteins with truth (proteins lacking
#                                    predictions contribute recall = 0)
#   F1 = 2*avg_p*avg_r / (avg_p+avg_r)
# F-max = max F1.
# This is the standard CAFA III/IV protocol (without IC-weighting; that's S-max).
# ============================================================

def _cafa_per_acc(predictions, truth, threshold, aspect):
    """Return per-protein (precision_or_None, recall) at this threshold."""
    out = []
    for acc in truth:
        true_set = set()
        for asp in truth[acc]:
            if aspect and asp != aspect:
                continue
            true_set |= truth[acc][asp]
        if not true_set:
            continue
        pred_set = set()
        for asp, terms in predictions.get(acc, {}).items():
            if aspect and asp != aspect:
                continue
            for g, s in terms.items():
                if s >= threshold:
                    pred_set.add(g)
        inter = len(true_set & pred_set)
        if pred_set:
            precision = inter / len(pred_set)
        else:
            precision = None  # this protein doesn't contribute to avg_precision
        recall = inter / len(true_set)
        out.append((acc, precision, recall))
    return out


def fmax_cafa_with_ci(predictions, truth, aspect=None, thresholds=None,
                      n_bootstrap=200, seed=42):
    if thresholds is None:
        thresholds = [i * 0.05 for i in range(1, 21)]
    per_t = {t: _cafa_per_acc(predictions, truth, t, aspect) for t in thresholds}

    def f1_from(rows):
        precisions = [r[1] for r in rows if r[1] is not None]
        recalls = [r[2] for r in rows]
        if not recalls:
            return 0.0
        avg_p = sum(precisions) / len(precisions) if precisions else 0.0
        avg_r = sum(recalls) / len(recalls)
        if avg_p + avg_r == 0:
            return 0.0
        return 2 * avg_p * avg_r / (avg_p + avg_r)

    best_f = 0.0
    best_t = 0.0
    for t in thresholds:
        f = f1_from(per_t[t])
        if f > best_f:
            best_f = f
            best_t = t

    rows_at_best = per_t[best_t] if best_t else (per_t[thresholds[0]] if thresholds else [])
    if not rows_at_best or n_bootstrap <= 0:
        return {'fmax': best_f, 'threshold': best_t, 'ci_low': best_f, 'ci_high': best_f}
    rng = random.Random(seed)
    n = len(rows_at_best)
    boots = []
    for _ in range(n_bootstrap):
        sample_idx = [rng.randrange(n) for _ in range(n)]
        best_b = 0.0
        for t in thresholds:
            rows = per_t[t]
            sample_rows = [rows[i] for i in sample_idx]
            f = f1_from(sample_rows)
            if f > best_b:
                best_b = f
        boots.append(best_b)
    boots.sort()
    lo = boots[int(0.025 * n_bootstrap)]
    hi = boots[int(0.975 * n_bootstrap) - 1]
    return {'fmax': best_f, 'threshold': best_t, 'ci_low': lo, 'ci_high': hi}


# ============================================================
# Other metrics
# ============================================================

def coverage(predictions, truth, threshold=0.5):
    denom = len(truth)
    num = 0
    for acc in truth:
        for asp, terms in predictions.get(acc, {}).items():
            if any(s >= threshold for s in terms.values()):
                num += 1
                break
    return num / denom if denom else 0.0


def unique_terms(predictions, threshold=0.5):
    seen = set()
    for acc_preds in predictions.values():
        for asp_preds in acc_preds.values():
            for g, s in asp_preds.items():
                if s >= threshold:
                    seen.add(g)
    return len(seen)


def compute_ic(truth_path):
    by_term = defaultdict(set)
    total = set()
    with open(truth_path) as fh:
        next(fh, None)
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            by_term[parts[2]].add(parts[0])
            total.add(parts[0])
    n = max(1, len(total))
    return {g: -math.log2(len(a) / n) for g, a in by_term.items() if a}


def ic_recall(predictions, truth, ic_map, threshold=0.5, aspect=None):
    total = 0.0
    recovered = 0.0
    for acc in truth:
        for asp in truth[acc]:
            if aspect and asp != aspect:
                continue
            true_set = set(truth[acc][asp])
            pred_set = {g for g, s in predictions.get(acc, {}).get(asp, {}).items() if s >= threshold}
            for g in true_set:
                ic = ic_map.get(g, 0.0)
                total += ic
                if g in pred_set:
                    recovered += ic
    return recovered / total if total else 0.0


def summarize(name, preds, truth, ic_map, n_bootstrap):
    row = {'method': name}
    # Per-genome micro-averaged F-max (existing)
    f = fmax_with_ci(preds, truth, n_bootstrap=n_bootstrap)
    row['fmax_overall'] = f['fmax']
    row['fmax_ci'] = [f['ci_low'], f['ci_high']]
    row['fmax_t'] = f['threshold']
    for asp in ('MF', 'BP', 'CC'):
        fa = fmax_with_ci(preds, truth, aspect=asp, n_bootstrap=n_bootstrap)
        row[f'fmax_{asp}'] = fa['fmax']
        row[f'fmax_{asp}_ci'] = [fa['ci_low'], fa['ci_high']]
    # CAFA-style protein-centric F-max (new)
    fc = fmax_cafa_with_ci(preds, truth, n_bootstrap=n_bootstrap)
    row['fmax_cafa_overall'] = fc['fmax']
    row['fmax_cafa_ci'] = [fc['ci_low'], fc['ci_high']]
    row['fmax_cafa_t'] = fc['threshold']
    for asp in ('MF', 'BP', 'CC'):
        fca = fmax_cafa_with_ci(preds, truth, aspect=asp, n_bootstrap=n_bootstrap)
        row[f'fmax_cafa_{asp}'] = fca['fmax']
        row[f'fmax_cafa_{asp}_ci'] = [fca['ci_low'], fca['ci_high']]
    row['coverage'] = coverage(preds, truth)
    row['unique_terms'] = unique_terms(preds)
    row['ic_recall'] = ic_recall(preds, truth, ic_map)
    return row


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--truth', action='append', required=True,
                   help='Repeatable: name:path (e.g., exp:truth_exp.tsv)')
    p.add_argument('--gspa', required=True)
    p.add_argument('--gspa-priors', default=None,
                   help='Optional second GSPA file (with-priors run)')
    p.add_argument('--pgap', default=None)
    p.add_argument('--gspa-key-map', default=None)
    p.add_argument('--tag', required=True)
    p.add_argument('--n-bootstrap', type=int, default=200)
    args = p.parse_args()

    acc_mapper = None
    if args.gspa_key_map:
        acc_mapper = {}
        with open(args.gspa_key_map) as fh:
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 2:
                    acc_mapper[parts[0]] = parts[1]

    gspa = load_gspa_integrated(args.gspa, acc_mapper=acc_mapper)
    gspa_priors = load_gspa_integrated(args.gspa_priors, acc_mapper=acc_mapper) if args.gspa_priors else None
    pgap = load_pgap_binary(args.pgap) if args.pgap else None

    out = {'tag': args.tag, 'by_truth': {}}

    for spec in args.truth:
        name, path = spec.split(':', 1)
        truth = load_truth(path)
        ic_map = compute_ic(path)
        total_annotations = sum(sum(len(terms) for terms in asps.values()) for asps in truth.values())
        block = {
            'truth_proteins': len(truth),
            'truth_annotations': total_annotations,
            'results': [summarize('GSPA', gspa, truth, ic_map, args.n_bootstrap)],
        }
        if gspa_priors is not None:
            block['results'].append(summarize('GSPA+priors', gspa_priors, truth, ic_map, args.n_bootstrap))
        if pgap is not None:
            block['results'].append(summarize('PGAP', pgap, truth, ic_map, args.n_bootstrap))
        out['by_truth'][name] = block

    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
