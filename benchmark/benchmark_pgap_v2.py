#!/usr/bin/env python3
"""Patched benchmark script: per-genome micro-averaged F-max + CAFA-style protein-centric F-max.

Supports two annotation tracks selected via ``--annotation-type``:
 - ``GO`` (default, existing behaviour): per-aspect (MF / BP / CC) F-max
   in addition to an overall number.
 - ``EC``: a single overall F-max; rows are treated as
   {(protein_id, '', 'EC:1.1.1.1')} in the truth file.

The loaders both filter the integrated / PGAP inputs by the requested
type; the F-max routines already accept an ``aspect`` argument and are
type-agnostic once the loaders emit consistent aspect strings.
"""
import argparse
import json
import math
import random
from collections import defaultdict


# Aspects reported in the per-aspect output when --annotation-type=GO.
# EC runs collapse to a single overall number.
GO_ASPECTS = ('MF', 'BP', 'CC')


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


def load_gspa_integrated(path, acc_mapper=None, ann_type='GO'):
    """Stream an integrated GSPA TSV; keep rows whose ``type`` column matches
    ``ann_type``. For EC rows the ``go_aspect`` column will usually be blank,
    so we collapse all EC rows into the single synthetic aspect ``''``.
    """
    out = defaultdict(lambda: defaultdict(dict))
    with open(path) as fh:
        headers = next(fh, '').rstrip('\n').split('\t')
        idx_pid = headers.index('protein_id')
        idx_type = headers.index('type')
        idx_func = headers.index('function_id')
        idx_aspect = headers.index('go_aspect')
        idx_prob = headers.index('posterior_prob')
        default_aspect = 'BP' if ann_type == 'GO' else ''
        for line in fh:
            fields = line.rstrip('\n').split('\t')
            if len(fields) <= idx_prob or fields[idx_type] != ann_type:
                continue
            try:
                prob = float(fields[idx_prob])
            except ValueError:
                continue
            pid = fields[idx_pid]
            if acc_mapper is not None:
                pid = acc_mapper.get(pid, pid)
            if ann_type == 'GO':
                aspect = fields[idx_aspect] or default_aspect
            else:
                aspect = default_aspect
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


# ============================================================
# CAFA S-min: IC-weighted semantic distance.
#   ru_p(t) = sum IC(t) over t in truth_p \ pred_p(t)      (missed truth)
#   mi_p(t) = sum IC(t) over t in pred_p(t) \ truth_p      (wrong extras)
#   ru(t)   = mean ru_p over proteins with truth
#   mi(t)   = mean mi_p over proteins with truth
#   S(t)    = sqrt(ru(t)^2 + mi(t)^2)
#   S-min   = min over thresholds t
# Lower is better. 0 = perfect.
# ============================================================

def smin_cafa(predictions, truth, ic_map, aspect=None, thresholds=None):
    if thresholds is None:
        thresholds = [i * 0.05 for i in range(0, 21)]  # include 0.0 (recall-max)
    best_s = float('inf')
    best_t = 0.0
    best_ru = 0.0
    best_mi = 0.0
    for t in thresholds:
        rus, mis = [], []
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
                    if s >= t:
                        pred_set.add(g)
            ru = sum(ic_map.get(g, 0.0) for g in (true_set - pred_set))
            mi = sum(ic_map.get(g, 0.0) for g in (pred_set - true_set))
            rus.append(ru)
            mis.append(mi)
        if not rus:
            continue
        ru_bar = sum(rus) / len(rus)
        mi_bar = sum(mis) / len(mis)
        s = math.sqrt(ru_bar * ru_bar + mi_bar * mi_bar)
        if s < best_s:
            best_s = s
            best_t = t
            best_ru = ru_bar
            best_mi = mi_bar
    if best_s == float('inf'):
        return {'smin': 0.0, 'threshold': 0.0, 'ru': 0.0, 'mi': 0.0}
    return {'smin': best_s, 'threshold': best_t, 'ru': best_ru, 'mi': best_mi}


# ============================================================
# EC hierarchical evaluation.
# EC IDs like "EC:1.1.1.1" or "1.1.1.1". We split on '.' and truncate
# to 1..4 levels, ancestor-closing by prefix. Level-4 = exact 4-digit
# match; level-3 = class.sub.sub-sub; level-1 = class only.
# Used by load_truth/load_gspa_integrated outputs (aspect='' for EC).
# ============================================================

def _ec_strip_prefix(term):
    t = term
    if t.startswith('EC:'):
        t = t[3:]
    return t


def _ec_truncate(term, level):
    t = _ec_strip_prefix(term)
    parts = t.split('.')
    if len(parts) < level:
        return None  # can't ancestor to this level (partial EC)
    return '.'.join(parts[:level])


def ec_level_view(preds_or_truth, level, is_truth=False):
    """Project an EC preds/truth dict onto `level` (1..4) by prefix
    truncation. Predictions keep the max score per truncated term.
    Returns a new nested dict of the same shape."""
    out = defaultdict(lambda: defaultdict(dict if not is_truth else set))
    for acc, by_asp in preds_or_truth.items():
        for asp, entries in by_asp.items():
            if is_truth:
                for term in entries:
                    trunc = _ec_truncate(term, level)
                    if trunc is None:
                        continue
                    out[acc][asp].add(trunc)
            else:
                for term, score in entries.items():
                    trunc = _ec_truncate(term, level)
                    if trunc is None:
                        continue
                    if trunc not in out[acc][asp] or out[acc][asp][trunc] < score:
                        out[acc][asp][trunc] = score
    # convert inner defaults to plain structures
    return {a: {s: (dict(v) if not is_truth else set(v)) for s, v in by_asp.items()}
            for a, by_asp in out.items()}


def summarize_ec_levels(name, preds, truth, n_bootstrap):
    """Per-level EC F-max (levels 1..4). Level 4 = exact-4-digit match."""
    rows = {}
    for lvl in (1, 2, 3, 4):
        p_lvl = ec_level_view(preds, lvl, is_truth=False)
        t_lvl = ec_level_view(truth, lvl, is_truth=True)
        # fmax_with_ci expects truth as dict[acc][asp]->set; t_lvl already matches.
        f = fmax_with_ci(p_lvl, t_lvl, n_bootstrap=n_bootstrap)
        fc = fmax_cafa_with_ci(p_lvl, t_lvl, n_bootstrap=n_bootstrap)
        rows[f'ec_level{lvl}_fmax'] = f['fmax']
        rows[f'ec_level{lvl}_fmax_ci'] = [f['ci_low'], f['ci_high']]
        rows[f'ec_level{lvl}_fmax_cafa'] = fc['fmax']
        rows[f'ec_level{lvl}_fmax_cafa_ci'] = [fc['ci_low'], fc['ci_high']]
    return rows


def summarize(name, preds, truth, ic_map, n_bootstrap, ann_type='GO'):
    row = {'method': name}
    # Per-genome micro-averaged F-max (existing)
    f = fmax_with_ci(preds, truth, n_bootstrap=n_bootstrap)
    row['fmax_overall'] = f['fmax']
    row['fmax_ci'] = [f['ci_low'], f['ci_high']]
    row['fmax_t'] = f['threshold']
    if ann_type == 'GO':
        for asp in GO_ASPECTS:
            fa = fmax_with_ci(preds, truth, aspect=asp, n_bootstrap=n_bootstrap)
            row[f'fmax_{asp}'] = fa['fmax']
            row[f'fmax_{asp}_ci'] = [fa['ci_low'], fa['ci_high']]
    # CAFA-style protein-centric F-max
    fc = fmax_cafa_with_ci(preds, truth, n_bootstrap=n_bootstrap)
    row['fmax_cafa_overall'] = fc['fmax']
    row['fmax_cafa_ci'] = [fc['ci_low'], fc['ci_high']]
    row['fmax_cafa_t'] = fc['threshold']
    if ann_type == 'GO':
        for asp in GO_ASPECTS:
            fca = fmax_cafa_with_ci(preds, truth, aspect=asp, n_bootstrap=n_bootstrap)
            row[f'fmax_cafa_{asp}'] = fca['fmax']
            row[f'fmax_cafa_{asp}_ci'] = [fca['ci_low'], fca['ci_high']]
    row['coverage'] = coverage(preds, truth)
    row['unique_terms'] = unique_terms(preds)
    row['ic_recall'] = ic_recall(preds, truth, ic_map)
    # CAFA S-min (IC-weighted semantic distance)
    s = smin_cafa(preds, truth, ic_map)
    row['smin_overall'] = s['smin']
    row['smin_t'] = s['threshold']
    row['smin_ru'] = s['ru']
    row['smin_mi'] = s['mi']
    if ann_type == 'GO':
        for asp in GO_ASPECTS:
            sa = smin_cafa(preds, truth, ic_map, aspect=asp)
            row[f'smin_{asp}'] = sa['smin']
    if ann_type == 'EC':
        row.update(summarize_ec_levels(name, preds, truth, n_bootstrap))
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
    p.add_argument('--annotation-type', default='GO', choices=['GO', 'EC'],
                   help='Which annotation type to evaluate (default GO).')
    args = p.parse_args()

    acc_mapper = None
    if args.gspa_key_map:
        acc_mapper = {}
        with open(args.gspa_key_map) as fh:
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 2:
                    acc_mapper[parts[0]] = parts[1]

    gspa = load_gspa_integrated(args.gspa, acc_mapper=acc_mapper,
                                ann_type=args.annotation_type)
    gspa_priors = (load_gspa_integrated(args.gspa_priors, acc_mapper=acc_mapper,
                                         ann_type=args.annotation_type)
                   if args.gspa_priors else None)
    pgap = load_pgap_binary(args.pgap) if args.pgap else None

    out = {'tag': args.tag, 'annotation_type': args.annotation_type, 'by_truth': {}}

    for spec in args.truth:
        name, path = spec.split(':', 1)
        truth = load_truth(path)
        ic_map = compute_ic(path)
        total_annotations = sum(sum(len(terms) for terms in asps.values()) for asps in truth.values())
        block = {
            'truth_proteins': len(truth),
            'truth_annotations': total_annotations,
            'results': [summarize('GSPA', gspa, truth, ic_map, args.n_bootstrap,
                                  ann_type=args.annotation_type)],
        }
        if gspa_priors is not None:
            block['results'].append(summarize('GSPA+priors', gspa_priors, truth, ic_map,
                                              args.n_bootstrap,
                                              ann_type=args.annotation_type))
        if pgap is not None:
            block['results'].append(summarize('PGAP', pgap, truth, ic_map, args.n_bootstrap,
                                              ann_type=args.annotation_type))
        out['by_truth'][name] = block

    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
