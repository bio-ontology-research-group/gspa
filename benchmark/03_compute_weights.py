#!/usr/bin/env python3
"""
GSPA Benchmark Step 3: Compute optimal combination weights.

Memory-efficient: only loads GO annotations for proteins referenced by predictors.
"""

import argparse
import gzip
import json
import os
import re
import sys
import math
from collections import defaultdict
from itertools import product as iterproduct

def collect_targets(*result_files):
    """Collect all UniProt accessions from similarity-based predictor outputs."""
    targets = set()
    for rf in result_files:
        if not os.path.exists(rf): continue
        with open(rf) as f:
            for line in f:
                fields = line.strip().split('\t', 2)
                if len(fields) < 2: continue
                target = fields[1]
                acc = target.split('|')[1] if '|' in target else target
                targets.add(acc)
    return targets

def load_goa_for_targets(goa_file, target_accs):
    """Load GO annotations only for proteins in target_accs."""
    annotations = defaultdict(lambda: defaultdict(set))
    if not target_accs:
        return annotations
    opener = gzip.open if goa_file.endswith('.gz') else open
    print(f"  Loading GO annotations for {len(target_accs)} targets...")
    n_lines = 0
    n_kept = 0
    with opener(goa_file, 'rt') as f:
        for line in f:
            n_lines += 1
            if n_lines % 10000000 == 0:
                print(f"    {n_lines:,} lines scanned, {n_kept:,} kept")
            if line.startswith('!'): continue
            fields = line.split('\t', 9)
            if len(fields) < 9: continue
            acc = fields[1]
            if acc not in target_accs: continue
            qualifier = fields[3]
            if 'NOT' in qualifier.upper(): continue
            go_id = fields[4]
            aspect = fields[8]
            aspect_name = {'F': 'MF', 'P': 'BP', 'C': 'CC'}.get(aspect, aspect)
            annotations[acc][aspect_name].add(go_id)
            n_kept += 1
    print(f"  Total: {n_lines:,} lines scanned, {n_kept:,} annotations kept")
    return annotations

def load_ground_truth(truth_file):
    truth = defaultdict(lambda: defaultdict(set))
    with open(truth_file) as f:
        next(f)
        for line in f:
            acc, aspect, go_term = line.strip().split('\t')
            truth[acc][aspect].add(go_term)
    return truth

def parse_similarity_results(results_file, goa_ref):
    """Parse DIAMOND or MMseqs2 output, transfer GO from hits."""
    predictions = defaultdict(lambda: defaultdict(dict))
    with open(results_file) as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 3: continue
            query = fields[0]
            target = fields[1]
            try:
                pident = float(fields[2])
            except ValueError:
                continue
            acc = target.split('|')[1] if '|' in target else target
            score = min(1.0, pident / 100.0)
            for aspect, terms in goa_ref.get(acc, {}).items():
                for term in terms:
                    cur = predictions[query][aspect].get(term, 0)
                    if score > cur:
                        predictions[query][aspect][term] = score
    return predictions

def parse_pfam_results(results_file, pfam2go_file=None):
    """Parse HMMER domtblout, map Pfam to GO via pfam2go."""
    pfam2go = {}
    if pfam2go_file and os.path.exists(pfam2go_file):
        with open(pfam2go_file) as f:
            for line in f:
                if line.startswith('!'): continue
                m = re.match(r'Pfam:(PF\d+)\s+\S+\s*>\s*GO:.*?;\s*(GO:\d+)', line)
                if m:
                    pfam2go.setdefault(m.group(1), set()).add(m.group(2))
        total_go = sum(len(v) for v in pfam2go.values())
        print(f"  Loaded {len(pfam2go)} Pfam families with {total_go} GO mappings")

    predictions = defaultdict(lambda: defaultdict(dict))
    with open(results_file) as f:
        for line in f:
            if line.startswith('#'): continue
            fields = line.split()
            if len(fields) < 22: continue
            query = fields[0]
            pfam_acc = re.sub(r'\.\d+$', '', fields[4])
            try:
                dom_score = float(fields[13])
            except ValueError:
                continue
            score = min(1.0, dom_score / 100.0)
            for go_term in pfam2go.get(pfam_acc, []):
                # Pfam2GO doesn't give aspect, put in all — will be filtered by ground truth
                for aspect in ['MF', 'BP', 'CC']:
                    cur = predictions[query][aspect].get(go_term, 0)
                    if score > cur:
                        predictions[query][aspect][go_term] = score
    return predictions

def parse_eggnog_results(results_file):
    """Parse eggNOG-mapper annotations for GO terms."""
    predictions = defaultdict(lambda: defaultdict(dict))
    with open(results_file) as f:
        for line in f:
            if line.startswith('#'): continue
            fields = line.strip().split('\t')
            if len(fields) < 10: continue
            query = fields[0]
            go_terms = fields[9] if len(fields) > 9 else '-'
            if go_terms == '-' or not go_terms: continue
            try:
                evalue = float(fields[2])
                if evalue > 0:
                    score = min(1.0, max(0.4, -math.log10(evalue) / 50.0))
                else:
                    score = 1.0
            except ValueError:
                score = 0.7
            for go_term in go_terms.split(','):
                go_term = go_term.strip()
                if go_term.startswith('GO:'):
                    for aspect in ['MF', 'BP', 'CC']:
                        cur = predictions[query][aspect].get(go_term, 0)
                        if score > cur:
                            predictions[query][aspect][go_term] = score
    return predictions

def compute_metrics(predictions, truth, threshold=0.0, aspect_filter=None):
    tp = fp = fn = 0
    for acc in truth:
        for aspect in truth[acc]:
            if aspect_filter and aspect != aspect_filter: continue
            true_terms = truth[acc][aspect]
            pred_terms = {t for t, s in predictions.get(acc, {}).get(aspect, {}).items() if s >= threshold}
            tp += len(true_terms & pred_terms)
            fp += len(pred_terms - true_terms)
            fn += len(true_terms - pred_terms)
    p = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    r = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * p * r / (p + r) if (p + r) > 0 else 0.0
    return {'precision': p, 'recall': r, 'f1': f1, 'tp': tp, 'fp': fp, 'fn': fn}

def compute_fmax(predictions, truth, aspect_filter=None):
    thresholds = [i * 0.05 for i in range(21)]
    best_f1 = 0
    best_t = 0
    for t in thresholds:
        m = compute_metrics(predictions, truth, t, aspect_filter)
        if m['f1'] > best_f1:
            best_f1 = m['f1']
            best_t = t
    return best_f1, best_t

def combine_predictions(pred_list, weights):
    combined = defaultdict(lambda: defaultdict(dict))
    for preds, w in zip(pred_list, weights):
        for acc in preds:
            for aspect in preds[acc]:
                for term, score in preds[acc][aspect].items():
                    weighted = score * w
                    cur = combined[acc][aspect].get(term, 0)
                    if weighted > cur:
                        combined[acc][aspect][term] = weighted
    return combined

def grid_search_weights(pred_list, truth, step=0.2, aspect_filter=None):
    n = len(pred_list)
    vals = [round(i * step, 2) for i in range(int(1.0/step) + 1)]
    best_fmax = 0
    best_weights = [1.0 / n] * n
    for combo in iterproduct(vals, repeat=n):
        total = sum(combo)
        if total < 0.1: continue
        weights = [w / total for w in combo]
        combined = combine_predictions(pred_list, weights)
        fmax, _ = compute_fmax(combined, truth, aspect_filter)
        if fmax > best_fmax:
            best_fmax = fmax
            best_weights = weights
    return best_weights, best_fmax

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--results-dir', required=True)
    p.add_argument('--data-dir', required=True)
    p.add_argument('--goa', required=True)
    p.add_argument('--pfam2go', default=None)
    p.add_argument('--output', default='weights.json')
    p.add_argument('--weight-step', type=float, default=0.2)
    args = p.parse_args()

    print("Loading ground truth...")
    truth = load_ground_truth(os.path.join(args.data_dir, 'test_go_annotations.tsv'))
    test_accs = set(truth.keys())
    print(f"  {len(test_accs)} test proteins")

    diamond_file = os.path.join(args.results_dir, 'diamond_results.tsv')
    mmseqs_file = os.path.join(args.results_dir, 'mmseqs2_results.tsv')
    pfam_file = os.path.join(args.results_dir, 'pfam_results.domtbl')
    eggnog_file = os.path.join(args.results_dir, 'eggnog_results.emapper.annotations')

    print("\nPass 1: Collecting target accessions from similarity searches...")
    target_accs = collect_targets(diamond_file, mmseqs_file)
    target_accs -= test_accs
    print(f"  {len(target_accs)} unique targets")

    print("\nPass 2: Loading GO annotations from GOA (filtered by targets)...")
    goa_ref = load_goa_for_targets(args.goa, target_accs)
    print(f"  {len(goa_ref)} targets have GO annotations")

    predictions = {}
    tool_names = []

    if os.path.exists(diamond_file) and os.path.getsize(diamond_file) > 0:
        print("\nParsing DIAMOND results...")
        predictions['diamond'] = parse_similarity_results(diamond_file, goa_ref)
        tool_names.append('diamond')
        print(f"  {len(predictions['diamond'])} proteins predicted")

    if os.path.exists(pfam_file) and os.path.getsize(pfam_file) > 0:
        print("\nParsing Pfam results...")
        predictions['pfam'] = parse_pfam_results(pfam_file, args.pfam2go)
        tool_names.append('pfam')
        print(f"  {len(predictions['pfam'])} proteins predicted")

    if os.path.exists(mmseqs_file) and os.path.getsize(mmseqs_file) > 0:
        print("\nParsing MMseqs2 results...")
        predictions['mmseqs2'] = parse_similarity_results(mmseqs_file, goa_ref)
        tool_names.append('mmseqs2')
        print(f"  {len(predictions['mmseqs2'])} proteins predicted")

    if os.path.exists(eggnog_file) and os.path.getsize(eggnog_file) > 0:
        print("\nParsing eggNOG results...")
        predictions['eggnog'] = parse_eggnog_results(eggnog_file)
        tool_names.append('eggnog')
        print(f"  {len(predictions['eggnog'])} proteins predicted")

    if not tool_names:
        print("ERROR: No predictor results found!")
        sys.exit(1)

    del goa_ref

    print("\n=== Per-tool overall performance ===")
    results = {}
    for name in tool_names:
        fmax, t = compute_fmax(predictions[name], truth)
        m = compute_metrics(predictions[name], truth, t)
        results[name] = {'fmax': fmax, 'best_threshold': t,
                         'precision': m['precision'], 'recall': m['recall'], 'f1': m['f1']}
        print(f"  {name:10s}  F-max: {fmax:.4f}  P: {m['precision']:.4f}  R: {m['recall']:.4f}  (t={t:.2f})")

    print("\n=== Per-aspect performance ===")
    per_aspect = {}
    for aspect in ['MF', 'BP', 'CC']:
        print(f"\n  --- {aspect} ---")
        per_aspect[aspect] = {}
        for name in tool_names:
            fmax, t = compute_fmax(predictions[name], truth, aspect_filter=aspect)
            per_aspect[aspect][name] = {'fmax': fmax, 'threshold': t}
            print(f"    {name:10s}  F-max: {fmax:.4f}  (t={t:.2f})")

    print(f"\n=== Weight optimization (step={args.weight_step}) ===")
    pred_list = [predictions[n] for n in tool_names]
    best_w, best_f = grid_search_weights(pred_list, truth, args.weight_step)
    weight_dict = {n: round(w, 3) for n, w in zip(tool_names, best_w)}
    print(f"\nOverall optimal weights: {weight_dict}")
    print(f"Overall combined F-max: {best_f:.4f}")

    per_aspect_weights = {}
    for aspect in ['MF', 'BP', 'CC']:
        print(f"\n  Optimizing for {aspect}...")
        bw, bf = grid_search_weights(pred_list, truth, args.weight_step, aspect_filter=aspect)
        per_aspect_weights[aspect] = {
            'weights': {n: round(w, 3) for n, w in zip(tool_names, bw)},
            'combined_fmax': bf,
        }
        print(f"    weights: {per_aspect_weights[aspect]['weights']}")
        print(f"    F-max: {bf:.4f}")

    individual_best = max(results[n]['fmax'] for n in tool_names)
    individual_best_name = max(tool_names, key=lambda n: results[n]['fmax'])
    improvement = (best_f - individual_best) / individual_best * 100 if individual_best > 0 else 0
    print(f"\nBest individual: {individual_best_name} ({individual_best:.4f})")
    print(f"Combined improvement: {improvement:+.1f}%")

    output = {
        'optimal_weights_overall': weight_dict,
        'combined_fmax_overall': best_f,
        'per_tool': results,
        'per_aspect_fmax': per_aspect,
        'per_aspect_weights': per_aspect_weights,
        'improvement_over_best_individual': f"{improvement:+.1f}%",
        'n_test_proteins': len(test_accs),
        'tools_evaluated': tool_names,
    }
    with open(args.output, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved: {args.output}")

if __name__ == '__main__':
    main()
