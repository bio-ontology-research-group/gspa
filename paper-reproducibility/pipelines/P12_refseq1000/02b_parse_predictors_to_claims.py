#!/usr/bin/env python3
"""
GSPA Benchmark Step 2b: Parse predictor outputs into a unified claims.jsonl
that the Phase 7 integrator (gspa-cli integrate) can consume.

Runs once per benchmark. The BO driver then invokes gspa-cli integrate
repeatedly with different theta.json files but the same claims.jsonl,
which is far cheaper than re-running DIAMOND / eggNOG / Pfam per iteration.

Reads the output of 02_run_predictors.sh:
  diamond_results.tsv                      (DIAMOND blastp outfmt 6)
  mmseqs2_results.tsv                      (MMseqs2 convertalis)
  pfam_results.domtbl                      (hmmsearch --domtblout)
  eggnog_results.emapper.annotations       (emapper.py annotations)

Plus the Swiss-Prot/GOA annotations file (filtered to hit proteins) so
similarity hits can be transferred into GO terms.

Output: claims.jsonl — one JSON object per line.
"""

import argparse
import gzip
import json
import math
import os
import re
import sys
from collections import defaultdict


def collect_targets(*result_files):
    """Collect all UniProt accessions from similarity-based predictor outputs."""
    targets = set()
    for rf in result_files:
        if not rf or not os.path.exists(rf):
            continue
        with open(rf) as f:
            for line in f:
                fields = line.strip().split('\t', 2)
                if len(fields) < 2:
                    continue
                target = fields[1]
                acc = target.split('|')[1] if '|' in target else target
                targets.add(acc)
    return targets


def load_goa_for_targets(goa_file, target_accs):
    """Load GO annotations only for proteins in target_accs."""
    annotations = defaultdict(lambda: defaultdict(set))  # acc -> aspect -> GO terms
    if not target_accs:
        return annotations
    opener = gzip.open if goa_file.endswith('.gz') else open
    print(f"  Loading GO annotations for {len(target_accs):,} targets from {goa_file}...", flush=True)
    n_lines = 0
    n_kept = 0
    with opener(goa_file, 'rt') as f:
        for line in f:
            n_lines += 1
            if n_lines % 10_000_000 == 0:
                print(f"    {n_lines:,} lines scanned, {n_kept:,} kept", flush=True)
            if line.startswith('!'):
                continue
            fields = line.split('\t', 9)
            if len(fields) < 9:
                continue
            acc = fields[1]
            if acc not in target_accs:
                continue
            qualifier = fields[3]
            if 'NOT' in qualifier.upper():
                continue
            go_id = fields[4]
            aspect = fields[8]
            aspect_name = {'F': 'MF', 'P': 'BP', 'C': 'CC'}.get(aspect, aspect)
            annotations[acc][aspect_name].add(go_id)
            n_kept += 1
    print(f"  Total: {n_lines:,} lines scanned, {n_kept:,} GOA entries kept", flush=True)
    return annotations


def emit(out_fh, claim):
    """Write one claim as a JSON line."""
    out_fh.write(json.dumps(claim, separators=(',', ':')))
    out_fh.write('\n')


def parse_similarity(results_file, source, goa_ref, out_fh):
    """Parse DIAMOND / MMseqs2 outfmt6, transfer GO from hits."""
    if not os.path.exists(results_file):
        return 0
    n = 0
    with open(results_file) as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue
            query = fields[0]
            target = fields[1]
            try:
                pident = float(fields[2])
            except ValueError:
                continue
            acc = target.split('|')[1] if '|' in target else target
            raw = min(1.0, pident / 100.0)
            for aspect, terms in goa_ref.get(acc, {}).items():
                for term in terms:
                    emit(out_fh, {
                        'protein_id': query,
                        'function_type': 'GO',
                        'function_id': term,
                        'go_aspect': aspect,
                        'source': source,
                        'raw_score': raw,
                        'metadata': {'target': target, 'pident': pident},
                    })
                    n += 1
    return n


def parse_pfam(results_file, pfam2go_file, out_fh):
    """Parse HMMER --domtblout, map Pfam → GO via pfam2go mapping."""
    if not os.path.exists(results_file):
        return 0
    pfam2go = {}
    if pfam2go_file and os.path.exists(pfam2go_file):
        with open(pfam2go_file) as f:
            for line in f:
                if line.startswith('!'):
                    continue
                m = re.match(r'Pfam:(PF\d+)\s+\S+\s*>\s*GO:.*?;\s*(GO:\d+)', line)
                if m:
                    pfam2go.setdefault(m.group(1), set()).add(m.group(2))
        total_go = sum(len(v) for v in pfam2go.values())
        print(f"  Loaded {len(pfam2go):,} Pfam families with {total_go:,} GO mappings")
    n = 0
    with open(results_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.split()
            if len(fields) < 22:
                continue
            query = fields[0]
            pfam_acc = re.sub(r'\.\d+$', '', fields[4])
            try:
                dom_score = float(fields[13])
            except ValueError:
                continue
            raw = min(1.0, dom_score / 100.0)
            for go_term in pfam2go.get(pfam_acc, []):
                # pfam2go doesn't give an aspect. Use '?' so downstream can decide.
                for aspect in ['MF', 'BP', 'CC']:
                    emit(out_fh, {
                        'protein_id': query,
                        'function_type': 'GO',
                        'function_id': go_term,
                        'go_aspect': aspect,
                        'source': 'pfam',
                        'raw_score': raw,
                        'metadata': {'pfam': pfam_acc, 'dom_score': dom_score},
                    })
                    n += 1
    return n


def parse_eggnog(results_file, out_fh):
    """Parse eggNOG-mapper annotations for GO terms."""
    if not os.path.exists(results_file):
        return 0
    n = 0
    with open(results_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 10:
                continue
            query = fields[0]
            go_terms = fields[9] if len(fields) > 9 else '-'
            if go_terms == '-' or not go_terms:
                continue
            try:
                evalue = float(fields[2])
                raw = 1.0 if evalue <= 0 else min(1.0, max(0.4, -math.log10(evalue) / 50.0))
            except ValueError:
                raw = 0.7
            for go_term in go_terms.split(','):
                go_term = go_term.strip()
                if not go_term.startswith('GO:'):
                    continue
                for aspect in ['MF', 'BP', 'CC']:
                    emit(out_fh, {
                        'protein_id': query,
                        'function_type': 'GO',
                        'function_id': go_term,
                        'go_aspect': aspect,
                        'source': 'eggnog-mapper',
                        'raw_score': raw,
                        'metadata': {'evalue': evalue if 'evalue' in dir() else None},
                    })
                    n += 1
    return n


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results-dir', required=True,
                   help='Directory with 02_run_predictors.sh outputs')
    p.add_argument('--goa', required=True,
                   help='Swiss-Prot/GOA annotation file (gzipped ok)')
    p.add_argument('--pfam2go', default=None,
                   help='Pfam → GO mapping (external2go format)')
    p.add_argument('--output', required=True,
                   help='Output claims.jsonl path')
    p.add_argument('--test-accs', default=None,
                   help='Optional: file with test protein accessions (one per line) to exclude '
                        'from similarity targets so we never transfer GO from the query itself.')
    args = p.parse_args()

    diamond_file = os.path.join(args.results_dir, 'diamond_results.tsv')
    mmseqs_file = os.path.join(args.results_dir, 'mmseqs2_results.tsv')
    pfam_file = os.path.join(args.results_dir, 'pfam_results.domtbl')
    eggnog_file = os.path.join(args.results_dir, 'eggnog_results.emapper.annotations')

    print('Collecting similarity targets…', flush=True)
    targets = collect_targets(diamond_file, mmseqs_file)
    if args.test_accs and os.path.exists(args.test_accs):
        with open(args.test_accs) as f:
            exclude = {line.strip() for line in f if line.strip()}
        targets -= exclude
    print(f'  {len(targets):,} unique target accessions', flush=True)

    print('Loading GOA for targets…', flush=True)
    goa_ref = load_goa_for_targets(args.goa, targets)
    print(f'  {len(goa_ref):,} targets have GO annotations', flush=True)

    os.makedirs(os.path.dirname(os.path.abspath(args.output)) or '.', exist_ok=True)
    n_total = 0
    with open(args.output, 'w') as out_fh:
        if os.path.exists(diamond_file):
            print('Parsing DIAMOND…', flush=True)
            n = parse_similarity(diamond_file, 'diamond', goa_ref, out_fh)
            print(f'  {n:,} DIAMOND claims', flush=True)
            n_total += n
        if os.path.exists(mmseqs_file):
            print('Parsing MMseqs2…', flush=True)
            n = parse_similarity(mmseqs_file, 'mmseqs2', goa_ref, out_fh)
            print(f'  {n:,} MMseqs2 claims', flush=True)
            n_total += n
        if os.path.exists(pfam_file):
            print('Parsing HMMER/Pfam…', flush=True)
            n = parse_pfam(pfam_file, args.pfam2go, out_fh)
            print(f'  {n:,} Pfam claims', flush=True)
            n_total += n
        if os.path.exists(eggnog_file):
            print('Parsing eggNOG-mapper…', flush=True)
            n = parse_eggnog(eggnog_file, out_fh)
            print(f'  {n:,} eggNOG claims', flush=True)
            n_total += n

    print(f'Wrote {n_total:,} claims to {args.output}', flush=True)
    if n_total == 0:
        print('ERROR: no claims produced. Check predictor result files exist.', file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
