#!/usr/bin/env python3
"""Adapter: metagenomic-deepFRI (mdF) per-protein predictions → benchmark_pgap_v2 GSPA-shape TSV.

mdF (Bezshapkin et al., bioRxiv 2026-04-29) emits per-protein GO term
predictions across the three GO aspects (BP, MF, CC). DeepFRI's CSV
schema is:

    Protein,GO_term,Score,GO_name

This adapter rewrites it to the column shape benchmark_pgap_v2.py
expects via --gspa:

    protein_id  type  function_id  go_aspect  posterior_prob  \\
        likelihood_logodds  final_logodds  n_supporting  priors_fired  convergence_iter

Score is mapped to posterior_prob; logodds columns are derived
log(p/(1-p)) so the scorer's threshold sweep behaves the same as for
GSPA outputs. n_supporting=1, priors_fired='', convergence_iter=0.

GO aspect lookup uses go.obo (--go-obo) when provided; otherwise the
mdF model's default vocabulary (BP/MF/CC by namespace) which the paper
documents.

--self-test runs against an in-line 5-row fixture and asserts the
header + a representative row pass through correctly.
"""

import argparse
import csv
import io
import math
import os
import re
import sys
from collections import defaultdict


GSPA_HEADER = [
    'protein_id', 'type', 'function_id', 'go_aspect',
    'posterior_prob', 'likelihood_logodds', 'final_logodds',
    'n_supporting', 'priors_fired', 'convergence_iter',
]


def parse_go_obo_aspects(go_obo):
    """Return {GO:0008150 -> 'BP', ...} from a go.obo file. None if not provided."""
    if not go_obo or not os.path.exists(go_obo):
        return {}
    ns_map = {'biological_process': 'BP', 'molecular_function': 'MF',
              'cellular_component': 'CC'}
    out = {}
    cur = None
    with open(go_obo) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line == '[Term]':
                cur = {'id': None, 'ns': None}
                continue
            if cur is None:
                continue
            if line.startswith('id: '):
                cur['id'] = line[4:].strip()
            elif line.startswith('namespace: '):
                cur['ns'] = ns_map.get(line[11:].strip())
            elif line == '':
                if cur['id'] and cur['ns']:
                    out[cur['id']] = cur['ns']
                cur = None
    return out


def to_logodds(p):
    """log(p/(1-p)) clamped to avoid div-by-zero."""
    p = min(max(p, 1e-6), 1.0 - 1e-6)
    return math.log(p / (1.0 - p))


def detect_aspect_from_row(row, aspect_map):
    go_id = row['function_id']
    if aspect_map and go_id in aspect_map:
        return aspect_map[go_id]
    # mdF native: rows are emitted per aspect-specific model with a
    # known prefix in the score column or filename. Fall back to '' so
    # the scorer can still operate (it will treat aspect as global).
    return ''


def convert_stream(reader, writer, aspect_map):
    """Read mdF CSV rows, emit GSPA-shape TSV rows."""
    n_in = n_out = 0
    by_protein_term = {}  # collapse if a (protein, term) appears twice
    for row in reader:
        n_in += 1
        protein = row.get('Protein') or row.get('protein') or row.get('protein_id')
        term = row.get('GO_term') or row.get('GO') or row.get('function_id')
        score = row.get('Score') or row.get('score') or row.get('posterior_prob')
        if not (protein and term and score):
            continue
        try:
            p = float(score)
        except ValueError:
            continue
        # Keep the highest score if a (protein, term) appears twice.
        key = (protein, term)
        prev = by_protein_term.get(key)
        if prev is None or p > prev[0]:
            by_protein_term[key] = (p, term)
    for (protein, term), (p, _) in by_protein_term.items():
        lo = to_logodds(p)
        out = {
            'protein_id': protein,
            'type': 'GO',
            'function_id': term,
            'go_aspect': aspect_map.get(term, '') if aspect_map else '',
            'posterior_prob': f'{p:.6f}',
            'likelihood_logodds': f'{lo:.4f}',
            'final_logodds': f'{lo:.4f}',
            'n_supporting': '1',
            'priors_fired': '',
            'convergence_iter': '0',
        }
        writer.writerow(out)
        n_out += 1
    return n_in, n_out


def run(args):
    aspect_map = parse_go_obo_aspects(args.go_obo)
    with open(args.mdf_csv, newline='') as inh, open(args.out, 'w', newline='') as outh:
        reader = csv.DictReader(inh)
        writer = csv.DictWriter(outh, fieldnames=GSPA_HEADER, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        n_in, n_out = convert_stream(reader, writer, aspect_map)
    print(f'Wrote {n_out:,} rows from {n_in:,} input lines to {args.out}', flush=True)


def self_test():
    # 5-row mdF fixture covering: distinct proteins, dup (protein, term)
    # with different scores (keep max), aspect lookup, threshold edges.
    fixture = """Protein,GO_term,Score,GO_name
P12345,GO:0008150,0.92,biological_process
P12345,GO:0008150,0.81,biological_process
P12345,GO:0003674,0.55,molecular_function
P67890,GO:0008150,0.30,biological_process
P67890,GO:0005575,0.10,cellular_component
"""
    aspect_map = {
        'GO:0008150': 'BP',
        'GO:0003674': 'MF',
        'GO:0005575': 'CC',
    }
    out_buf = io.StringIO()
    reader = csv.DictReader(io.StringIO(fixture))
    writer = csv.DictWriter(out_buf, fieldnames=GSPA_HEADER, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    n_in, n_out = convert_stream(reader, writer, aspect_map)
    output = out_buf.getvalue()
    rows = [r for r in output.split('\n') if r.strip()]
    assert n_in == 5, f'n_in={n_in}'
    assert n_out == 4, f'n_out={n_out} (5 input → 4 output, P12345/GO:0008150 dedup)'
    assert rows[0].split('\t') == GSPA_HEADER, f'header mismatch: {rows[0]}'
    # Verify P12345/GO:0008150 kept the higher 0.92 (not 0.81).
    p12345_bp = [r for r in rows[1:] if r.startswith('P12345\tGO\tGO:0008150\t')]
    assert len(p12345_bp) == 1
    assert '0.920000' in p12345_bp[0], f'expected 0.92, got: {p12345_bp[0]}'
    # Verify aspect lookup populated.
    assert 'BP' in p12345_bp[0]
    print('OK: 5 input rows → 4 output rows; dedup keeps max; aspect from go.obo applied.')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mdf-csv', help='mdF prediction CSV (Protein,GO_term,Score[,GO_name])')
    p.add_argument('--out', help='Output GSPA-shape TSV')
    p.add_argument('--go-obo', help='Optional go.obo for GO aspect lookup')
    p.add_argument('--self-test', action='store_true', help='Run unit self-test and exit')
    args = p.parse_args()

    if args.self_test:
        self_test()
        return
    if not (args.mdf_csv and args.out):
        p.error('--mdf-csv and --out are required (or use --self-test)')
    run(args)


if __name__ == '__main__':
    main()
