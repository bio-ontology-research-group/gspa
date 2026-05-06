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


_MDF_ASPECT_FROM_MODE = {
    'GO Biological Process': 'BP',
    'GO Molecular Function': 'MF',
    'GO Cellular Component': 'CC',
    # 'Enzyme Commission' rows carry EC numbers, not GO terms — skipped
    # below because benchmark_pgap_v2.py reads --gspa as GO-only here.
}


def convert_stream(reader, writer, aspect_map):
    """Read mdF results.tsv rows, emit GSPA-shape TSV rows.

    Accepts both the canonical mdF v1.0 schema:
      protein  network_type  prediction_mode  go_term  score  go_name  ...
    and the legacy DeepFRI CSV schema:
      Protein  GO_term  Score  GO_name
    """
    n_in = n_out = n_skipped = 0
    by_protein_term = {}  # collapse if a (protein, term) appears twice
    for row in reader:
        n_in += 1
        protein = row.get('protein') or row.get('Protein') or row.get('protein_id')
        term = row.get('go_term') or row.get('GO_term') or row.get('GO') or row.get('function_id')
        score = row.get('score') or row.get('Score') or row.get('posterior_prob')
        mode = row.get('prediction_mode')
        # Skip EC predictions — benchmark_pgap_v2.py via --gspa expects GO only.
        if mode and mode == 'Enzyme Commission':
            n_skipped += 1
            continue
        if not (protein and term and score):
            continue
        try:
            p = float(score)
        except ValueError:
            continue
        # Resolve aspect: prefer mdF prediction_mode, then go.obo lookup.
        aspect = _MDF_ASPECT_FROM_MODE.get(mode or '')
        if not aspect and aspect_map:
            aspect = aspect_map.get(term, '')
        # Keep the highest score if a (protein, term) appears twice.
        key = (protein, term)
        prev = by_protein_term.get(key)
        if prev is None or p > prev[0]:
            by_protein_term[key] = (p, aspect or '')
    for (protein, term), (p, aspect) in by_protein_term.items():
        lo = to_logodds(p)
        out = {
            'protein_id': protein,
            'type': 'GO',
            'function_id': term,
            'go_aspect': aspect,
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


def _sniff_delimiter(path):
    """Return '\\t' for tab-separated input, ',' for comma. Defaults to tab
    for mdF results.tsv (the canonical mdF v1.0 output)."""
    with open(path) as fh:
        first = fh.readline()
    return '\t' if first.count('\t') > first.count(',') else ','


def run(args):
    aspect_map = parse_go_obo_aspects(args.go_obo)
    delim = _sniff_delimiter(args.mdf_csv)
    with open(args.mdf_csv, newline='') as inh, open(args.out, 'w', newline='') as outh:
        reader = csv.DictReader(inh, delimiter=delim)
        writer = csv.DictWriter(outh, fieldnames=GSPA_HEADER, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        n_in, n_out = convert_stream(reader, writer, aspect_map)
    print(f'Wrote {n_out:,} rows from {n_in:,} input lines to {args.out}', flush=True)


def self_test():
    # Fixture covers the canonical mdF v1.0 results.tsv schema:
    # - aspect inferred from prediction_mode column
    # - EC rows skipped
    # - dup (protein, term) keeps higher score
    fixture = (
        "protein\tnetwork_type\tprediction_mode\tgo_term\tscore\tgo_name\n"
        "P12345\tcnn\tGO Biological Process\tGO:0008150\t0.92\tbiological_process\n"
        "P12345\tcnn\tGO Biological Process\tGO:0008150\t0.81\tbiological_process\n"
        "P12345\tcnn\tGO Molecular Function\tGO:0003674\t0.55\tmolecular_function\n"
        "P67890\tcnn\tGO Biological Process\tGO:0008150\t0.30\tbiological_process\n"
        "P67890\tcnn\tGO Cellular Component\tGO:0005575\t0.10\tcellular_component\n"
        "P67890\tcnn\tEnzyme Commission\tEC:1.1.1.1\t0.40\talcohol_dehydrogenase\n"
    )
    out_buf = io.StringIO()
    reader = csv.DictReader(io.StringIO(fixture), delimiter='\t')
    writer = csv.DictWriter(out_buf, fieldnames=GSPA_HEADER, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    n_in, n_out = convert_stream(reader, writer, aspect_map={})
    output = out_buf.getvalue()
    rows = [r for r in output.split('\n') if r.strip()]
    assert n_in == 6, f'n_in={n_in}'
    assert n_out == 4, f'n_out={n_out} (6 input -> 4 output: 1 EC skipped + 1 dup deduped)'
    assert rows[0].split('\t') == GSPA_HEADER, f'header mismatch: {rows[0]}'
    # P12345/GO:0008150 kept the higher 0.92.
    p12345_bp = [r for r in rows[1:] if r.startswith('P12345\tGO\tGO:0008150\t')]
    assert len(p12345_bp) == 1
    assert '0.920000' in p12345_bp[0]
    # Aspect inferred from prediction_mode.
    assert '\tBP\t' in p12345_bp[0], f'expected BP aspect, got: {p12345_bp[0]}'
    # No EC rows in output.
    assert not any('EC:' in r for r in rows[1:])
    print('OK: 6 mdF input rows -> 4 output rows; dedup keeps max; aspect from prediction_mode; EC skipped.')


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
