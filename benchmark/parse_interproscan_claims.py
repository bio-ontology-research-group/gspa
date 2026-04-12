#!/usr/bin/env python3
"""
Parse InterProScan TSV output into GSPA claims.jsonl.

InterProScan TSV columns (0-indexed):
  0: protein_id
  1: MD5
  2: length
  3: analysis (Pfam, TIGRFAM, CDD, SUPERFAMILY, etc.)
  4: signature accession
  5: signature description
  6: start
  7: stop
  8: e-value
  9: match status (T/F)
 10: date
 11: InterPro accession
 12: InterPro description
 13: GO annotations (pipe-separated, e.g. "GO:0003677(InterPro)|GO:0006260(InterPro)")
 14: Pathways (optional)

Emits one claim per (protein, GO term) pair with source='interproscan',
evidence_type='SEQUENCE_DOMAIN', and raw_score derived from the e-value.
"""
import argparse
import json
import math
import sys


def evalue_to_score(evalue_str):
    """Convert e-value to a [0,1] confidence score.

    InterProScan hits that pass the gathering threshold are
    high-confidence. Calibration:
      e-value ≤ 1e-50 → 0.95
      e-value   1e-20 → 0.85
      e-value   1e-5  → 0.60
      e-value   1     → 0.30
    """
    try:
        ev = float(evalue_str)
    except (ValueError, TypeError):
        return 0.7  # matched but no e-value → assume decent
    if ev <= 0:
        return 0.95
    neg_log = -math.log10(max(ev, 1e-300))
    # Sigmoid-like mapping: score = 0.95 / (1 + exp(-0.1*(neg_log - 10)))
    score = 0.95 / (1.0 + math.exp(-0.1 * (neg_log - 10.0)))
    return max(0.1, min(0.95, score))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input', required=True, help='InterProScan TSV')
    p.add_argument('--output', required=True, help='Output claims JSONL')
    args = p.parse_args()

    seen = set()
    n_claims = 0
    with open(args.input) as fin, open(args.output, 'w') as fout:
        for line in fin:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 14:
                continue
            protein_id = fields[0]
            analysis = fields[3]
            evalue = fields[8]
            go_field = fields[13] if len(fields) > 13 else ''

            if not go_field or go_field == '-':
                continue

            score = evalue_to_score(evalue)

            # Parse GO terms: "GO:0003677(InterPro)|GO:0006260(InterPro)"
            for go_entry in go_field.split('|'):
                go_entry = go_entry.strip()
                if not go_entry.startswith('GO:'):
                    continue
                go_id = go_entry.split('(')[0]
                key = (protein_id, go_id)
                if key in seen:
                    continue
                seen.add(key)

                claim = {
                    'protein_id': protein_id,
                    'function_type': 'GO',
                    'function_id': go_id,
                    'go_aspect': '',  # InterProScan doesn't provide aspect in TSV
                    'source': 'interproscan',
                    'raw_score': round(score, 4),
                    'metadata': {
                        'analysis': analysis,
                        'evalue': evalue,
                    }
                }
                fout.write(json.dumps(claim) + '\n')
                n_claims += 1

    print(f'  {n_claims:,} InterProScan claims from {len(seen):,} unique (protein, GO) pairs')


if __name__ == '__main__':
    main()
