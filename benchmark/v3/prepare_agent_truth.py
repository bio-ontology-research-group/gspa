#!/usr/bin/env python3
"""Join experimental Swiss-Prot positives to already-frozen public cases.

This is a private evaluation phase. Missing assertions do not become negatives.
GO ancestry propagates positives only; no exhaustive locus truth is claimed.
"""
import argparse
import csv
from pathlib import Path
from collections import defaultdict, Counter

import evaluate_stage3_agent as evaluation
from v3_metrics import load_obo, GO_ROOTS
from v3_runner import hash_record


def prepare(public, labels_path, go_path):
    cases = evaluation.public_cases(public)
    split_of = {c['genome_id']: c['split'] for c in cases}
    # Include zero-case genomes in the evidence record, too.
    split_of.update({r['genome_id']: r['split'] for r in public.get('coverage', [])})
    relevant = defaultdict(set)
    onto = load_obo(str(go_path))
    for c in cases:
        a = c['anomaly']
        target = a['target_function_id']
        terms = onto.ancestors(target) - GO_ROOTS if target.startswith('GO:') else {target}
        for candidate in a['candidates']:
            relevant[c['genome_id'], candidate['protein_id']].update(terms)
    by_split = {'calibration': [], 'test': []}
    counts = Counter()
    with Path(labels_path).open() as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            gid, pid, fid = row['genome_id'], row['protein_id'], row['function_id']
            counts['input_rows'] += 1
            if gid not in split_of or row['tier'] != 'experimental':
                continue
            counts['experimental_rows_in_evaluation_genomes'] += 1
            wanted = relevant[gid, pid]
            terms = onto.ancestors(fid) - GO_ROOTS if fid.startswith('GO:') else {fid}
            for term in sorted(terms & wanted):
                by_split[split_of[gid]].append({'genome_id': gid, 'protein_id': pid,
                    'function_id': term, 'status': 'positive',
                    'source': {'accession': row['swissprot_accession'], 'asserted_function': fid,
                               'evidence': row['evidence'], 'match_basis': row['match_basis'],
                               'propagated': term != fid}})
    return by_split, {'schema': 'gspa.agent-truth-provenance/2',
                     'public_sha256': evaluation.digest(public),
                     'inputs': [hash_record(str(labels_path)), hash_record(str(go_path))],
                     'counts': dict(counts), 'policy': 'experimental positives only; every unlabelled pair stays unknown'}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--public', required=True, type=Path)
    p.add_argument('--labels', required=True, type=Path)
    p.add_argument('--go-obo', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    a = p.parse_args()
    labels, provenance = prepare(evaluation.read(a.public), a.labels, a.go_obo)
    for split, rows in labels.items():
        evaluation.write(a.out / f'{split}.json', rows)
    evaluation.write(a.out / 'provenance.json', provenance)
