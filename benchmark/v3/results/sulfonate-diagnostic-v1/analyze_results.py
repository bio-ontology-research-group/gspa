#!/usr/bin/env python3
"""Describe all supported responses without calling ordinal yields EC truth."""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def rows(path):
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def analyze(job):
    work = ROOT / 'jobs' / str(job) / 'results'
    plans = {(r['case_id'], r['arm']): r for r in rows(work / 'plan/plans.jsonl')}
    predictions = rows(work / 'run/predictions.jsonl')
    truth = {r['case_id']: r for r in rows(ROOT / 'dataset/observations.jsonl')}
    profiles = defaultdict(lambda: defaultdict(dict))
    disagreements, counts = [], defaultdict(Counter)
    for row in predictions:
        plan = plans[row['case_id'], row['arm']]
        counts[row['arm']][str(row['prediction'])] += 1
        if plan['status'] != 'ready':
            continue
        profiles[row['arm']][plan['genotype_label']][plan['sulfur_source']] = row['prediction']
        measured = truth[row['case_id']]
        if row['prediction'] is False and measured['symbol'] in ('+', '++', '+++'):
            disagreements.append({'case_id': row['case_id'], 'arm': row['arm'],
                'genotype': plan['genotype_label'], 'sulfur_source': plan['sulfur_source'],
                'published_symbol': measured['symbol'], 'source_cell': measured['source']})
    classes = {}
    for arm, profiles_by_genotype in sorted(profiles.items()):
        groups = defaultdict(list)
        for genotype, profile in profiles_by_genotype.items():
            groups[tuple(sorted(profile.items()))].append(genotype)
        classes[arm] = [{'genotypes': sorted(v), 'predictions': dict(k)} for k, v in sorted(groups.items())]
    initial = ROOT / 'jobs/10190/results'
    before = {(r['case_id'], r['arm']): r for r in rows(initial / 'run/predictions.jsonl')}
    def signature(row):
        return {'status': row['status'], 'prediction': row['prediction'],
                'records': [{k: record[k] for k in ('kind', 'status', 'outcome', 'measurements')}
                            for record in row.get('production', {}).get('experiments', [])]}
    return {'schema': 'gspa.sulfonate_summary/1', 'job': str(job),
        'inputs': {str(p.relative_to(ROOT)): sha(p) for p in (work / 'run/run.json', work / 'plan/protocol.json', ROOT / 'dataset/observations.jsonl')},
        'prediction_counts_by_arm': {k: dict(v) for k, v in sorted(counts.items())},
        'positive_yield_with_predicted_no_growth': disagreements,
        'genotype_response_classes': classes,
        'initial_replay_comparison': {'same_plans_bytes': sha(initial / 'plan/plans.jsonl') == sha(work / 'plan/plans.jsonl'),
            'same_case_population': set(before) == {(r['case_id'], r['arm']) for r in predictions},
            'changed_predictions_or_LP_records': sum(signature(r) != signature(before[r['case_id'], r['arm']]) for r in predictions)},
        'interpretation': 'Descriptive known-genotype diagnostic. Identical response patterns do not distinguish genotypes; neither model agreement nor growth feasibility establishes unknown protein function.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--job', required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():parser.error('output exists')
    args.out.write_text(json.dumps(analyze(args.job), sort_keys=True, indent=2) + '\n')
