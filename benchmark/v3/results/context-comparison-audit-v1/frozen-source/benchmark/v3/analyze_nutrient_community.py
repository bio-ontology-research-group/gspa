#!/usr/bin/env python3
"""Freeze and apply paired, dependence-aware descriptive coculture comparisons."""
import argparse
from copy import deepcopy
from collections import defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
import run_nutrient_community as B

ANALYSIS = {
    'schema': 'gspa.nutrient_community_analysis/1',
    'seed': 20261027, 'resamples': 2000, 'percentiles': [2.5, 97.5],
    'groups': ['species', 'medium_oxygen_arm', 'primary_or_no_nitrogen'],
    'bulk_methods': ['permissive', 'native_unsupplied', 'no_unsupplied_uptake', 'always_growth', 'never_growth'],
    'rescue_methods': ['permissive', 'always_growth', 'never_growth'],
    'paired_comparisons': {'bulk_growth': [
        ['permissive', 'native_unsupplied'], ['permissive', 'no_unsupplied_uptake'],
        *[[method, baseline] for method in B.PROTOCOL['policies'] for baseline in ['always_growth', 'never_growth']]],
        'paired_rescue_proxy': [['permissive', 'always_growth'], ['permissive', 'never_growth']]},
    'unknowns': 'Report all unknowns and coverage; compare accuracy on common complete cases only. Undefined resamples remain undefined.',
    'environment_resampling': 'Resample whole environments with replacement, preserving all pairs together; conditional on the fixed strain panel.',
    'environment_and_strain_resampling': 'Independently resample environments and strain vertices; each observed distinct pair receives environment multiplicity times both strain multiplicities. No synthetic same-strain pairs are created.',
    'leave_one_strain_out': 'Drop all pairs containing each strain; descriptive sensitivity without fitting or threshold adjustment.',
    'interpretation': 'Descriptive percentile ranges, not calibrated confidence intervals, significance tests, or performance on unseen strains. Six selected strains and four no-nitrogen environments especially limit inference.',
    'selection': 'Report every registered arm, method and comparison; no outcome-based arm selection or promotion.',
}


def source_hashes():
    return {'analysis': B.N.sha(Path(__file__)), 'execution': B.N.sha(Path(B.__file__))}


def catalogue(root):
    B.verify_plan(root)
    cases, identities = {}, set()
    for plan in B.N.read_rows(root / 'plans.jsonl'):
        identity = plan['case_id'], plan['arm']
        if identity in identities: raise ValueError('duplicate plan')
        identities.add(identity)
        record = {k: plan[k] for k in ('case_id', 'species', 'environment', 'primary_environment')}
        record['members'] = sorted(m['id'] for m in plan['members'])
        if len(set(record['members'])) != 2: raise ValueError('distinct strain pair required')
        previous = cases.setdefault(record['case_id'], record)
        if previous != record: raise ValueError('case metadata changes across arms')
    return [cases[k] for k in sorted(cases)], identities


def freeze(root, out):
    if out.exists(): raise ValueError('analysis freeze already exists')
    cases, identities = catalogue(root)
    record = {'schema': ANALYSIS['schema'], 'analysis': deepcopy(ANALYSIS), 'implementation': source_hashes(),
        'created_utc': datetime.now(timezone.utc).isoformat(), 'protocol_sha256': B.N.sha(root / 'protocol.json'),
        'catalogue_sha256': B.N.digest(cases), 'cases': len(cases), 'plans': len(identities),
        'inputs_read': ['frozen execution protocol and outcome-free plans'],
        'outcomes_or_predictions_read': False}
    B.N.write_json(out, record)
    return record


def validate_freeze(root, frozen_path):
    frozen = json.loads(frozen_path.read_text())
    cases, identities = catalogue(root)
    if (frozen['analysis'] != ANALYSIS or frozen['implementation'] != source_hashes()
            or frozen['protocol_sha256'] != B.N.sha(root / 'protocol.json')
            or frozen['catalogue_sha256'] != B.N.digest(cases) or frozen['plans'] != len(identities)
            or frozen['cases'] != len(cases) or frozen['outcomes_or_predictions_read'] is not False):
        raise ValueError('analysis freeze or execution identity mismatch')
    return cases, identities


def validate_observations(cases, observations):
    rows = B.N.read_rows(observations)
    by_id, mono = {}, set()
    for row in rows:
        if row['case_id'] in by_id: raise ValueError('duplicate observation identity')
        by_id[row['case_id']] = row
        if row['kind'] == 'monoculture':
            if len(row['strains']) != 1: raise ValueError('invalid monoculture identity')
            key = row['species'], row['environment'], row['strains'][0]
            if key in mono: raise ValueError('duplicate matched monoculture')
            mono.add(key)
    for case in cases:
        row = by_id.get(case['case_id'])
        if row is None: continue  # missing remains unknown, never a negative endpoint
        if (row['kind'] != 'coculture' or row['species'] != case['species'] or row['environment'] != case['environment']
                or sorted(row['species'] + ':' + s for s in row['strains']) != case['members']):
            raise ValueError('coculture observation identity mismatch')


def resampling_weights(rows, mode, replicates=None):
    """One weight vector applies to both methods in every paired comparison."""
    if mode not in ('environment', 'environment_and_strain'): raise ValueError('unknown resampling mode')
    n = ANALYSIS['resamples'] if replicates is None else replicates
    environments = sorted({r['environment'] for r in rows})
    strains = sorted({s for r in rows for s in r['members']})
    if len(environments) < 2 or (mode == 'environment_and_strain' and len(strains) < 3):
        return np.empty((0, len(rows))), environments, strains
    ei, si = {x: i for i, x in enumerate(environments)}, {x: i for i, x in enumerate(strains)}
    e = np.array([ei[r['environment']] for r in rows])
    a = np.array([si[r['members'][0]] for r in rows])
    b = np.array([si[r['members'][1]] for r in rows])
    rng = np.random.default_rng(ANALYSIS['seed'])
    weights = []
    for _ in range(n):
        ec = np.bincount(rng.integers(len(environments), size=len(environments)), minlength=len(environments))
        w = ec[e]
        if mode == 'environment_and_strain':
            sc = np.bincount(rng.integers(len(strains), size=len(strains)), minlength=len(strains))
            w = w * sc[a] * sc[b]
        weights.append(w)
    return np.asarray(weights), environments, strains


def interval(values, eligible, weights):
    point = float(np.sum(values) / np.sum(eligible)) if np.sum(eligible) else None
    if not len(weights): return {'estimate': point, 'percentiles': [None, None], 'defined_resamples': 0, 'resamples': 0}
    totals, denominators = weights @ values, weights @ eligible
    valid = denominators > 0
    draws = totals[valid] / denominators[valid]
    return {'estimate': point, 'percentiles': np.percentile(draws, ANALYSIS['percentiles']).tolist() if len(draws) else [None, None],
            'defined_resamples': int(valid.sum()), 'resamples': len(weights)}


def compare(rows, left, right, weights):
    valid, differences, corrections, harms = [], [], 0, 0
    unknown = dict(observation=0, left=0, right=0)
    for row in rows:
        y, a, b = row['observed'], row['methods'][left], row['methods'][right]
        if any(v is not None and not isinstance(v, bool) for v in (y, a, b)): raise ValueError('boolean or unknown required')
        unknown['observation'] += y is None
        unknown['left'] += a is None
        unknown['right'] += b is None
        known = all(v is not None for v in (y, a, b))
        valid.append(float(known))
        differences.append(float((a == y) - (b == y)) if known else 0.)
        corrections += known and a == y and b != y
        harms += known and a != y and b == y
    valid, differences = np.asarray(valid), np.asarray(differences)
    summary = {'left': left, 'right': right, 'n_total': len(rows), 'n_complete_pairs': int(valid.sum()),
        'unknown': unknown, 'corrections': int(corrections), 'new_errors': int(harms),
        'accuracy_difference': float(differences.sum() / valid.sum()) if valid.sum() else None,
        'resampling': {mode: interval(differences, valid, matrix) for mode, matrix in weights.items()}}
    strains = sorted({s for row in rows for s in row['members']})
    summary['leave_one_strain_out'] = []
    for strain in strains:
        retain = np.array([strain not in r['members'] for r in rows], dtype=float)
        n = float(retain @ valid)
        summary['leave_one_strain_out'].append({'excluded_strain': strain, 'n_complete_pairs': int(n),
            'accuracy_difference': float(retain @ differences / n) if n else None})
    return summary


def summarize(cases, base):
    by_id = {r['case_id']: r for r in cases}
    grouped, seen = defaultdict(list), set()
    for row in base['cases']:
        identity = row['case_id'], row['arm']
        if identity in seen: raise ValueError('duplicated scored case')
        seen.add(identity)
        meta = by_id[row['case_id']]
        group = meta['species'], row['arm'], 'primary' if meta['primary_environment'] else 'no_nitrogen_control'
        p = row['predictions']
        for endpoint in ('bulk_growth', 'paired_rescue_proxy'):
            methods = {'always_growth': True, 'never_growth': False}
            if endpoint == 'bulk_growth':
                methods.update({method: p[method] for method in ANALYSIS['bulk_methods'][:3]})
                observed = row['observed_bulk_growth']
            else:
                methods['permissive'] = B.rescue(p['permissive'], p['monoculture_0'], p['monoculture_1'])
                observed = row['observed_rescue_proxy']
            grouped[group + (endpoint,)].append({'case_id': row['case_id'], 'environment': meta['environment'],
                'members': meta['members'], 'observed': observed, 'methods': methods})
    output = []
    for group, raw in sorted(grouped.items()):
        rows = sorted(raw, key=lambda r: r['case_id'])
        weights = {}
        for mode in ('environment', 'environment_and_strain'):
            weights[mode], environments, strains = resampling_weights(rows, mode)
        methods = []
        for name in sorted(rows[0]['methods']):
            known = np.array([r['observed'] is not None and r['methods'][name] is not None for r in rows], dtype=float)
            correct = np.array([r['observed'] == r['methods'][name] for r in rows], dtype=float) * known
            methods.append({'method': name, **B.N.confusion([(r['observed'], r['methods'][name]) for r in rows]),
                'coverage': float(known.sum() / len(rows)),
                'resampling': {mode: interval(correct, known, matrix) for mode, matrix in weights.items()}})
        output.append({'species': group[0], 'arm': group[1], 'scope': group[2], 'endpoint': group[3],
            'n_environments': len(environments), 'n_strains': len(strains), 'methods': methods,
            'paired': [compare(rows, a, b, weights) for a, b in ANALYSIS['paired_comparisons'][group[3]]]})
    return output


def analyze(root, run, observations, frozen_path, out):
    cases, identities = validate_freeze(root, frozen_path)
    validate_observations(cases, observations)
    if out.exists(): raise ValueError('analysis output exists')
    out.mkdir(parents=True)
    base = B.evaluate(root, run, observations, out / 'endpoint-evaluation.json')
    if {(r['case_id'], r['arm']) for r in base['cases']} != identities or len(base['cases']) != len(identities):
        raise ValueError('scored catalogue incomplete')
    report = {'schema': 'gspa.nutrient_community_dependence_analysis/1', 'analysis_sha256': B.N.sha(frozen_path),
        'endpoint_evaluation_sha256': B.N.sha(out / 'endpoint-evaluation.json'), 'implementation': source_hashes(),
        'interpretation': ANALYSIS['interpretation'], 'groups': summarize(cases, base), 'promotion': 'none'}
    B.N.write_json(out / 'analysis.json', report)
    return report


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='command', required=True)
    for command, keys in [('freeze', ['plan', 'out']), ('analyze', ['plan', 'run', 'observations', 'analysis', 'out'])]:
        c = sub.add_parser(command)
        for k in keys: c.add_argument('--' + k, type=Path, required=True)
    a = p.parse_args()
    report = freeze(a.plan, a.out) if a.command == 'freeze' else analyze(a.plan, a.run, a.observations, a.analysis, a.out)
    print(json.dumps({'schema': report['schema'], 'groups': len(report.get('groups', []))}))
