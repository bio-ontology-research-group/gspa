#!/usr/bin/env python3
"""Verify retained evidence and computations without running a model or solver."""
import collections
import gzip
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def obj(path):
    return json.loads(path.read_text())


def rows(path):
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def check(job='10192', source='source-v2', control_job='10192'):
    manifest = dict(line.split('  ', 1)[::-1] for line in (ROOT / 'SHA256SUMS').read_text().splitlines())
    files = {str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if p.is_file() and p.name != 'SHA256SUMS'}
    assert files == set(manifest), 'archive file population differs'
    for name, expected in manifest.items():
        assert sha(ROOT / name) == expected, name
    spec = importlib.util.spec_from_file_location('frozen_sulfonate', ROOT / 'source/benchmark/v3/prepare_sulfonate_benchmark.py')
    extraction = importlib.util.module_from_spec(spec);spec.loader.exec_module(extraction)
    inputs, observations, report = extraction.extract(ROOT / 'sources/publisher-web-extraction.json', ROOT / 'sources/publisher-methods-web-extraction.json')
    assert inputs == rows(ROOT / 'dataset/inputs.jsonl')
    assert observations == rows(ROOT / 'dataset/observations.jsonl')
    assert (report['reported_cells'], report['missing_cells'], report['uncertain_genotype_cases']) == (247, 13, 20)
    work = ROOT / 'jobs' / job / 'results'
    protocol = obj(work / 'plan/protocol.json')
    for name, expected in protocol['implementation'].items():
        assert sha(ROOT / source / name) == expected, name
    assert protocol['inputs_sha256'] == sha(ROOT / 'dataset/inputs.jsonl')
    if job != '10190':assert protocol['observations_sha256'] == sha(ROOT / 'dataset/observations.jsonl')
    assert hashlib.sha256(gzip.decompress((work / 'iML1515.xml.gz').read_bytes())).hexdigest() == protocol['model_sha256']
    plans = rows(work / 'plan/plans.jsonl')
    assert len(plans) == 780 and digest(plans) == protocol['plans_sha256']
    keyed = {(p['case_id'], p['arm']): p for p in plans}
    assert len(keyed) == 780
    assert collections.Counter(p['status'] for p in plans) == {'ready': 216, 'unsupported': 564}
    run = obj(work / 'run/run.json')
    assert run['protocol_sha256'] == sha(work / 'plan/protocol.json')
    assert run['predictions_sha256'] == sha(work / 'run/predictions.jsonl')
    assert run['model_sha256'] == protocol['model_sha256']
    predictions = rows(work / 'run/predictions.jsonl')
    predicted = {(p['case_id'], p['arm']): p for p in predictions}
    assert len(predictions) == len(predicted) == 780 and predicted.keys() == keyed.keys()
    oracle = obj(work / 'oracle.json')
    assert oracle['run_sha256'] == sha(work / 'run/run.json')
    assert oracle['model_sha256'] == protocol['model_sha256']
    oracles = {(r['case_id'], r['arm']): r for r in oracle['results']}
    assert len(oracle['results']) == len(oracles) == 780 and oracles.keys() == keyed.keys()
    count = 0
    for key, plan in keyed.items():
        row, other = predicted[key], oracles[key]
        assert row['plan_sha256'] == digest(plan)
        if plan['status'] == 'unsupported':
            assert row['status'] == 'unsupported' and row['prediction'] is None
            assert other['agrees'] is None and other['checks'] == []
            continue
        assert row['status'] == 'executed' and row['prediction'] is not None
        records = row['production']['experiments']
        kinds = ['baseline', 'knockout'] if plan['genes'] else ['baseline']
        assert [r['kind'] for r in records] == kinds
        calls = row['production'].get('solver_log', row['production'].get('solver', {}).get('calls'))
        assert len(calls) == len(kinds)
        assert all(call['status'] in ('optimal', 'infeasible') for call in calls)
        assert [r['kind'] for r in other['checks']] == kinds
        assert other['agrees'] is True
        for record, comparison in zip(records, other['checks']):
            count += 1
            assert record['status'] == 'completed'
            expected = comparison['oracle']
            if record['outcome'] == 'infeasible_with_maintenance':
                assert expected['status'] == 'infeasible'
            else:
                assert expected['status'] == 'optimal'
                assert math.isclose(record['measurements']['max_growth'], expected['growth'], rel_tol=1e-6, abs_tol=1e-6)
                assert expected['max_balance_residual'] <= 1e-6
                assert expected['max_bound_violation'] <= 1e-6
            assert record['measurements'].get('knocked_out_reactions') == expected['blocked_reactions']
            assert comparison['flux_agrees'] is True and comparison['gpr_agrees'] is True
        target = records[-1]
        expected_growth = target['outcome'] != 'infeasible_with_maintenance' and target['measurements']['max_growth'] >= 1e-6
        assert row['prediction'] == expected_growth
    assert count == oracle['comparisons'] == 414 and oracle['disagreements'] == 0
    truth = {r['case_id']: r for r in observations}
    cross = collections.Counter((r['arm'], truth[r['case_id']]['observation_role'], truth[r['case_id']]['symbol'] or 'unreported', str(r['prediction'])) + ((r['status'],) if job != '10190' else ()) for r in predictions)
    evaluation = obj(work / 'evaluation.json')
    assert evaluation['observations_sha256'] == sha(ROOT / 'dataset/observations.jsonl')
    assert evaluation['run_sha256'] == sha(work / 'run/run.json')
    reported_cross = {(r['arm'], r['role'], r['observed_symbol'], r['predicted_growth']) + ((r['prediction_status'],) if job != '10190' else ()): r['count'] for r in evaluation['cross_tabulation']}
    assert dict(cross) == reported_cross
    assert len(evaluation['records']) == 780
    controls = obj(ROOT / 'jobs' / control_job / 'results/medium-controls.json')
    assert controls['implementation_sha256'] == sha(ROOT / source / 'benchmark/v3/check_sulfonate_medium.py')
    assert controls['protocol_sha256'] == sha(work / 'plan/protocol.json')
    assert controls['model_sha256'] == protocol['model_sha256']
    assert controls['passes'] and len(controls['results']) == 2
    for control in controls['results']:
        assert control['growth_detected'] is False and control['agrees'] is True
        assert control['plan']['manifest']['baseline']['medium']['R_EX_so4_e'] == 0
        assert control['oracle']['status'] == 'optimal' and abs(control['oracle']['growth']) < 1e-6
        record = control['production']['production']['experiments'][0]
        assert record['status'] == 'completed' and abs(record['measurements']['max_growth']) < 1e-6
    if job != '10190':
        assert run['status_counts'] == {'executed': 216, 'unsupported': 564} and run['execution_complete']
        assert evaluation['execution_complete'] and evaluation['status_counts'] == run['status_counts']
        analysis_spec = importlib.util.spec_from_file_location('analysis', ROOT / 'analyze_results.py')
        analysis = importlib.util.module_from_spec(analysis_spec);analysis_spec.loader.exec_module(analysis)
        summary = obj(ROOT / 'summary.json')
        assert analysis.analyze(job) == summary
        assert summary['initial_replay_comparison'] == {'same_plans_bytes': True, 'same_case_population': True, 'changed_predictions_or_LP_records': 0}
        assert all(len(groups) == 4 for groups in summary['genotype_response_classes'].values())
        assert len(summary['positive_yield_with_predicted_no_growth']) == 18
    print(json.dumps({'job': job, 'files_verified': len(manifest), 'study_cases': len(inputs), 'reported_cells': 247,
                      'plans_verified': 780, 'production_oracle_LP_pairs_verified': count, 'medium_controls_verified': 2}))


if __name__ == '__main__':
    check('10190', 'source', '10191')
    check()
