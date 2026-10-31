"""Fault injection for complete populations, honest failures and frozen labels.

No model or solver runs in these tests; the real-model oracle is a separate gate.
"""
import copy
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import run_sulfonate_benchmark as S


def plan(case, status='ready'):
    return {'case_id': case, 'arm': 'test', 'status': status, 'genes': [],
            'reasons': [] if status == 'ready' else ['sulfur_exchange_unrepresented']}


def prediction(entry, status):
    row = {'case_id': entry['case_id'], 'arm': entry['arm'], 'plan_sha256': S.N.digest(entry),
           'status': status, 'prediction': None}
    if status == 'unsupported':
        row['reasons'] = entry['reasons']
    elif status == 'failed':
        row['error'] = 'injected production failure'
    else:
        row['production'] = {'experiments': [{'kind': 'baseline',
            'status': 'completed' if status == 'executed' else 'solver_error',
            'outcome': 'test', 'measurements': {'max_growth': 1.0}}]}
        if status == 'executed':row['prediction'] = True
    return row


class SulfonateRunnerTest(unittest.TestCase):
    def setUp(self):
        base = Path.home() / 'tmp';base.mkdir(exist_ok=True)
        self.temp = tempfile.TemporaryDirectory(dir=base)
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.plans = [plan('supported'), plan('unsupported', 'unsupported')]
        self.rows = [prediction(self.plans[0], 'executed'), prediction(self.plans[1], 'unsupported')]
        self.truth = [{'case_id': p['case_id'], 'symbol': '+', 'observation_role': 'deletion_response',
                       'genotype_status': 'table_reported'} for p in self.plans]
        self.plan_dir, self.run_dir = self.root / 'plan', self.root / 'run'
        self.plan_dir.mkdir();self.run_dir.mkdir()
        self.observations = self.root / 'observations.jsonl'
        self.freeze()

    def freeze(self):
        def write(path, rows):
            path.write_text(''.join(json.dumps(r, sort_keys=True) + '\n' for r in rows))
        write(self.observations, self.truth)
        write(self.plan_dir / 'plans.jsonl', self.plans)
        S.N.write_json(self.plan_dir / 'protocol.json', {'protocol': S.PROTOCOL, 'implementation': S.implementation(),
            'model_sha256': '0' * 64, 'plans_sha256': S.N.digest(self.plans),
            'observations_sha256': S.N.sha(self.observations)})
        write(self.run_dir / 'predictions.jsonl', self.rows)
        counts = dict(S.Counter(r['status'] for r in self.rows))
        S.N.write_json(self.run_dir / 'run.json', {'model_sha256': '0' * 64,
            'protocol_sha256': S.N.sha(self.plan_dir / 'protocol.json'),
            'predictions_sha256': S.N.sha(self.run_dir / 'predictions.jsonl'),
            'status_counts': counts, 'execution_complete': not any(counts.get(s, 0) for s in ('failed', 'inconclusive'))})

    def test_missing_duplicate_and_wrong_plan_results_rejected(self):
        original = copy.deepcopy(self.rows)
        fixtures = [original[:1], [original[0], original[0]], copy.deepcopy(original)]
        fixtures[-1][0]['plan_sha256'] = 'bad'
        for fixture in fixtures:
            with self.subTest(fixture=fixture):
                self.rows = fixture;self.freeze()
                with self.assertRaises(ValueError):S.verify_run(self.plan_dir, self.run_dir)

    def test_unsupported_cannot_silently_acquire_prediction_or_execution(self):
        good = self.rows[1]
        fixtures = []
        for key, value in [('prediction', False), ('production', {'experiments': []}), ('status', 'executed')]:
            bad = copy.deepcopy(good);bad[key] = value;fixtures.append(bad)
        for bad in fixtures:
            with self.assertRaises(ValueError):S.validate_prediction(self.plans[1], bad)

    def test_incomplete_solver_result_cannot_be_called_executed(self):
        row = prediction(self.plans[0], 'inconclusive')
        S.validate_prediction(self.plans[0], row)
        row['status'] = 'executed'
        with self.assertRaises(ValueError):S.validate_prediction(self.plans[0], row)

    def test_missing_experiment_and_fabricated_binary_prediction_rejected(self):
        bad = copy.deepcopy(self.rows[0]);bad['production']['experiments'] = []
        with self.assertRaises(ValueError):S.validate_prediction(self.plans[0], bad)
        bad = copy.deepcopy(self.rows[0]);bad['prediction'] = False
        with self.assertRaises(ValueError):S.validate_prediction(self.plans[0], bad)

    def test_observations_changed_after_registration_rejected(self):
        self.observations.write_text(self.observations.read_text().replace('"+"', '"−"'))
        with self.assertRaisesRegex(ValueError, 'observations differ'):
            S.evaluate(self.plan_dir, self.run_dir, self.observations, self.root / 'evaluation.json')

    def test_unsupported_failed_and_inconclusive_are_separate_in_summary(self):
        self.plans = [plan('unsupported', 'unsupported'), plan('failed'), plan('inconclusive')]
        self.rows = [prediction(p, p['case_id']) for p in self.plans]
        self.truth = [{'case_id': p['case_id'], 'symbol': '+', 'observation_role': 'deletion_response',
                       'genotype_status': 'table_reported'} for p in self.plans]
        self.freeze()
        out = self.root / 'evaluation.json'
        S.evaluate(self.plan_dir, self.run_dir, self.observations, out)
        report = json.loads(out.read_text())
        self.assertFalse(report['execution_complete'])
        self.assertEqual({r['prediction_status'] for r in report['cross_tabulation']}, {'failed', 'inconclusive', 'unsupported'})
        self.assertEqual(len(report['records']), 3)

    def test_run_preserves_all_cases_and_reports_failure(self):
        with patch.object(S, 'verify', return_value=({}, self.plans)), \
             patch.object(S.ME, 'load_sbml_model', return_value=SimpleNamespace(sha256='0' * 64)), \
             patch.object(S.N, 'simulate', side_effect=RuntimeError('injected')):
            report = S.run(self.plan_dir, self.root / 'model.xml', self.root / 'failed_run')
        self.assertFalse(report['execution_complete'])
        self.assertEqual(report['status_counts'], {'failed': 1, 'unsupported': 1})
        records = S.N.read_rows(self.root / 'failed_run/predictions.jsonl')
        self.assertEqual(len(records), 2)
        self.assertTrue(all(r['prediction'] is None for r in records))

    def test_failed_ready_case_cannot_pass_oracle_with_zero_comparisons(self):
        import cobra
        self.rows[0] = prediction(self.plans[0], 'failed');self.freeze()
        model_path = self.root / 'model.xml';model_path.write_text('not parsed')
        frozen, plans, predictions = S.verify_run(self.plan_dir, self.run_dir)
        frozen['model_sha256'] = S.N.sha(model_path)
        with patch.object(S, 'verify_run', return_value=(frozen, plans, predictions)), \
             patch.object(cobra.io, 'read_sbml_model', return_value=SimpleNamespace(solver=None)):
            report = S.oracle(self.plan_dir, self.run_dir, model_path, self.root / 'oracle.json')
        self.assertEqual(report['comparisons'], 0)
        self.assertEqual(report['disagreements'], 1)
        self.assertEqual(report['unsupported'], 1)


if __name__ == '__main__':unittest.main()
