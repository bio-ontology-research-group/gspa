"""Truth separation and actual-change scoring, not biological validation."""
import copy
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import evaluate_stage3_agent as E


class Guard:
    def check_go_term(self, term, taxon):
        return {'valid': term != 'GO:0005634'}

    def go_ancestors(self, term):
        return {'GO:0000002'}


def bundle():
    cases = []
    for split, gid in [('calibration', 'G1'), ('test', 'G2')]:
        cases.append({'split': split, 'genome_id': gid, 'anomaly': {
            'anomaly_id': 'A_' + gid, 'target_function_id': 'GO:0000001',
            'organism_taxon': 'NCBITaxon_2', 'candidates': [
                {'protein_id': 'p1', 'baseline_probability': .4, 'raw_sequence_score': .6,
                 'model_probabilities': {'learned_logistic': .8, 'learned_gbdt': .7}},
                {'protein_id': 'p2', 'baseline_probability': .3, 'raw_sequence_score': .5,
                 'model_probabilities': {'learned_logistic': .2, 'learned_gbdt': .3}}]}})
    return {'schema': 'gspa.agent-public/2', 'cases': cases,
            'shared_production_wall_seconds': 12, 'limitations': []}


def truth(gid, status='positive', pid='p1', fid='GO:0000001'):
    return [{'genome_id': gid, 'protein_id': pid, 'function_id': fid,
             'status': status, 'source': 'fixture:independent-assay'}]


class EvaluationTests(unittest.TestCase):
    def setUp(self):
        self.public = bundle()
        self.guard = Guard()

    def proposals(self, mode='rank_by_raw_seq'):
        return E.suggest(self.public, mode, self.guard)

    def test_missing_truth_stays_unknown(self):
        p = self.proposals()
        policy = E.calibrate(self.public, p, truth('G1'), min_decisions=1)
        report = E.evaluate(self.public, p, policy, [])
        self.assertEqual(report['counts'], {'unknown': 1})
        self.assertIsNone(report['known_decision_precision'])
        self.assertEqual(report['demonstrated_precision_lower_bound'], 0)

    def test_explicit_truth_distinguishes_correct_and_harmful(self):
        p = self.proposals()
        policy = E.calibrate(self.public, p, truth('G1'), min_decisions=1)
        self.assertEqual(E.evaluate(self.public, p, policy, truth('G2'))['counts'], {'correct': 1})
        self.assertEqual(E.evaluate(self.public, p, policy, truth('G2', 'negative'))['counts'], {'harmful': 1})

    def test_unknown_calibration_cannot_create_confident_policy(self):
        p = self.proposals()
        policy = E.calibrate(self.public, p, [], min_decisions=1)
        self.assertIsNone(policy['threshold'])
        self.assertEqual(E.evaluate(self.public, p, policy, truth('G2'))['counts'], {'abstained': 1})

    def test_test_truth_rejected_by_calibration(self):
        with self.assertRaisesRegex(ValueError, 'another split'):
            E.calibrate(self.public, self.proposals(), truth('G2'), min_decisions=1)

    def test_nested_truth_rejected_before_agent(self):
        self.public['cases'][0]['anomaly']['context'] = {'nested': [{'is_true_target': True}]}
        with self.assertRaisesRegex(ValueError, 'truth fields'):
            E.suggest(self.public, 'agent', self.guard, solve=lambda a: self.fail('leaked input'))

    def test_identical_guardrails_for_all_rankers_and_agent(self):
        for c in self.public['cases']:
            c['anomaly']['target_function_id'] = 'GO:0005634'
        for mode in E.MODES:
            self.assertTrue(all(r['action'] == 'ABSTAIN' for r in self.proposals(mode)['resolutions']))
        def solve(a):
            return {'anomaly_id': a['anomaly_id'], 'action': 'PROMOTE', 'protein_id': 'p1',
                    'function_id': a['target_function_id'], 'confidence': .9}
        self.assertTrue(all(r['action'] == 'ABSTAIN' for r in E.suggest(self.public, 'agent', self.guard, solve)['resolutions']))

    def test_missing_duplicate_or_modified_predictions_rejected(self):
        p = self.proposals()
        policy = E.calibrate(self.public, p, truth('G1'), min_decisions=1)
        for mutate in (lambda x: x['resolutions'].pop(),
                       lambda x: x['resolutions'].append(x['resolutions'][0]),
                       lambda x: x.update(public_sha256='wrong')):
            q = copy.deepcopy(p)
            mutate(q)
            with self.assertRaises(ValueError):
                E.evaluate(self.public, q, policy, truth('G2'))

    def test_policy_bound_to_calibration_predictions(self):
        p = self.proposals()
        policy = E.calibrate(self.public, p, truth('G1'), min_decisions=1)
        p['resolutions'][0]['confidence'] = .99
        with self.assertRaisesRegex(ValueError, 'changed after'):
            E.evaluate(self.public, p, policy, truth('G2'))

    def test_already_called_function_not_a_recovery(self):
        c = self.public['cases'][1]
        c['anomaly']['candidates'][0]['baseline_probability'] = .8
        r = E.rank('rank_by_raw_seq', c['anomaly'])
        t = E.truth_map(truth('G2'), self.public['cases'], 'test')
        self.assertEqual(E.outcome(c, r, t), 'unchanged')

    def test_demoting_true_function_harmful_and_unlabelled_unknown(self):
        c = self.public['cases'][1]
        c['anomaly']['candidates'][0]['baseline_probability'] = .8
        r = E.rank('rank_by_raw_seq', c['anomaly'])
        r['action'] = 'DEMOTE'
        t = E.truth_map(truth('G2'), self.public['cases'], 'test')
        self.assertEqual(E.outcome(c, r, t), 'harmful')
        self.assertEqual(E.outcome(c, r, {}), 'unknown')

    def test_relax_requires_original_negative_and_ancestor_positive(self):
        c = self.public['cases'][1]
        c['anomaly']['candidates'][0]['baseline_probability'] = .8
        r = E.rank('rank_by_raw_seq', c['anomaly'])
        r.update(action='RELAX', function_id='GO:0000002')
        labels = truth('G2', fid='GO:0000002')
        self.assertEqual(E.outcome(c, r, E.truth_map(labels, self.public['cases'], 'test')), 'unknown')
        labels += truth('G2', 'negative')
        self.assertEqual(E.outcome(c, r, E.truth_map(labels, self.public['cases'], 'test')), 'correct')

    def test_missing_model_score_not_invented(self):
        a = self.public['cases'][0]['anomaly']
        for c in a['candidates']:
            c['model_probabilities']['learned_logistic'] = None
        self.assertEqual(E.rank('learned_logistic', a)['action'], 'ABSTAIN')

    def test_agent_error_survives_threshold(self):
        p = self.proposals()
        p['resolutions'][1] = E.abstain(self.public['cases'][1]['anomaly'], guardrail_status='ERROR')
        policy = E.calibrate(self.public, p, [], min_decisions=1)
        self.assertEqual(E.evaluate(self.public, p, policy, truth('G2'))['counts'], {'error': 1})

    def test_artifacts_immutable(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'policy.json'
            E.write(path, {'value': 1})
            with self.assertRaises(FileExistsError):
                E.write(path, {'value': 2})
            self.assertEqual(E.read(path), {'value': 1})

    def test_legacy_runs_cannot_supply_cases(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'run_record.json'
            E.write(path, {'status': 'executed'})
            with self.assertRaisesRegex(ValueError, 'input identity'):
                E.verified_run(path)

    def test_live_adapter_sees_public_cases_and_counts_reported_cost(self):
        def call(**kwargs):
            text = kwargs['messages'][1]['content']
            self.assertNotIn('expected_protein_id', text)
            self.assertIn('model_probabilities', text)
            return {'usage': {'prompt_tokens': 20, 'completion_tokens': 4},
                    'choices': [{'message': {'content': '{"action":"PROMOTE","protein_id":"p1","function_id":"GO:0000001","confidence":0.8}'}}]}
        with patch.object(E.agent, 'call_chat_completion', side_effect=call), patch.object(E.agent, 'load_api_key', return_value='fixture'):
            p = E.suggest_agent(self.public, self.guard, api_base='http://fixture.invalid/v1', workers=1)
        self.assertEqual(p['api_usage']['completed_calls'], 2)
        self.assertEqual(p['api_usage']['prompt_tokens'], 40)
        self.assertEqual(p['api_usage']['completion_tokens'], 8)
        self.assertTrue(all(r['action'] == 'PROMOTE' for r in p['resolutions']))
        self.assertIsNone(p['billing_cost'])


if __name__ == '__main__':
    unittest.main()
