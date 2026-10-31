"""Paired denominators, shared-environment/strain dependence, and frozen analysis identity."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import analyze_nutrient_community as A
from test_nutrient_community import fixture, case


def rows():
    return [dict(case_id=f'{e}:{a}:{b}', environment=e, members=[a, b], observed=True,
                 methods={'left': e == 'E1', 'right': e == 'E2'})
            for e in ('E1', 'E2') for a, b in [('a', 'b'), ('a', 'c'), ('b', 'c')]]


class CommunityAnalysisTests(unittest.TestCase):
    def test_same_environment_pairs_share_bootstrap_weight_and_repeated_rows_do_not_create_clusters(self):
        raw = rows()
        weights, environments, strains = A.resampling_weights(raw, 'environment', 200)
        self.assertEqual((environments, strains), (['E1', 'E2'], ['a', 'b', 'c']))
        np.testing.assert_array_equal(weights[:, 0], weights[:, 1])
        np.testing.assert_array_equal(weights[:, 1], weights[:, 2])
        np.testing.assert_array_equal(weights.sum(axis=1), np.full(200, 6))
        duplicated = [r for r in raw for _ in range(10)]
        w2, _, _ = A.resampling_weights(duplicated, 'environment', 200)
        first = A.compare(raw, 'left', 'right', {'environment': weights})
        second = A.compare(duplicated, 'left', 'right', {'environment': w2})
        self.assertEqual(first['resampling'], second['resampling'])

    def test_strain_weights_are_shared_across_both_pair_positions_and_zero_edge_draws_stay_undefined(self):
        raw = rows()
        weights, _, _ = A.resampling_weights(raw, 'environment_and_strain', 500)
        reverse = [dict(r, members=r['members'][::-1]) for r in raw]
        w2, _, _ = A.resampling_weights(reverse, 'environment_and_strain', 500)
        np.testing.assert_array_equal(weights, w2)
        self.assertTrue(np.any(weights.sum(axis=1) == 0))
        result = A.compare(raw, 'left', 'right', {'crossed': weights})
        self.assertLess(result['resampling']['crossed']['defined_resamples'], 500)
        self.assertEqual(len(result['leave_one_strain_out']), 3)
        self.assertTrue(all(r['n_complete_pairs'] == 2 for r in result['leave_one_strain_out']))

    def test_paired_comparison_uses_common_cases_and_reports_corrections_harm_and_unknowns(self):
        raw = rows()[:4]
        raw[0]['methods'] = dict(left=True, right=False)
        raw[1]['methods'] = dict(left=False, right=True)
        raw[2]['methods'] = dict(left=None, right=True)
        raw[3]['observed'] = None
        result = A.compare(raw, 'left', 'right', {'identity': np.ones((1, len(raw)))})
        self.assertEqual((result['n_complete_pairs'], result['corrections'], result['new_errors']), (2, 1, 1))
        self.assertEqual(result['accuracy_difference'], 0.)
        self.assertEqual(result['unknown'], dict(observation=1, left=1, right=0))
        for r in raw: r['methods']['left'] = None
        result = A.compare(raw, 'left', 'right', {'identity': np.ones((1, len(raw)))})
        self.assertIsNone(result['accuracy_difference'])
        self.assertEqual(result['resampling']['identity']['percentiles'], [None, None])

    def test_sparse_design_has_no_invented_interval_and_comparison_sign_reverses(self):
        raw = rows()[:3]
        w, _, _ = A.resampling_weights(raw, 'environment')
        self.assertEqual(len(w), 0)
        result = A.compare(raw, 'left', 'right', {'environment': w})
        self.assertEqual(result['accuracy_difference'], 1.)
        self.assertEqual(result['resampling']['environment']['percentiles'], [None, None])
        opposite = A.compare(raw, 'right', 'left', {'environment': w})
        self.assertEqual(opposite['accuracy_difference'], -1.)
        self.assertEqual(opposite['new_errors'], result['corrections'])

    def test_scored_order_does_not_change_ranges(self):
        cases, scored = [], []
        for r in rows():
            cases.append({k: r[k] for k in ('case_id', 'environment', 'members')} | dict(species='Ec', primary_environment=True))
            scored.append(dict(case_id=r['case_id'], arm='fixed', observed_bulk_growth=True,
                observed_rescue_proxy=True, predictions=dict(monoculture_0=False, monoculture_1=False,
                permissive=r['methods']['left'], native_unsupplied=r['methods']['right'], no_unsupplied_uptake=False)))
        self.assertEqual(A.summarize(cases, {'cases': scored}), A.summarize(cases[::-1], {'cases': scored[::-1]}))

    def test_identity_mismatches_and_duplicate_matched_monocultures_fail(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'obs.jsonl'
            expected = [dict(case_id='p', species='Ec', environment='e', members=['Ec:A', 'Ec:B'])]
            obs = [dict(case_id='p', species='Ec', environment='e', kind='coculture', strains=['A', 'C'])]
            p.write_text(json.dumps(obs[0]) + '\n')
            with self.assertRaisesRegex(ValueError, 'identity'): A.validate_observations(expected, p)
            obs = [dict(case_id=str(i), species='Ec', environment='e', kind='monoculture', strains=['A']) for i in range(2)]
            p.write_text(''.join(json.dumps(r) + '\n' for r in obs))
            with self.assertRaisesRegex(ValueError, 'duplicate matched'): A.validate_observations(expected, p)

    def test_frozen_analysis_requires_complete_verified_execution_but_no_outcomes_at_freeze(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            model, meta, context = fixture(root / 'toy.xml')
            inputs = root / 'inputs.jsonl'
            inputs.write_text(json.dumps(case()) + '\n')
            registration = root / 'registration.json'
            A.B.N.write_json(registration, {'inputs_sha256': A.B.N.sha(inputs), 'models': {'Ec': meta}})
            A.B.prepare(inputs, registration, root / 'plan')
            frozen = A.freeze(root / 'plan', root / 'analysis-plan.json')
            self.assertFalse(frozen['outcomes_or_predictions_read'])
            with self.assertRaisesRegex(ValueError, 'already exists'): A.freeze(root / 'plan', root / 'analysis-plan.json')
            A.B.run(root / 'plan', root / 'run', {'Ec': root / 'toy.xml'})
            obs = root / 'observations.jsonl'
            obs.write_text(json.dumps(dict(case_id='test', species='Ec', environment='test medium', kind='coculture',
                                           strains=['V', 'M'], mean_od600=.1)) + '\n')
            report = A.analyze(root / 'plan', root / 'run', obs, root / 'analysis-plan.json', root / 'out')
            self.assertEqual(len(report['groups']), 12)
            self.assertEqual(report['promotion'], 'none')
            rescue = next(r for r in report['groups'] if r['endpoint'] == 'paired_rescue_proxy')
            self.assertTrue(all(r['n_complete_pairs'] == 0 for r in rescue['paired']))
            frozen['analysis']['seed'] += 1
            self.assertNotEqual(frozen['analysis']['seed'], A.ANALYSIS['seed'])
            A.B.N.write_json(root / 'altered.json', frozen)
            with self.assertRaisesRegex(ValueError, 'identity'): A.validate_freeze(root / 'plan', root / 'altered.json')


if __name__ == '__main__':
    unittest.main()
