"""Separate genotypes, bounded pools, uptake controls, unknowns and outcome isolation."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / 'scripts/tests'))
import run_nutrient_community as B
import check_nutrient_community_cobra as O
from test_community_metabolism import write_sbml


def fixture(path, reverse=False):
    reactions = {
        'R_EX_glc__D_e': ({'glc_e': 1}, {}, -1000, 1000, None),
        'R_EX_vit_e': ({'vit_e': 1}, {}, 0, 1000, None),
        'R_EX_met__L_e': ({'met_e': 1}, {}, 0, 1000, None),
        'T_glc': ({'glc_e': 1}, {'pyr_c': 4, 'atp_c': 8}, 0, 1000, None),
        'VITS': ({'pyr_c': 1}, {'vit_c': 1}, 0, 1000, 'gP'),
        'METS': ({'pyr_c': 1}, {'met_c': 1}, 0, 1000, 'gQ'),
        'T_vit': ({'vit_c': 1}, {'vit_e': 1}, -1000, 1000, None),
        'T_met': ({'met_c': 1}, {'met_e': 1}, -1000, 1000, None),
        'R_ATPM': ({'atp_c': 1}, {}, 1, 1000, None),
        'BIO': ({'vit_c': 1, 'met_c': 1, 'pyr_c': 1, 'atp_c': 1}, {}, 0, 1000, None),
    }
    if reverse:
        reactions['R_EX_vit_e'] = ({}, {'vit_e': 1}, -1000, 0, None)
    write_sbml(path, reactions, ['glc_e', 'vit_e', 'met_e', 'pyr_c', 'atp_c', 'vit_c', 'met_c'], genes=['gP', 'gQ'])
    model = B.N.ME.load_sbml_model(path)
    meta = {'path': str(path), 'sha256': model.sha256,
            'genes': [{'symbol': 'vitA', 'canonical': 'gP'}, {'symbol': 'metA', 'canonical': 'gQ'}]}
    return model, meta, B.context_for(model, meta)


def case():
    return {'case_id': 'test', 'species': 'Ec', 'kind': 'coculture', 'pair': 'V+M', 'environment': 'test medium',
        'primary_environment': True, 'assay': {'coculture_od_ratio': [1, 1]},
        'organisms': [{'strain_id': 'Ec:V', 'species': 'Ec', 'deleted_genes': ['vitA']},
                      {'strain_id': 'Ec:M', 'species': 'Ec', 'deleted_genes': ['metA']}],
        'medium': {'base': {k: 'declared' for k in B.N.BASE}, 'carbon_source': 'D-glucose',
                   'nitrogen_source': None, 'supplement': None}}


class NutrientCommunityTests(unittest.TestCase):
    def test_real_production_reciprocal_rescue_requires_unsupplied_uptake_and_separate_genotypes(self):
        for reverse in (False, True):
            with self.subTest(reverse=reverse), tempfile.TemporaryDirectory() as d:
                path = Path(d) / 'toy.xml'
                model, meta, context = fixture(path, reverse)
                plan = B.make_plan(case(), model, meta, context, 0., 'printed')
                self.assertEqual([m['blocked_reactions'] for m in plan['members']], [['VITS'], ['METS']])
                before = (model.lb.copy(), model.ub.copy())
                result = B.simulate(plan, context, model, path)
                self.assertEqual(result['predictions'], dict(monoculture_0=False, monoculture_1=False,
                    permissive=True, native_unsupplied=False, no_unsupplied_uptake=False))
                self.assertEqual(len(result['solver_calls']), 5)
                self.assertEqual([r['id'] for r in result['production']], B.PROTOCOL['solve_order'])
                self.assertEqual(result, B.simulate(plan, context, model, path))
                self.assertTrue((model.lb == before[0]).all() and (model.ub == before[1]).all())
                positive = result['production'][2]
                self.assertEqual(set(positive['per_organism_growth']), {'Ec:V', 'Ec:M'})
                self.assertTrue(positive['cross_fed'])
                self.assertEqual(result['promotion'], 'none')
                oracle = O.NutrientOracle(path, context, B.PROTOCOL)
                for member in plan['members']:
                    self.assertEqual(oracle.blocked(member['genes']), member['blocked_reactions'])
                for record in result['production']:
                    with self.subTest(kind=record['id']):
                        expected = oracle.simulate(plan, record['id'])
                        self.assertTrue(O.agrees(record, expected), (record, expected))
                        self.assertLessEqual(expected.get('max_balance_residual', 0), 1e-6)

    def test_outcomes_coverage_and_input_order_do_not_change_a_plan(self):
        with tempfile.TemporaryDirectory() as d:
            model, meta, context = fixture(Path(d) / 'toy.xml')
            c = case()
            expected = B.make_plan(c, model, meta, context, 0., 'printed')
            c.update(mean_od600=100, n_wells_measured=0, observed_growth=True)
            self.assertEqual(expected, B.make_plan(c, model, meta, context, 0., 'printed'))
            self.assertNotIn('observed_growth', json.dumps(expected))

    def test_invalid_pair_species_ratio_or_genotype_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            model, meta, context = fixture(Path(d) / 'toy.xml')
            for change in ('duplicate', 'species', 'ratio', 'gene'):
                c = case()
                if change == 'duplicate': c['organisms'][1]['strain_id'] = 'Ec:V'
                elif change == 'species': c['organisms'][1]['species'] = 'Bs'
                elif change == 'ratio': c['assay']['coculture_od_ratio'] = [2, 1]
                else: c['organisms'][0]['deleted_genes'] = ['missing']
                with self.subTest(change=change), self.assertRaises(ValueError):
                    B.make_plan(c, model, meta, context, 0., 'printed')

    def test_tampered_knockout_and_model_are_rejected_before_solving(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'toy.xml'
            model, meta, context = fixture(path)
            plan = B.make_plan(case(), model, meta, context, 0., 'printed')
            plan['members'][0]['blocked_reactions'] = []
            with self.assertRaisesRegex(ValueError, 'genotype'):
                B.simulate(plan, context, model, path)
            context['model_sha256'] = 'changed'
            with self.assertRaisesRegex(ValueError, 'model'):
                B.simulate(plan, context, model, path)

    def test_solver_unknowns_and_missing_endpoints_are_not_negative_labels(self):
        for status in ('inconclusive', 'skipped_budget'):
            self.assertIsNone(B.prediction({'status': status, 'outcome': 'error', 'mu': 0.}))
        self.assertFalse(B.prediction({'status': 'completed', 'outcome': 'infeasible_with_maintenance'}))
        self.assertIsNone(B.prediction({'status': 'completed', 'outcome': 'growth', 'mu': float('nan')}))
        self.assertIsNone(B.endpoint({}))
        self.assertIsNone(B.rescue(True, None, False))
        self.assertFalse(B.rescue(True, True, False))
        self.assertTrue(B.rescue(True, False, False))
        for invalid in (True, float('nan'), '0.1'):
            with self.assertRaises(ValueError): B.endpoint({'mean_od600': invalid})

    def test_freezing_and_execution_never_require_observations(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            model, meta, context = fixture(root / 'toy.xml')
            inputs = root / 'inputs.jsonl'
            inputs.write_text(json.dumps(case()) + '\n')
            registration = root / 'registration.json'
            B.N.write_json(registration, {'inputs_sha256': B.N.sha(inputs), 'models': {'Ec': meta}})
            frozen = B.prepare(inputs, registration, root / 'plan')
            self.assertEqual((frozen['cases'], frozen['plans']), (1, 6))
            receipt = B.run(root / 'plan', root / 'run', {'Ec': root / 'toy.xml'})
            self.assertEqual(sum(receipt['counts'].values()), 6)
            obs = root / 'obs.jsonl'
            observations = [dict(case_id='test', species='Ec', kind='coculture', environment='test medium', strains=['V', 'M'], mean_od600=.1)]
            for letter in ('V', 'M'):
                observations.append(dict(case_id=letter, species='Ec', kind='monoculture', environment='test medium', strains=[letter], mean_od600=.005))
            obs.write_text(''.join(json.dumps(r) + '\n' for r in observations))
            report = B.evaluate(root / 'plan', root / 'run', obs, root / 'evaluation.json')
            self.assertEqual(len(report['cases']), 6)
            self.assertTrue(all(r['observed_rescue_proxy'] for r in report['cases']))
            with (root / 'run/predictions.jsonl').open('a') as f: f.write('{}\n')
            with self.assertRaisesRegex(ValueError, 'provenance'):
                B.evaluate(root / 'plan', root / 'run', obs, root / 'bad.json')

    def test_shards_merge_to_exact_serial_results_and_missing_or_duplicate_shards_fail(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            model, meta, context = fixture(root / 'toy.xml')
            inputs = root / 'inputs.jsonl'
            inputs.write_text(json.dumps(case()) + '\n')
            registration = root / 'registration.json'
            B.N.write_json(registration, {'inputs_sha256': B.N.sha(inputs), 'models': {'Ec': meta}})
            B.prepare(inputs, registration, root / 'plan')
            B.run(root / 'plan', root / 'serial', {'Ec': root / 'toy.xml'})
            shards = []
            for i in range(2):
                folder = root / str(i)
                B.run(root / 'plan', folder, {'Ec': root / 'toy.xml'}, i, 2)
                shards.append(folder)
            B.merge(root / 'plan', shards[::-1], root / 'merged')
            self.assertEqual((root / 'serial/predictions.jsonl').read_bytes(), (root / 'merged/predictions.jsonl').read_bytes())
            report = O.run(root / 'plan', shards[0], root / 'oracle', {'Ec': root / 'toy.xml'})
            self.assertEqual((report['comparisons'], report['agree']), (15, 15))
            for incomplete in [shards[:1], [shards[0], shards[0]]]:
                with self.assertRaises(ValueError):
                    B.merge(root / 'plan', incomplete, root / 'bad_merge')
            with (shards[0] / 'predictions.jsonl').open('a') as f: f.write('{}\n')
            with self.assertRaisesRegex(ValueError, 'provenance'):
                B.merge(root / 'plan', shards, root / 'bad_merge')


if __name__ == '__main__':
    unittest.main()
