"""Secretion capacity under perturbations, including nongrowing controls."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from test_metabolic_experiments import ME, REACTIONS, manifest, write_sbml
import masked_assignment_experiments as A


class TestExtracellularProducts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.path = Path(cls.tmp.name) / 'secreted.xml'
        write_sbml(cls.path)
        cls.model = ME.load_sbml_model(cls.path)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def standard(self, **cfg):
        m = manifest({'hypotheses': [{'id': 'H', 'gene': 'g6'}], 'evidence': [
            {'id': 'bgc', 'type': 'bgc', 'hypothesis_id': 'H', 'product_exchange': 'EX_sec_e',
             'condition_hint': 'supplement_removal', 'exchange': 'EX_nh4_e'}]},
            strategy='fixed', max_experiments=20, max_solver_calls=100,
            max_focus_reactions=0, max_baseline_focus_reactions=0,
            product_requirement='maintenance_only')
        m['config'].update(cfg)
        m['model'] = str(self.path)
        return m

    def masked(self, **cfg):
        m = self.standard(**cfg)
        m['assignment_recovery'] = {'mask_gene': 'g6',
            'hypotheses': [{'id': 'H0', 'assignments': {}},
                           {'id': 'H1', 'assignments': {'SEC': [[A.TARGET]]}},
                           {'id': 'H2', 'assignments': {'ACS': [[A.TARGET]]}}],
            'conditions': [{'id': 'C0', 'mods': [{'op': 'close', 'exchange': 'EX_nh4_e'}],
                            'product_exchange': 'EX_sec_e', 'evidence_channels': ['bgc']}]}
        return m

    def test_standard_measures_secretion_under_linked_knockout_and_nutrient_removal(self):
        out = ME.run_experiments(self.standard(), self.model)
        rows = {r['experiment_id']: r for r in out['experiments']}
        control = rows['product:EX_sec_e']['measurements']['product']
        no_nitrogen = rows['product:EX_sec_e@supplement_removal:EX_nh4_e']
        deleted = rows['product:EX_sec_e@knockout:H']['measurements']['product']
        self.assertGreater(control['max_secretion'], 0)
        self.assertEqual(no_nitrogen['measurements']['max_growth'], 0)
        self.assertGreater(no_nitrogen['measurements']['product']['max_secretion'], 0)
        self.assertEqual(deleted['max_secretion'], 0)
        self.assertNotIn('max_synthesis', control)
        self.assertEqual(out['promotion'], 'none')

    def test_unlinked_condition_is_not_combined_with_product(self):
        m = self.standard()
        m['context']['evidence'].append({'id': 'unrelated', 'type': 'transport',
                                        'substrate_exchange': 'EX_ac_e', 'role': 'carbon'})
        out = ME.run_experiments(m, self.model)
        ids = {r['experiment_id'] for r in out['experiments']}
        self.assertIn('product:EX_sec_e@knockout:H', ids)
        self.assertNotIn('product:EX_sec_e@carbon_source:EX_ac_e', ids)

    def test_masked_secretion_is_informative_without_growth_and_oracle_agrees(self):
        import sys
        sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'benchmark/v3'))
        import check_masked_assignment_cobra as oracle
        m = self.masked()
        out = ME.run_experiments(m, self.model)
        cells = {r['hypothesis_id']: r for r in out['cells']}
        self.assertEqual(cells['H1']['product_contrast'], 'secretion_capacity_lost')
        self.assertEqual(cells['H0']['product_contrast'], 'secretion_capacity_retained')
        self.assertEqual(cells['H2']['product_contrast'], 'secretion_capacity_retained')
        self.assertTrue(all(c['growth_status'] == 'inconclusive' and c['product_status'] == 'completed'
                            and c['status'] == 'completed' for c in cells.values()))
        self.assertEqual(cells['H1']['phenotypes']['wild_type']['product_target']['kind'], 'extracellular_secretion')
        report = oracle.compare(self.path, m, out)
        self.assertEqual((report['n_agree'], report['n_comparisons']), (6, 6))
        self.assertEqual((report['n_cell_agree'], report['n_cell_comparisons']), (3, 3))
        broken = copy.deepcopy(out)
        broken['cells'][0]['product_contrast'] = 'secretion_capacity_lost'
        self.assertLess(oracle.compare(self.path, m, broken)['n_cell_agree'], 3)

    def test_transport_loss_does_not_become_synthesis_evidence(self):
        model = copy.deepcopy(self.model)
        model.ub[model.rxn_index['T_sec']] = 0
        out = ME.run_experiments(self.masked(), model)
        self.assertTrue(all(c['product_contrast'] == 'no_secretion_capacity_in_control' for c in out['cells']))
        self.assertTrue(all(c['phenotypes']['wild_type']['product_interval']['maximum'] == 0 for c in out['cells']))

    def test_growth_requirement_and_solver_failure_remain_unknown(self):
        out = ME.run_experiments(self.masked(product_requirement='growth_fraction'), self.model)
        self.assertTrue(all(c['product_contrast'] == 'not_assessed_growth_requirement'
                            and c['product_status'] == 'inconclusive' for c in out['cells']))
        with patch.object(ME, 'linprog', side_effect=RuntimeError('injected LP failure')):
            out = ME.run_experiments(self.masked(), self.model)
        self.assertTrue(all(c['product_contrast'] == 'inconclusive' and c['status'] == 'inconclusive'
                            for c in out['cells']))
        self.assertTrue(all('product_interval' not in p for c in out['cells'] for p in c['phenotypes'].values()))

    def test_budget_and_agent_feedback_preserve_the_paired_secretion_result(self):
        short = ME.run_experiments(self.masked(max_solver_calls=21), self.model)
        self.assertEqual(short['cells'], [])
        self.assertEqual(short['solver_calls'], [])
        m = self.masked(strategy='qwen', max_experiments=2, max_solver_calls=44)
        m['assignment_recovery']['conditions'].append({'id': 'C1', 'mods': [], 'product_exchange': 'EX_sec_e'})
        prompts = []
        def chat(messages):
            prompts.append(json.loads(messages[1]['content']))
            return json.dumps({'experiment_id': 'C0' if len(prompts) == 1 else 'C1'})
        out = ME.run_experiments(m, self.model, chat=chat)
        self.assertEqual(len(prompts), 2)
        lost = next(c for c in prompts[1]['previous_results'] if c['hypothesis_id'] == 'H1')
        self.assertEqual(lost['product_contrast'], 'secretion_capacity_lost')
        self.assertEqual(lost['product_status'], 'completed')
        self.assertLessEqual(len(out['solver_calls']), 44)
        self.assertEqual(out['promotion'], 'none')

    def test_standard_growth_requirement_is_unassessed_not_zero_secretion(self):
        out = ME.run_experiments(self.standard(product_requirement='growth_fraction'), self.model)
        row = next(r for r in out['experiments'] if r['experiment_id'] ==
                   'product:EX_sec_e@supplement_removal:EX_nh4_e')
        product = row['measurements']['product']
        self.assertEqual(product['status'], 'not_assessed_no_growth')
        self.assertEqual(product['exchange'], 'EX_sec_e')
        self.assertEqual(product['target_kind'], 'extracellular_secretion')
        self.assertNotIn('max_secretion', product)

    def test_reversed_exchange_preserves_absolute_capacity_and_knockout_loss(self):
        reactions = copy.deepcopy(REACTIONS)
        reactions['EX_sec_e'] = ({}, {'sec_e': 1}, -1000, 1000, None)
        path = Path(self.tmp.name) / 'reversed.xml'
        write_sbml(path, reactions=reactions)
        model = ME.load_sbml_model(path)
        out = ME.run_experiments(self.standard(), model)
        product = next(r for r in out['experiments'] if r['experiment_id'] == 'product:EX_sec_e')
        # 10 glucose -> 20 pyruvate + 20 ATP. Respiration consumes 1/9
        # pyruvate to supply the ATP shortfall for secretion plus maintenance.
        self.assertAlmostEqual(product['measurements']['product']['max_secretion'], 179 / 9)
        out = ME.run_experiments(self.masked(), model)
        cell = next(c for c in out['cells'] if c['hypothesis_id'] == 'H1')
        self.assertAlmostEqual(cell['phenotypes']['wild_type']['product_interval']['maximum'], 179 / 9)
        self.assertEqual(cell['phenotypes']['target_knockout']['product_interval']['maximum'], 0)
        self.assertEqual(cell['product_contrast'], 'secretion_capacity_lost')

    def test_mixed_product_targets_keep_synthesis_and_secretion_separate(self):
        reactions = copy.deepcopy(REACTIONS)
        reactions['DM_sec'] = ({'sec_c': 1}, {}, 0, 1000, None)
        path = Path(self.tmp.name) / 'mixed.xml'
        write_sbml(path, reactions=reactions)
        model = ME.load_sbml_model(path)
        m = self.standard()
        m['baseline']['product_demands'] = ['DM_sec']
        m['context']['evidence'].append({'id': 'demand', 'type': 'bgc',
            'hypothesis_id': 'H', 'product_demand': 'DM_sec'})
        out = ME.run_experiments(m, model)
        menu = {r['id']: r for r in out['menu']}
        self.assertEqual(menu['product:EX_sec_e@knockout:H']['product_exchange'], 'EX_sec_e')
        self.assertNotIn('product_demand', menu['product:EX_sec_e@knockout:H'])
        self.assertEqual(menu['synthesis:DM_sec@knockout:H']['product_demand'], 'DM_sec')
        self.assertIsNone(menu['synthesis:DM_sec@knockout:H']['product_exchange'])
        m = self.masked(max_solver_calls=50)
        m['baseline']['product_demands'] = ['DM_sec']
        m['assignment_recovery']['conditions'].append({'id': 'C1', 'mods': [], 'product_demand': 'DM_sec'})
        out = ME.run_experiments(m, model)
        for cell in out['cells']:
            expected = 'secretion' if cell['condition_id'] == 'C0' else 'synthesis'
            self.assertTrue(cell['product_contrast'].startswith(expected + '_capacity_'))


if __name__ == '__main__':
    unittest.main()
