"""Intracellular product experiments through the production and masked loops."""
import copy
import json
from pathlib import Path
import tempfile
import unittest

from test_metabolic_experiments import ME, REACTIONS, manifest, write_sbml
import masked_assignment_experiments as A


class TestIntracellularProducts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.path = Path(cls.tmp.name) / 'intracellular.xml'
        rxns = copy.deepcopy(REACTIONS)
        del rxns['T_sec']  # synthesis is possible; extracellular secretion is not
        rxns['DM_sec'] = ({'sec_c': 1}, {}, 0, 1000, None)
        rxns['DM_rev'] = ({}, {'sec_c': 1}, -1000, 1000, None)
        rxns['DM_double'] = ({'sec_c': 2}, {}, 0, 1000, None)
        write_sbml(cls.path, reactions=rxns)
        cls.model = ME.load_sbml_model(cls.path)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def standard(self, **cfg):
        m = manifest({'hypotheses': [{'id': 'H', 'gene': 'g6'}], 'evidence': [
            {'id': 'bgc', 'type': 'bgc', 'hypothesis_id': 'H', 'product_demand': 'DM_sec',
             'condition_hint': 'supplement_removal', 'exchange': 'EX_nh4_e'}]},
            strategy='fixed', max_experiments=20, max_solver_calls=100,
            max_focus_reactions=0, max_baseline_focus_reactions=0, product_requirement='maintenance_only')
        m['config'].update(cfg)
        m['model'] = str(self.path)
        m['baseline']['product_demands'] = ['DM_sec', 'DM_rev']
        return m

    def masked(self, **cfg):
        m = self.standard(**cfg)
        m['assignment_recovery'] = {'mask_gene': 'g6',
            'hypotheses': [{'id': 'H0', 'assignments': {}},
                           {'id': 'H1', 'assignments': {'SEC': [[A.TARGET]]}},
                           {'id': 'H2', 'assignments': {'ACS': [[A.TARGET]]}}],
            'conditions': [{'id': 'C0', 'mods': [{'op': 'close', 'exchange': 'EX_nh4_e'}],
                            'product_demand': 'DM_sec', 'evidence_channels': ['bgc']}]}
        return m

    def test_synthesis_without_export_and_informative_condition_and_knockout(self):
        out = ME.run_experiments(self.standard(), self.model)
        rows = {r['experiment_id']: r for r in out['experiments']}
        base = rows['synthesis:DM_sec']['measurements']['product']
        closed = rows['synthesis:DM_sec@supplement_removal:EX_nh4_e']
        deleted = rows['synthesis:DM_sec@knockout:H']['measurements']['product']
        self.assertGreater(base['max_synthesis'], 0)
        self.assertEqual(base['min_synthesis_at_required_state'], 0)
        self.assertEqual(base['target_kind'], 'intracellular_synthesis_capacity')
        self.assertNotIn('max_secretion', base)
        self.assertEqual(closed['measurements']['max_growth'], 0)
        self.assertGreater(closed['measurements']['product']['max_synthesis'], 0)
        self.assertEqual(deleted['max_synthesis'], 0)
        self.assertEqual(out['promotion'], 'none')
        external = self.standard()
        external['context']['evidence'] = [{'id': 'ex', 'type': 'bgc', 'product_exchange': 'EX_sec_e'}]
        ex_out = ME.run_experiments(external, self.model)
        product = next(r for r in ex_out['experiments'] if r['experiment_id'] == 'product:EX_sec_e')
        self.assertEqual(product['measurements']['product']['max_secretion'], 0)

    def test_reversed_demand_uses_the_drain_direction_without_opening_uptake(self):
        m = self.standard()
        m['context']['evidence'][0]['product_demand'] = 'DM_rev'
        normal = ME.run_experiments(self.standard(), self.model)
        reverse = ME.run_experiments(m, self.model)
        n = next(r for r in normal['experiments'] if r['experiment_id'] == 'synthesis:DM_sec')
        r = next(r for r in reverse['experiments'] if r['experiment_id'] == 'synthesis:DM_rev')
        self.assertAlmostEqual(n['measurements']['product']['max_synthesis'], r['measurements']['product']['max_synthesis'])
        lb, ub = ME.baseline_bounds(self.model, m, ME.resolve_config(m))
        self.assertEqual(lb[self.model.rxn_index['DM_sec']], 0)
        self.assertEqual(ub[self.model.rxn_index['DM_rev']], 0)

    def test_invalid_targets_fail_before_any_masked_solver_calls(self):
        from unittest.mock import patch
        changes = [{'product_demand': 'EX_sec_e'}, {'product_demand': 'SEC'},
                   {'product_demand': 'ATPM'}, {'product_demand': 'DM_double'},
                   {'product_demand': 'MISSING'}, {'product_demand': 'DM_sec', 'product_exchange': 'EX_sec_e'}]
        for target in changes:
            with self.subTest(target=target):
                m = self.masked()
                m['baseline']['product_demands'] += ['EX_sec_e', 'SEC', 'ATPM', 'DM_double', 'MISSING']
                m['assignment_recovery']['conditions'][0].update(target)
                with patch.object(ME.BudgetedSolver, 'optimize') as solve:
                    with self.assertRaises(ME.ManifestError):
                        ME.run_experiments(m, self.model)
                solve.assert_not_called()

    def test_undeclared_or_deleted_or_gene_associated_reporter_is_not_admissible(self):
        cases = []
        m = self.masked(); m['baseline']['product_demands'] = []; cases.append((m, self.model))
        m = self.masked(); m['assignment_recovery']['conditions'][0]['mods'].append(
            {'op': 'knockout', 'reactions': ['DM_sec'], 'genes': []}); cases.append((m, self.model))
        model = copy.deepcopy(self.model); model.gpr['DM_sec'] = [['g2']]; cases.append((self.masked(), model))
        for m, model in cases:
            with self.subTest(case=len(m['baseline']['product_demands'])):
                with self.assertRaises(ME.ManifestError):
                    ME.run_experiments(m, model)

    def test_growth_requirement_is_unassessed_not_zero_product(self):
        out = ME.run_experiments(self.standard(product_requirement='growth_fraction'), self.model)
        row = next(r for r in out['experiments'] if r['experiment_id'] == 'synthesis:DM_sec@supplement_removal:EX_nh4_e')
        self.assertEqual(row['measurements']['product']['status'], 'not_assessed_no_growth')
        self.assertNotIn('max_synthesis', row['measurements']['product'])

    def test_standard_demand_reserves_two_product_calls(self):
        for budget, status, used in [(3, 'skipped_budget', 1), (4, 'completed', 4)]:
            out = ME.run_experiments(self.standard(max_solver_calls=budget, max_experiments=1,
                panel=['synthesis:DM_sec']), self.model)
            row = next(r for r in out['experiments'] if r['experiment_id'] == 'synthesis:DM_sec')
            self.assertEqual(row['status'], status)
            self.assertEqual(out['solver']['used'], used)

    def test_masked_demand_distinguishes_assignments_without_growth(self):
        out = ME.run_experiments(self.masked(), self.model)
        cells = {r['hypothesis_id']: r for r in out['cells']}
        self.assertEqual(cells['H1']['product_contrast'], 'synthesis_capacity_lost')
        self.assertEqual(cells['H0']['product_contrast'], 'synthesis_capacity_retained')
        self.assertEqual(cells['H2']['product_contrast'], 'synthesis_capacity_retained')
        self.assertTrue(all(r['status'] == 'completed' and r['contrast'] == 'no_growing_control' for r in cells.values()))
        self.assertEqual(cells['H1']['phenotypes']['wild_type']['product_target']['kind'], 'intracellular_synthesis_capacity')
        self.assertEqual(out['promotion'], 'none')
        self.assertIsNone(out['calibrated_probability'])

    def test_masked_budget_preserves_full_comparison_rectangle(self):
        out = ME.run_experiments(self.masked(max_solver_calls=21), self.model)
        self.assertEqual(out['cells'], [])
        self.assertEqual(out['solver_calls'], [])
        out = ME.run_experiments(self.masked(max_solver_calls=22), self.model)
        self.assertEqual(len(out['cells']), 3)
        self.assertEqual(len(out['solver_calls']), 18)

    def test_failed_masked_synthesis_is_unknown_and_never_a_zero_capacity(self):
        from unittest.mock import patch
        with patch.object(ME, 'linprog', side_effect=RuntimeError('injected LP failure')):
            out = ME.run_experiments(self.masked(), self.model)
        self.assertEqual(len(out['solver_calls']), 6)
        for cell in out['cells']:
            self.assertEqual(cell['status'], 'inconclusive')
            self.assertEqual(cell['product_contrast'], 'inconclusive')
            for phenotype in cell['phenotypes'].values():
                self.assertNotIn('product_interval', phenotype)
                self.assertEqual(phenotype['product_target']['kind'], 'intracellular_synthesis_capacity')

    def test_masked_agent_knows_the_product_requirement_before_selecting(self):
        prompts = []
        def chat(messages):
            prompts.append(json.loads(messages[1]['content']))
            return json.dumps({'experiment_id': 'C0'})
        out = ME.run_experiments(self.masked(strategy='qwen'), self.model, chat=chat)
        self.assertEqual(prompts[0]['product_requirement'], 'maintenance_only')
        self.assertEqual(prompts[0]['conditions'][0]['product_demand'], 'DM_sec')
        self.assertEqual(out['selection'][0]['status'], 'agent_selected')

    def test_informative_growth_loss_survives_an_unassessed_product_test(self):
        m = self.masked(product_requirement='growth_fraction')
        m['assignment_recovery']['conditions'][0]['mods'] = []
        m['assignment_recovery']['hypotheses'][1]['assignments'] = {'BIOMASS': [[A.TARGET]]}
        out = ME.run_experiments(m, self.model)
        cell = next(c for c in out['cells'] if c['hypothesis_id'] == 'H1')
        self.assertEqual(cell['contrast'], 'knockout_no_growth')
        self.assertEqual(cell['growth_status'], 'completed')
        self.assertEqual(cell['product_status'], 'inconclusive')
        self.assertEqual(cell['product_contrast'], 'not_assessed_growth_requirement')
        self.assertEqual(cell['status'], 'completed')

    def test_blocked_and_forced_reporter_bounds_are_visible_without_rewriting_them(self):
        for lo, hi in [(0, 0), (1, 1000)]:
            model = copy.deepcopy(self.model)
            col = model.rxn_index['DM_sec']; model.lb[col], model.ub[col] = lo, hi
            out = ME.run_experiments(self.standard(panel=['synthesis:DM_sec']), model)
            row = next(r for r in out['experiments'] if r['experiment_id'] == 'synthesis:DM_sec')
            p = row['measurements']['product']
            self.assertEqual(p['reporter_drain_bounds'], [lo, hi])
            self.assertEqual(p['min_synthesis_at_required_state'], lo)
            if hi == 0:
                self.assertEqual(p['max_synthesis'], 0)
                self.assertTrue(p['capacity_at_reporter_bound'])

    def test_growth_can_force_byproduct_synthesis_with_zero_reporter_lower_bound(self):
        reactions = copy.deepcopy(REACTIONS)
        del reactions['T_coup']
        reactions['DM_coup'] = ({'coup_c': 1}, {}, 0, 1000, None)
        path = Path(self.tmp.name) / 'coupled.xml'; write_sbml(path, reactions=reactions)
        model = ME.load_sbml_model(path)
        m = self.standard(product_requirement='growth_fraction', panel=['synthesis:DM_coup'])
        m['baseline']['product_demands'] = ['DM_coup']
        m['context']['evidence'] = [{'id': 'coupled', 'type': 'bgc', 'product_demand': 'DM_coup'}]
        out = ME.run_experiments(m, model)
        p = next(r for r in out['experiments'] if r['experiment_id'] == 'synthesis:DM_coup')['measurements']['product']
        self.assertEqual(p['reporter_drain_bounds'][0], 0)
        self.assertGreater(p['min_synthesis_at_required_state'], 0)
        self.assertTrue(p['synthesis_unavoidable_in_model'])

    def test_oracle_detects_missing_intervals_and_wrong_endpoint_classification(self):
        import sys
        sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'benchmark/v3'))
        import check_masked_assignment_cobra as oracle
        m = self.masked(); out = ME.run_experiments(m, self.model)
        correct = oracle.compare(self.path, m, out)
        self.assertEqual(correct['n_agree'], correct['n_comparisons'])
        self.assertEqual(correct['n_cell_agree'], correct['n_cell_comparisons'])
        broken = copy.deepcopy(out)
        del broken['cells'][0]['phenotypes']['wild_type']['product_interval']
        report = oracle.compare(self.path, m, broken)
        self.assertLess(report['n_agree'], report['n_comparisons'])
        broken = copy.deepcopy(out)
        broken['cells'][0]['product_contrast'] = 'synthesis_capacity_lost'
        report = oracle.compare(self.path, m, broken)
        self.assertLess(report['n_cell_agree'], report['n_cell_comparisons'])

    def test_qwen_sees_typed_synthesis_feedback(self):
        prompts = []
        def chat(messages):
            prompts.append(json.loads(messages[1]['content']))
            return json.dumps({'experiment_id': 'synthesis:DM_sec' if len(prompts) == 1 else 'synthesis:DM_sec@knockout:H'})
        out = ME.run_experiments(self.standard(strategy='qwen', max_experiments=2), self.model, chat=chat)
        self.assertEqual(len(prompts), 2)
        self.assertTrue(all(x.get('valid') for x in out['selection_log']))
        result = next(r for r in prompts[1]['completed_results'] if r['experiment_id'] == 'synthesis:DM_sec')
        self.assertGreater(result['measurements']['product']['max_synthesis'], 0)
        self.assertNotIn('max_secretion', result['measurements']['product'])
        self.assertEqual(prompts[1]['schema'], ME.QWEN_PROMPT_SCHEMA)


if __name__ == '__main__':
    unittest.main()
