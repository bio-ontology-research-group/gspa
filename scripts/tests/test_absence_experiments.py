"""Real LP mechanism checks; these synthetic cases are not biological validation."""
import copy
import json
import subprocess
import sys
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import absence_experiments as A
import metabolic_experiments as ME
from test_metabolic_experiments import Base, manifest, write_sbml, REACTIONS


def problem(**config):
    m = manifest({}, strategy='fixed', max_experiments=3, max_solver_calls=24)
    m['config'].update(config)
    m['baseline']['medium']['EX_vit_e'] = 0
    m['absence_comparison'] = {
        'reaction': 'VITS', 'candidate_function_key': 'locus_g5|EC|EC:2.3.3.3',
        'assembly': {'status': 'COMPLETE', 'completeness_fraction': 1.0, 'source': 'synthetic closed genome'},
        'function_search': {'status': 'COMPLETE', 'source': 'synthetic function search'},
        'external_supply': {'exchange': 'EX_vit_e', 'uptake_limit': 2.0,
                            'transport_reactions': ['T_vit'], 'source': 'explicit synthetic medium'},
        'conditions': [{'id': 'control', 'mods': []},
                       {'id': 'withdrawal', 'mods': [{'op': 'close', 'exchange': 'EX_vit_e'}]},
                       {'id': 'transport_ko', 'mods': [{'op': 'knockout', 'reactions': ['T_vit']}]}]}
    return m


class AbsenceTests(Base):
    def run_problem(self, request=None, model=None):
        return ME.run_experiments(request or problem(), model or self.model)

    def cells(self, result):
        return {(c['condition_id'], c['hypothesis']): c for c in result['cells']}

    def test_annotation_sequence_loss_and_supply_have_distinct_intervention_patterns(self):
        out = self.run_problem()
        cells = self.cells(out)
        wt = lambda c, h: cells[c, h]['phenotypes']['wild_type']['class']
        ko = lambda c, h: cells[c, h]['phenotypes']['candidate_assignment_knockout']['class']
        for condition in ('control', 'withdrawal', 'transport_ko'):
            self.assertEqual(wt(condition, 'missing_annotation'), 'growth')
            self.assertEqual(ko(condition, 'missing_annotation'), 'no_growth')
            self.assertEqual(wt(condition, 'missing_sequence'), 'growth')
            self.assertEqual(ko(condition, 'missing_sequence'), 'growth')
            self.assertEqual(wt(condition, 'pathway_loss'), 'no_growth')
        self.assertEqual(wt('control', 'external_supply'), 'growth')
        self.assertEqual(ko('control', 'external_supply'), 'growth')
        self.assertEqual(wt('withdrawal', 'external_supply'), 'no_growth')
        self.assertEqual(wt('transport_ko', 'external_supply'), 'no_growth')
        self.assertEqual(len(out['solver_calls']), 24)
        self.assertTrue(all(c['status'] == 'completed' for c in out['cells']))
        self.assertEqual(out['gene_attribution'], 'unverified')
        self.assertIsNone(out['selected_explanation'])
        self.assertIsNone(out['calibrated_probability'])
        self.assertEqual(out['promotion'], 'none')

    def test_supply_has_exact_finite_cap_and_can_only_rescue_with_transport(self):
        m = problem()
        m['absence_comparison']['external_supply']['uptake_limit'] = .03
        cell = self.cells(self.run_problem(m))['control', 'external_supply']
        self.assertAlmostEqual(cell['phenotypes']['wild_type']['growth'], .03)
        closed = copy.deepcopy(self.model)
        closed.ub[closed.rxn_index['T_vit']] = 0
        cell = self.cells(self.run_problem(m, closed))['control', 'external_supply']
        self.assertEqual(cell['phenotypes']['wild_type']['class'], 'no_growth')

    def test_common_supplied_medium_keeps_loss_and_external_supply_indistinguishable(self):
        m = problem()
        m['baseline']['medium']['EX_vit_e'] = 2
        out = self.run_problem(m)
        pair = next(c for c in out['comparisons'] if c['hypotheses'] == ['pathway_loss', 'external_supply'])
        self.assertEqual(pair['status'], 'same_evaluated_growth_classes')
        self.assertEqual(pair['distinguishing_interventions'], [])

    def test_complete_draft_and_unknown_metadata_do_not_fabricate_biological_absence(self):
        results = []
        for status, fraction in [('COMPLETE', 1), ('DRAFT', .5), ('UNKNOWN', None)]:
            m = problem()
            m['absence_comparison']['assembly'].update(status=status, completeness_fraction=fraction)
            m['absence_comparison']['function_search']['status'] = 'UNKNOWN'
            out = self.run_problem(m)
            self.assertEqual(out['problem']['absence_comparison']['assembly']['status'], status)
            results.append(out['cells'])
        self.assertEqual(results[0], results[1])
        self.assertEqual(results[0], results[2])

    def test_removed_candidate_remains_unavailable_while_sequence_activity_is_a_hypothesis(self):
        m = problem()
        m['absence_comparison']['candidate_function_key'] = None
        m['absence_comparison']['assembly']['status'] = 'DRAFT'
        out = self.run_problem(m)
        self.assertEqual(len(out['solver_calls']), 9)
        for cell in out['cells']:
            if cell['hypothesis'] == 'missing_annotation':
                self.assertEqual(cell['status'], 'unavailable')
                self.assertEqual(cell['reason'], 'no_retained_candidate')
            else:
                self.assertEqual(cell['phenotypes']['candidate_assignment_knockout']['status'], 'unavailable')

    def test_missing_supply_is_unavailable_not_negative(self):
        m = problem()
        m['absence_comparison']['external_supply'] = None
        out = self.run_problem(m)
        self.assertEqual(len(out['solver_calls']), 18)
        self.assertTrue(all(c['status'] == 'unavailable' for c in out['cells'] if c['hypothesis'] == 'external_supply'))
        self.assertTrue(all(c['status'] == 'unavailable' for c in out['comparisons'] if 'external_supply' in c['hypotheses']))

    def test_budget_reserves_whole_condition_and_never_privileges_one_hypothesis(self):
        for budget, expected in [(0, 0), (7, 0), (8, 4), (15, 4), (16, 8)]:
            out = self.run_problem(problem(max_solver_calls=budget))
            self.assertEqual(len(out['cells']), expected)
            self.assertLessEqual(len(out['solver_calls']), budget)
            self.assertEqual(len(out['visited_conditions']) + len(out['unvisited_conditions']), 3)

    def test_original_gpr_and_hidden_context_do_not_influence_comparisons(self):
        m = problem()
        alternate = copy.deepcopy(self.model)
        alternate.gpr = {'VITS': [['g3', 'g4']], 'RESP': [['g5']]}
        first = self.run_problem(m)
        m['context'] = {'truth': 'pathway_loss', 'expected': 'external_supply'}
        second = self.run_problem(m, alternate)
        self.assertEqual(first['cells'], second['cells'])
        self.assertEqual(first['problem_sha256'], second['problem_sha256'])
        self.assertNotIn('truth', json.dumps(second['problem']))

    def test_positive_requirement_cannot_be_relaxed_by_focal_or_maintenance_knockout(self):
        model = copy.deepcopy(self.model)
        model.lb[model.rxn_index['VITS']] = 1
        out = self.run_problem(model=model)
        cell = self.cells(out)['control', 'pathway_loss']
        self.assertEqual(cell['phenotypes']['wild_type']['class'], 'infeasible_under_constraints')
        m = problem()
        m['absence_comparison']['conditions'] = [{'id': 'maintenance_ko', 'mods': [{'op': 'knockout', 'reactions': ['ATPM']}]}]
        out = self.run_problem(m)
        self.assertTrue(all(c['phenotypes']['wild_type']['class'] == 'infeasible_under_constraints' for c in out['cells']))

    def test_solver_errors_remain_unknown_not_negative_or_equal(self):
        with patch.object(ME.BudgetedSolver, 'optimize', return_value=ME.LPResult('error', None, None, None, 'failure')):
            out = self.run_problem()
        self.assertTrue(all(c['status'] == 'inconclusive' for c in out['cells']))
        self.assertTrue(all(c['status'] == 'unavailable' for c in out['comparisons']))

    def test_unrelated_internal_or_undeclared_supply_ports_rejected(self):
        for ports in ([], ['VITS'], ['T_glc'], ['T_vit', 'T_vit'], ['UNKNOWN']):
            m = problem()
            m['absence_comparison']['external_supply']['transport_reactions'] = ports
            with self.assertRaises(ME.ManifestError):
                self.run_problem(m)
        m = problem()
        m['absence_comparison']['external_supply']['exchange'] = 'ATPM'
        with self.assertRaises(ME.ManifestError):
            self.run_problem(m)

    def test_undeclared_alternative_transport_is_closed_and_declared_alternative_is_retained(self):
        reactions = dict(REACTIONS)
        reactions['T_vit_alternative'] = reactions['T_vit']
        path = Path(self.tmp.name) / 'alternative-port.xml'
        write_sbml(path, reactions=reactions)
        model = ME.load_sbml_model(path)
        m = problem()
        out = self.run_problem(m, model)
        self.assertEqual(out['closed_undeclared_transport_reactions'], ['T_vit_alternative'])
        self.assertEqual(self.cells(out)['transport_ko', 'external_supply']['phenotypes']['wild_type']['class'], 'no_growth')
        m['absence_comparison']['external_supply']['transport_reactions'].append('T_vit_alternative')
        out = self.run_problem(m, model)
        self.assertEqual(out['closed_undeclared_transport_reactions'], [])
        self.assertEqual(self.cells(out)['transport_ko', 'external_supply']['phenotypes']['wild_type']['class'], 'growth')

    def test_opposite_exchange_and_transport_signs_preserve_supply_results(self):
        reactions = dict(REACTIONS)
        for rid in ('EX_vit_e', 'T_vit'):
            reactants, products, lo, hi, gpr = reactions[rid]
            reactions[rid] = (products, reactants, -hi, -lo, gpr)
        path = Path(self.tmp.name) / 'reverse-port.xml'
        write_sbml(path, reactions=reactions)
        out = self.run_problem(model=ME.load_sbml_model(path))
        self.assertEqual(out['cells'], self.run_problem()['cells'])

    def test_lp_infeasibility_and_zero_growth_do_not_count_as_distinct_observed_phenotypes(self):
        out = self.run_problem()
        for cell in out['cells']:
            if cell['hypothesis'] == 'external_supply':
                cell['phenotypes'] = {p: {'class': 'infeasible_under_constraints'} for p in A.PERTURBATIONS}
        pair = next(p for p in A.compare(out['cells']) if p['hypotheses'] == ['pathway_loss', 'external_supply'])
        self.assertEqual(pair['status'], 'same_evaluated_growth_classes')

    def test_invalid_rates_metadata_objectives_and_modes_fail_even_at_zero_budget(self):
        mutations = [lambda p: p['external_supply'].update(uptake_limit=float('inf')),
                     lambda p: p['external_supply'].update(uptake_limit=True),
                     lambda p: p['assembly'].update(completeness_fraction=2),
                     lambda p: p['assembly'].update(source=''),
                     lambda p: p.update(reaction='BIOMASS'),
                     lambda p: p.update(reaction='ATPM'),
                     lambda p: p['conditions'].append(p['conditions'][0]),
                     lambda p: p['conditions'][0].update(mods=[{'op': 'knockout', 'genes': ['g5'], 'reactions': []}]),
                     lambda p: p.update(truth='missing_annotation')]
        for mutation in mutations:
            m = problem(max_solver_calls=0)
            mutation(m['absence_comparison'])
            with self.assertRaises(ME.ManifestError):
                self.run_problem(m)
        for config in ({'strategy': 'evidence'}, {'max_solver_calls': True}, {'max_experiments': -1}):
            with self.assertRaises(ME.ManifestError):
                self.run_problem(problem(**config))

    def test_actual_python_cli_emits_comparison(self):
        m = problem()
        m['model'] = str(self.path)
        path = Path(self.tmp.name) / 'absence.json'
        out = path.with_name('absence-result.json')
        path.write_text(json.dumps(m))
        proc = subprocess.run([sys.executable, str(Path(ME.__file__)), '--manifest', str(path), '--out', str(out)], capture_output=True, text=True)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        result = json.loads(out.read_text())
        self.assertEqual(result['schema'], A.SCHEMA)
        self.assertEqual(len(result['cells']), 12)
