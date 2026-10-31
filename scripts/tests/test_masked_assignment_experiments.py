"""Actual LP regressions for masked function hypotheses; synthetic mechanism data."""
import copy
import json
import io
import subprocess
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import masked_assignment_experiments as A
import metabolic_experiments as ME
from test_metabolic_experiments import Base, manifest


def problem(**config):
    m = manifest({}, max_solver_calls=60, max_experiments=2, strategy='fixed')
    m['config'].update(config)
    m['assignment_recovery'] = {
        'mask_gene': 'locus_g3',
        'hypotheses': [
            {'id': 'H0', 'assignments': {}, 'evidence': {}},
            {'id': 'H1', 'assignments': {'ACS': [[A.TARGET]]}, 'evidence': {'sequence': .3}},
            {'id': 'H2', 'assignments': {'VITS': [[A.TARGET]]}, 'evidence': {'domain': .5}}],
        'conditions': [
            {'id': 'C0', 'mods': [], 'evidence_channels': ['sequence']},
            {'id': 'C1', 'mods': [{'op': 'carbon_swap', 'exchange': 'EX_ac_e'}], 'evidence_channels': ['domain']}]}
    return m


class MaskingTests(Base):
    def test_entire_target_and_clause_removed_but_or_isozyme_survives(self):
        model = copy.deepcopy(self.model)
        model.gpr['ACS'] = [['g3', 'g4'], ['g2']]
        before = copy.deepcopy(model.gpr)
        masked = A.mask_model(model, 'locus_g3')
        self.assertEqual(masked.gpr['ACS'], [['g2']])
        self.assertNotIn('g3', masked.gene_alias.values())
        self.assertNotIn('locus_g3', masked.gene_alias)
        self.assertEqual(model.gpr, before)
        self.assertEqual(masked.rxns if hasattr(masked, 'rxns') else masked.rxn_ids, model.rxn_ids)
        self.assertTrue((masked.S != model.S).nnz == 0)

    def test_public_problem_identical_when_hidden_correct_assignment_changes(self):
        alternate = copy.deepcopy(self.model)
        alternate.gpr['ACS'] = []
        alternate.gpr['RESP'] = alternate.gpr.get('RESP', []) + [['g3']]
        alternate.sha256 = 'different private source'
        m = problem()
        m['context'] = {'expected': 'ACS', 'true_protein_id': 'locus_g3'}
        a, ma = A.public_problem(m, self.model)
        m['context']['expected'] = 'RESP'
        b, mb = A.public_problem(m, alternate)
        self.assertEqual(a, b)
        self.assertEqual(ma.sha256, mb.sha256)
        self.assertNotIn('locus_g3', json.dumps(a))
        self.assertNotIn('expected', json.dumps(a))
        first = ME.run_experiments(m, self.model)
        second = ME.run_experiments(m, alternate)
        self.assertEqual(first['cells'], second['cells'])

    def test_independent_partners_preserved_in_candidate_and_clauses(self):
        public, masked = A.public_problem(problem(), self.model)
        h = {'assignments': {'ACS': [[A.TARGET, 'g4']]}}
        model = A.assigned_model(masked, h)
        self.assertIn('ACS', model.blocked_by_genes({'g4'}))
        self.assertIn('ACS', model.blocked_by_genes({A.TARGET}))
        self.assertEqual(masked.gpr['ACS'], [])

    def test_missing_null_private_labels_aliases_and_unmapped_reactions_fail(self):
        mutations = [lambda p: p['hypotheses'].pop(0),
                     lambda p: p['hypotheses'][1].update(expected=True),
                     lambda p: p['hypotheses'][1].update(id='H_locus_g3'),
                     lambda p: p['hypotheses'][1].update(assignments={'NO_RXN': [[A.TARGET]]})]
        for mutate in mutations:
            m = problem()
            mutate(m['assignment_recovery'])
            with self.assertRaises(ME.ManifestError):
                A.public_problem(m, self.model)


class SimulationTests(Base):
    def run_problem(self, m=None, chat=None):
        return ME.run_experiments(m or problem(), self.model, chat=chat)

    def test_condition_and_target_knockouts_jointly_disable_or_alternatives(self):
        self.model = copy.deepcopy(self.model)
        self.model.gpr['ACS'] = [['g3'], ['g2']]
        m = problem()
        m['assignment_recovery']['hypotheses'][1]['assignments']['ACS'] = [[A.TARGET, 'g4']]
        m['assignment_recovery']['conditions'] = [{'id': 'C0', 'mods': [
            {'op': 'carbon_swap', 'exchange': 'EX_ac_e'},
            {'op': 'knockout', 'reactions': [], 'genes': ['g2']}], 'evidence_channels': ['domain']}]
        result = self.run_problem(m)
        cell = next(c for c in result['cells'] if c['hypothesis_id'] == 'H1')
        self.assertIn('ACS', cell['phenotypes']['target_knockout']['knocked_out_reactions'])
        self.assertGreater(cell['phenotypes']['wild_type']['growth'], 0)
        self.assertIsNone(cell['phenotypes']['target_knockout']['growth'])

    def test_multiple_knockout_modifications_are_one_combined_gene_intervention(self):
        self.model = copy.deepcopy(self.model)
        self.model.gpr['ACS'] = [['g3'], ['g2']]
        m = problem()
        lb, ub = self.model.lb.copy(), self.model.ub.copy()
        killed = ME.apply_mods(self.model, m, ME.DEFAULT_CONFIG, lb, ub, [
            {'op': 'knockout', 'reactions': [], 'genes': ['g2']},
            {'op': 'knockout', 'reactions': [], 'genes': ['g3']}])
        self.assertIn('ACS', killed)
        self.assertEqual(ub[self.model.rxn_index['ACS']], 0)

    def test_carbon_specific_assignment_separates_null_and_wrong_function(self):
        result = self.run_problem()
        cells = {(c['condition_id'], c['hypothesis_id']): c for c in result['cells']}
        self.assertGreater(cells['C0', 'H1']['phenotypes']['target_knockout']['growth'], 0)
        self.assertIsNone(cells['C1', 'H1']['phenotypes']['target_knockout']['growth'])
        self.assertEqual(cells['C1', 'H1']['contrast'], 'knockout_infeasible_with_maintenance')
        self.assertGreater(cells['C1', 'H0']['phenotypes']['target_knockout']['growth'], 0)
        self.assertGreater(cells['C1', 'H2']['phenotypes']['target_knockout']['growth'], 0)
        self.assertEqual(result['promotion'], 'none')
        self.assertIsNone(result['calibrated_probability'])
        self.assertLessEqual(len(result['solver_calls']), result['max_solver_calls'])

    def test_alternative_flux_intervals_are_computed_at_explicit_growth_floor(self):
        result = self.run_problem()
        cell = next(c for c in result['cells'] if c['condition_id'] == 'C1' and c['hypothesis_id'] == 'H1')
        fva = cell['phenotypes']['wild_type']['flux_intervals']['ACS']
        self.assertEqual(fva['statuses'], ['optimal', 'optimal'])
        self.assertGreater(fva['minimum'], 0)
        self.assertGreaterEqual(fva['maximum'], fva['minimum'])

    def test_complete_candidate_comparison_reserved_before_spending_budget(self):
        m = problem(max_solver_calls=9)
        result = self.run_problem(m)
        self.assertEqual(result['visited_conditions'], [])
        self.assertEqual(result['solver_calls'], [])
        result = self.run_problem(problem(max_solver_calls=10))
        self.assertEqual(len(result['visited_conditions']), 1)
        self.assertEqual(len(result['cells']), 3)

    def test_candidate_order_does_not_change_predictions_or_budget_allocation(self):
        m = problem(max_solver_calls=10)
        first = self.run_problem(m)
        m['assignment_recovery']['hypotheses'].reverse()
        second = self.run_problem(m)
        self.assertEqual(first['cells'], second['cells'])
        self.assertEqual(first['solver_calls'], second['solver_calls'])

    def test_same_budget_enforced_for_fixed_random_evidence_and_agent(self):
        for strategy in ('fixed', 'random', 'evidence', 'qwen'):
            result = self.run_problem(problem(strategy=strategy, max_solver_calls=10), chat=lambda _: '{"experiment_id":"C1"}')
            self.assertLessEqual(len(result['solver_calls']), 10)
            self.assertEqual(len(result['cells']), 3)

    def test_agent_sees_only_masked_candidates_and_results(self):
        messages = []
        def chat(prompt):
            messages.append(prompt)
            p = json.loads(prompt[-1]['content'])
            return json.dumps({'experiment_id': p['conditions'][0]['id']})
        m = problem(strategy='qwen')
        m['context'] = {'expected_function': 'private answer', 'gene': 'locus_g3'}
        result = self.run_problem(m, chat=chat)
        self.assertEqual(len(messages), 2)
        self.assertNotIn('locus_g3', json.dumps(messages))
        self.assertNotIn('private answer', json.dumps(messages))
        self.assertTrue(json.loads(messages[1][-1]['content'])['previous_results'])
        self.assertTrue(all(s['status'] == 'agent_selected' for s in result['selection']))

    def test_agent_failure_never_credited_as_fallback_success(self):
        result = self.run_problem(problem(strategy='qwen'), chat=lambda _: '{"experiment_id":"not-a-condition"}')
        self.assertEqual(result['selection'][0]['status'], 'agent_failed')
        self.assertEqual(result['solver_calls'], [])

    def test_live_adapter_accounts_for_usage_without_changing_response_contract(self):
        payload = {'choices': [{'message': {'content': '{"experiment_id":"C0"}'}}],
                   'usage': {'prompt_tokens': 17, 'completion_tokens': 6}}
        with patch.object(ME.urllib.request, 'urlopen', return_value=io.BytesIO(json.dumps(payload).encode())):
            chat = ME.make_qwen_chat(ME.resolve_config(problem()))
            self.assertEqual(chat([]), '{"experiment_id":"C0"}')
        self.assertEqual(chat.usage['prompt_tokens'], 17)
        self.assertEqual(chat.usage['completion_tokens'], 6)
        with patch.object(ME.urllib.request, 'urlopen', side_effect=OSError('injected')):
            with self.assertRaises(OSError):
                chat([])
        self.assertEqual(chat.usage['http_attempts'], 2)
        self.assertEqual(chat.usage['failed_calls'], 1)

    def test_product_flux_range_is_conditional_not_expression(self):
        m = problem(max_experiments=1)
        m['assignment_recovery']['conditions'][0]['product_exchange'] = 'EX_sec_e'
        result = self.run_problem(m)
        product = result['cells'][0]['phenotypes']['wild_type']['product_interval']
        self.assertEqual(product['minimum'], 0)
        self.assertGreater(product['maximum'], 0)
        self.assertEqual(result['promotion'], 'none')

    def test_product_maintenance_only_mode_can_test_a_nongrowing_feasible_state(self):
        model = copy.deepcopy(self.model)
        model.ub[model.rxn_index['BIOMASS']] = 0
        m = problem(max_experiments=1, product_requirement='maintenance_only')
        m['assignment_recovery']['conditions'][0]['product_exchange'] = 'EX_sec_e'
        result = ME.run_experiments(m, model)
        wt = result['cells'][0]['phenotypes']['wild_type']
        self.assertEqual(wt['growth'], 0)
        self.assertGreater(wt['product_interval']['maximum'], 0)
        m['config']['product_requirement'] = 'growth_fraction'
        growth_required = ME.run_experiments(m, model)
        self.assertNotIn('product_interval', growth_required['cells'][0]['phenotypes']['wild_type'])

    def test_actual_cli_accepts_assignment_manifest(self):
        m = problem(max_solver_calls=10)
        m['model'] = str(self.path)
        input_path = Path(self.tmp.name) / 'assignment.json'
        output_path = Path(self.tmp.name) / 'assignment-result.json'
        input_path.write_text(json.dumps(m))
        proc = subprocess.run([sys.executable, str(Path(ME.__file__)), '--manifest', str(input_path), '--out', str(output_path)], capture_output=True, text=True)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(json.loads(output_path.read_text())['schema'], A.SCHEMA)

    def test_solver_exceptions_consume_budget_and_remain_inconclusive(self):
        with patch.object(ME, 'linprog', side_effect=RuntimeError('injected solver failure')):
            result = self.run_problem(problem(max_solver_calls=10))
        self.assertEqual(len(result['solver_calls']), 6)
        self.assertTrue(all(c['status'] == 'inconclusive' for c in result['cells']))
        self.assertTrue(all(c['status'] == 'error' for c in result['solver_calls']))
        self.assertTrue(all(c['phenotypes']['target_knockout']['growth'] is None for c in result['cells']))

    def test_invalid_budget_and_growth_floor_rejected_at_zero_budget(self):
        for overrides in ({'max_solver_calls': 1.5}, {'max_solver_calls': True}, {'max_solver_calls': 0, 'growth_fraction': 0}):
            with self.assertRaises(ME.ManifestError):
                self.run_problem(problem(**overrides))

    def test_independent_cobra_parser_and_glpk_agree_with_growth_and_fva(self):
        sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'benchmark/v3'))
        import check_masked_assignment_cobra as oracle
        m = problem()
        m['assignment_recovery']['conditions'][0]['product_exchange'] = 'EX_sec_e'
        result = self.run_problem(m)
        report = oracle.compare(self.path, m, result)
        self.assertEqual(report['n_comparisons'], 12)
        self.assertEqual(report['n_agree'], 12, report)


if __name__ == '__main__':
    unittest.main()
