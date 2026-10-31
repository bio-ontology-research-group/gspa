"""Known genotype, assay matching, bounded media, independent oracle and unknowns."""
import copy
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / 'scripts/tests'))
import run_nutrient_phenotypes as N
import check_nutrient_cobra as O
from test_metabolic_experiments import write_sbml, manifest, REACTIONS, libsbml


def case():
    return {'case_id': 'test', 'species': 'Ec', 'kind': 'monoculture', 'environment': 'NoN + NoC',
            'primary_environment': False, 'organisms': [{'strain_id': 'Ec:M', 'deleted_genes': ['metA']}],
            'medium': {'base': {k: 'printed concentration' for k in N.BASE},
                       'carbon_source': None, 'nitrogen_source': None, 'supplement': None}}


class NutrientPhenotypeTests(unittest.TestCase):
    def test_no_added_carbon_or_nitrogen_retains_basal_organics(self):
        model = SimpleNamespace(rxn_index={'R_EX_trp__L_e': 0, 'R_EX_thm_e': 1})
        medium, bindings = N.medium_for(case(), model, 0., 'printed')
        self.assertEqual(medium, {'R_EX_trp__L_e': .1, 'R_EX_thm_e': .01})
        self.assertTrue(any(b['status'] == 'unrepresented' for b in bindings))

    def test_missing_thiamine_is_not_replaced_by_thiamine_phosphate(self):
        medium, bindings = N.medium_for(case(), SimpleNamespace(rxn_index={'R_EX_tmp_e': 1}), 20., 'printed')
        self.assertEqual(medium, {})
        self.assertIn({'component': 'thiamine', 'exchange': 'R_EX_thm_e', 'status': 'unrepresented',
                       'assumed_rate': .01}, bindings)

    def test_assay_supplement_and_unprinted_sulfate_are_separate(self):
        c = case()
        c['medium']['supplement'] = 'serine'
        model = SimpleNamespace(rxn_index={'R_EX_ser__L_e': 0, 'R_EX_so4_e': 1, 'R_EX_nh4_e': 2})
        medium, _ = N.medium_for(c, model, 2., 'printed')
        self.assertEqual(medium, {'R_EX_ser__L_e': 5.})
        c['medium']['supplement'] = None
        self.assertEqual(N.medium_for(c, model, 2., 'printed')[0], {})
        self.assertEqual(N.medium_for(c, model, 2., 'additional_inorganic')[0], {'R_EX_so4_e': 1000.})
        c['medium']['nitrogen_source'] = 'ammonium sulfate'
        self.assertEqual(N.medium_for(c, model, 2., 'printed')[0], {'R_EX_nh4_e': 10., 'R_EX_so4_e': 10.})

    def test_unregistered_medium_and_arms_rejected(self):
        c = case()
        with self.assertRaisesRegex(ValueError, 'arm'):
            N.medium_for(c, SimpleNamespace(rxn_index={}), 9., 'printed')
        c['medium']['base']['extra'] = 'unknown'
        with self.assertRaisesRegex(ValueError, 'ingredient'):
            N.medium_for(c, SimpleNamespace(rxn_index={}), 0., 'printed')

    def test_label_changes_do_not_change_plan_and_pairs_are_rejected(self):
        model = SimpleNamespace(gene_alias={'g5': 'g5'}, objective_terms=[('bio', 1.)], objective_sense='maximize',
                                rxn_index={'R_ATPM': 0, 'R_EX_h_e': 1}, lb=[1., -1000.])
        meta = {'genes': [{'symbol': 'metA', 'canonical': 'g5'}], 'path': 'model', 'sha256': 'hash'}
        from unittest.mock import patch
        with patch.object(N.ME, 'validate_boundary'):
            c = case()
            first = N.make_plan(c, model, meta, 0., 'printed')
            c['mean_od600'], c['above_author_monoculture_threshold'] = 100., True
            self.assertEqual(first, N.make_plan(c, model, meta, 0., 'printed'))
            c['organisms'].append(c['organisms'][0])
            with self.assertRaisesRegex(ValueError, 'single-organism'):
                N.make_plan(c, model, meta, 0., 'printed')

    def test_failed_solver_is_unknown_not_negative(self):
        for status in ('inconclusive', 'skipped_budget'):
            self.assertIsNone(N.prediction({'status': status, 'outcome': 'solver_error', 'measurements': {'max_growth': 0.}}))
        self.assertFalse(N.prediction({'status': 'completed', 'outcome': 'infeasible_with_maintenance'}))
        self.assertIsNone(N.prediction({'status': 'completed', 'outcome': 'growth', 'measurements': {'max_growth': float('nan')}}))
        self.assertEqual(N.confusion([(None, False), (True, None)])['unknown'], 2)
        self.assertFalse(O.agreement({'status': 'inconclusive', 'outcome': 'solver_error'},
                                     {'status': 'optimal', 'growth': 0.}))

    def test_oracle_result_is_independent_of_condition_and_knockout_history(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'model.xml'
            write_sbml(path)
            conditions = []
            for supplied, oxygen, genes in [(False, 10., []), (False, 10., ['g5']),
                                            (True, 10., ['g5']), (False, 0., []),
                                            (True, 0., ['g1', 'g2'])]:
                m = manifest()
                m['baseline']['medium']['EX_o2_e'] = oxygen
                if not supplied:
                    m['baseline']['medium'].pop('EX_vit_e')
                conditions.append((m, genes))
            expected = []
            for m, genes in conditions:
                fresh = O.cobra.io.read_sbml_model(path, f_replace={})
                fresh.solver = 'glpk'
                expected.append(O.simulate(fresh, m, genes))
            reused = O.cobra.io.read_sbml_model(path, f_replace={})
            reused.solver = 'glpk'
            for index in (4, 0, 3, 1, 2, 0, 4, 2, 1, 3):
                result = O.simulate(reused, *conditions[index])
                self.assertEqual(result['status'], expected[index]['status'])
                self.assertEqual(result['blocked_reactions'], expected[index]['blocked_reactions'])
                if result['status'] == 'optimal':
                    self.assertAlmostEqual(result['growth'], expected[index]['growth'], places=9)
                    self.assertLessEqual(result['max_balance_residual'], 1e-6)
                    self.assertLessEqual(result['max_bound_violation'], 1e-6)

    def test_actual_production_pair_rescue_reversed_exchange_and_independent_gpr(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'model.xml'
            write_sbml(path)
            model = N.ME.load_sbml_model(path)
            oracle = O.cobra.io.read_sbml_model(path, f_replace={})
            oracle.solver = 'glpk'
            before = {r.id: r.bounds for r in oracle.reactions}
            for supplied in (False, True):
                m = manifest({'hypotheses': [{'id': 'observed_genotype', 'genes': ['g5']}]}, strategy='fixed',
                             max_experiments=1, max_solver_calls=2, max_focus_reactions=0, max_baseline_focus_reactions=0)
                m['model_sha256_expected'] = model.sha256
                if not supplied:
                    m['baseline']['medium'].pop('EX_vit_e')
                plan = {'case_id': 'test', 'arm': 'test', 'genes': ['g5'], 'manifest': m}
                result = N.simulate(plan, model)
                self.assertEqual(result['prediction'], supplied)
                for rec in result['production']['experiments']:
                    genes = [] if rec['kind'] == 'baseline' else ['g5']
                    independent = O.simulate(oracle, m, genes)
                    self.assertTrue(O.agreement(rec, independent))
                    self.assertEqual(independent['blocked_reactions'], rec['measurements']['knocked_out_reactions'])
            self.assertEqual(before, {r.id: r.bounds for r in oracle.reactions})

    def test_no_effect_isozyme_deletion_stays_no_effect_and_wt_needs_no_fake_gene(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'model.xml'
            write_sbml(path)
            doc = libsbml.readSBML(str(path))
            gpr = doc.getModel().getReaction('VITS').getPlugin('fbc').getGeneProductAssociation()
            gpr.unsetAssociation()
            either = gpr.createOr()
            for gene in ('g5', 'g6'):
                either.createGeneProductRef().setGeneProduct(gene)
            libsbml.writeSBMLToFile(doc, str(path))
            model = N.ME.load_sbml_model(path)
            m = manifest({'hypotheses': [{'id': 'observed_genotype', 'genes': ['g5']}]}, strategy='fixed',
                         max_experiments=1, max_solver_calls=2, max_focus_reactions=0, max_baseline_focus_reactions=0)
            m['model_sha256_expected'] = model.sha256
            m['baseline']['medium'].pop('EX_vit_e')
            plan = {'case_id': 'test', 'arm': 'test', 'genes': ['g5'], 'manifest': m}
            result = N.simulate(plan, model)
            self.assertTrue(result['prediction'])
            self.assertEqual(result['production']['experiments'][1]['measurements']['knocked_out_reactions'], [])
            m['context']['hypotheses'] = []
            plan['genes'] = []
            wt = N.simulate(plan, model)
            self.assertTrue(wt['prediction'])
            self.assertEqual(len(wt['production']['experiments']), 1)
            self.assertIn('Runner.run', wt['entrypoint'])
            plan['manifest']['model_sha256_expected'] = 'tampered'
            with self.assertRaisesRegex(ValueError, 'hash mismatch'):
                N.simulate(plan, model)

    def evaluation_fixture(self, root):
        plan_root, run_root = root / 'plan', root / 'run'
        plan_root.mkdir()
        run_root.mkdir()
        plans, predictions, observations = [], [], []
        for species in ('Ec', 'Bs'):
            for kind, value in [('monoculture', False), ('supplemented_monoculture', True)]:
                cid = species + ':' + kind
                plan = {'case_id': cid, 'arm': 'printed:oxygen=20', 'species': species, 'kind': kind,
                        'environment': 'Glut + Gluc', 'strain_id': species + ':M', 'primary_environment': True}
                plans.append(plan)
                predictions.append({'case_id': cid, 'arm': plan['arm'], 'plan_sha256': N.digest(plan),
                                    'prediction': value if species == 'Ec' else None})
                # Deliberately omit Bs supplementation: missing assay stays unknown.
                if species == 'Ec' or not value:
                    observations.append({'case_id': cid, 'mean_od600': .1 if value else .005})
        (plan_root / 'plans.jsonl').write_text(''.join(json.dumps(r) + '\n' for r in plans))
        (run_root / 'predictions.jsonl').write_text(''.join(json.dumps(r) + '\n' for r in predictions))
        N.write_json(plan_root / 'protocol.json', {'protocol': N.PROTOCOL, 'implementation': N.code_hashes(),
                                                'plans_sha256': N.sha(plan_root / 'plans.jsonl')})
        N.write_json(run_root / 'run.json', {'protocol_sha256': N.sha(plan_root / 'protocol.json'),
                                          'plans_sha256': N.sha(plan_root / 'plans.jsonl'),
                                          'predictions_sha256': N.sha(run_root / 'predictions.jsonl')})
        obs = root / 'observations.jsonl'
        obs.write_text(''.join(json.dumps(r) + '\n' for r in observations))
        return plan_root, run_root, obs

    def test_scoring_preserves_missing_assays_and_pairs_only_matched_species_strains(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            plan, run, obs = self.evaluation_fixture(root)
            scored = N.evaluate(plan, run, obs, root / 'evaluation.json')
            results = {r['species']: r for r in scored['paired_rescue']}
            self.assertEqual(results['Ec']['tp'], 1)
            self.assertEqual(results['Bs']['unknown'], 1)
            self.assertIsNone(results['Bs']['accuracy'])
            self.assertEqual(len(scored['cases']), 4)

    def test_evaluation_rejects_tampered_receipt_and_incomplete_prediction_catalogue(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            plan, run, obs = self.evaluation_fixture(root)
            predictions = run / 'predictions.jsonl'
            predictions.write_text(predictions.read_text().splitlines()[0] + '\n')
            with self.assertRaisesRegex(ValueError, 'provenance'):
                N.evaluate(plan, run, obs, root / 'bad.json')
            receipt = json.loads((run / 'run.json').read_text())
            receipt['predictions_sha256'] = N.sha(predictions)
            N.write_json(run / 'run.json', receipt)
            with self.assertRaisesRegex(ValueError, 'catalogue'):
                N.evaluate(plan, run, obs, root / 'bad.json')

    def test_frozen_execution_rejects_source_or_plan_change(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            plan, _, _ = self.evaluation_fixture(root)
            N.verify_plan(plan)
            frozen = json.loads((plan / 'protocol.json').read_text())
            frozen['implementation'] = {}
            N.write_json(plan / 'protocol.json', frozen)
            with self.assertRaisesRegex(ValueError, 'source changed'):
                N.verify_plan(plan)


if __name__ == '__main__':
    unittest.main()
