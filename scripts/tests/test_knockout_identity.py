"""Explicit interventions must resolve identity before changing any bounds."""
from dataclasses import replace
import copy
import unittest
import libsbml

import numpy as np
from test_metabolic_experiments import Base, manifest, ME


class TestKnockoutIdentity(Base):
    def modify(self, mods, model=None):
        model = model or self.model
        lb, ub = model.lb.copy(), model.ub.copy()
        killed = ME.apply_mods(model, manifest(), ME.DEFAULT_CONFIG, lb, ub, mods)
        return killed, lb, ub

    def test_label_and_id_knockouts_have_identical_effect(self):
        direct = self.modify([{"op": "knockout", "genes": ["g3"], "reactions": []}])
        label = self.modify([{"op": "knockout", "genes": ["locus_g3"], "reactions": []}])
        self.assertEqual(label[0], ["ACS"])
        self.assertEqual(direct[0], label[0])
        np.testing.assert_array_equal(direct[1], label[1])
        np.testing.assert_array_equal(direct[2], label[2])

    def test_mixed_aliases_are_resolved_before_combining_isozyme_deletions(self):
        gpr = copy.deepcopy(self.model.gpr)
        gpr["GLYC1"] = [["g1"], ["g2"]]
        model = replace(self.model, gpr=gpr)
        one = self.modify([{"op": "knockout", "genes": ["locus_g1"], "reactions": []}], model)
        self.assertNotIn("GLYC1", one[0])
        both = self.modify([{"op": "knockout", "genes": ["locus_g1"], "reactions": []},
                            {"op": "knockout", "genes": ["g2"], "reactions": []}], model)
        self.assertIn("GLYC1", both[0])

    def test_complex_subunit_and_explicit_reaction_are_both_removed(self):
        model = copy.deepcopy(self.model)
        model.gpr['GLYC1'] = [['g1', 'g2']]
        killed, _, _ = self.modify([{'op': 'knockout', 'genes': ['locus_g2', 'g2'],
                                    'reactions': ['ACS']}], model)
        self.assertTrue({'GLYC1', 'GLYC2', 'ACS'} <= set(killed))

    def test_ambiguous_labels_and_id_collisions_are_rejected_at_load(self):
        for label in ('locus_g1', 'g1'):
            with self.subTest(label=label):
                doc = libsbml.readSBMLFromFile(str(self.path))
                doc.getModel().getPlugin('fbc').getGeneProduct('g2').setLabel(label)
                path = self.path.parent / ('collision-' + label + '.xml')
                libsbml.writeSBMLToFile(doc, str(path))
                with self.assertRaisesRegex(ME.ManifestError, 'ambiguous SBML gene label'):
                    ME.load_sbml_model(path)

    def test_unknown_gene_fails_before_any_intervention_changes_bounds(self):
        lb, ub = self.model.lb.copy(), self.model.ub.copy()
        with self.assertRaisesRegex(ME.ManifestError, "unknown knockout gene"):
            ME.apply_mods(self.model, manifest(), ME.DEFAULT_CONFIG, lb, ub,
                [{"op": "close", "exchange": "EX_glc_e"},
                 {"op": "knockout", "genes": ["typo_g3"], "reactions": []}])
        np.testing.assert_array_equal(lb, self.model.lb)
        np.testing.assert_array_equal(ub, self.model.ub)

    def test_unknown_reaction_is_a_manifest_error_before_mutation(self):
        lb, ub = self.model.lb.copy(), self.model.ub.copy()
        with self.assertRaisesRegex(ME.ManifestError, "unknown knockout reaction"):
            ME.apply_mods(self.model, manifest(), ME.DEFAULT_CONFIG, lb, ub,
                [{"op": "close", "exchange": "EX_glc_e"},
                 {"op": "knockout", "genes": [], "reactions": ["NOT_IN_MODEL"]}])
        np.testing.assert_array_equal(lb, self.model.lb)
        np.testing.assert_array_equal(ub, self.model.ub)

    def test_malformed_intervention_identifiers_are_rejected(self):
        for field, value in (("genes", "g3"), ("genes", [None]), ("genes", [""]),
                             ("reactions", "ACS"), ("reactions", [False])):
            with self.subTest(field=field, value=value):
                mod = {"op": "knockout", "genes": [], "reactions": [], field: value}
                with self.assertRaises(ME.ManifestError):
                    self.modify([mod])

    def test_unmapped_hypothesis_annotation_is_not_an_executed_gene_deletion(self):
        m = manifest({"hypotheses": [{"id": "H", "gene": "unknown_candidate", "reactions": ["ACS"]}],
                      "evidence": []})
        menu, rejected = ME.build_menu(self.model, m)
        self.assertFalse(rejected)
        self.assertEqual(menu[0]["unmapped_candidate_genes"], ["unknown_candidate"])
        self.assertEqual(self.modify(menu[0]["mods"])[0], ["ACS"])

    def test_masked_assignment_loop_combines_label_condition_with_target_deletion(self):
        from test_masked_assignment_experiments import problem
        model = copy.deepcopy(self.model)
        model.gpr["ACS"] = [["g3"], ["g2"]]
        m = problem()
        m["assignment_recovery"]["conditions"] = [{"id": "C0", "mods": [
            {"op": "carbon_swap", "exchange": "EX_ac_e"},
            {"op": "knockout", "genes": ["locus_g2"], "reactions": []}], "evidence_channels": []}]
        out = ME.run_experiments(m, model)
        cell = next(c for c in out["cells"] if c["hypothesis_id"] == "H1")
        self.assertGreater(cell["phenotypes"]["wild_type"]["growth"], 0)
        self.assertIn("ACS", cell["phenotypes"]["target_knockout"]["knocked_out_reactions"])
        self.assertEqual(cell["phenotypes"]["target_knockout"]["status"], "infeasible")

    def test_masked_assignment_rejects_unknown_condition_gene_before_solver(self):
        from test_masked_assignment_experiments import problem
        from unittest.mock import patch
        m = problem()
        m["assignment_recovery"]["conditions"][0]["mods"] = [
            {"op": "knockout", "genes": ["unknown_condition_gene"], "reactions": []}]
        with patch.object(ME.BudgetedSolver, "optimize") as solve:
            with self.assertRaisesRegex(ME.ManifestError, "unknown knockout gene"):
                ME.run_experiments(m, self.model)
        solve.assert_not_called()


if __name__ == "__main__":
    unittest.main()
