"""Unit checks for the absence-explanation scoring and the false-assignment injection."""
import json
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / 'benchmark/v3'))
sys.path.insert(0, str(REPO / 'benchmark/v3/absence_crossfeeding'))
sys.path.insert(0, str(REPO / 'scripts'))
sys.path.insert(0, str(REPO / 'scripts/tests'))
try:
    import absence_truth as AT
    import false_assignments as FA
    from test_metabolic_experiments import write_sbml, ME
except ImportError:  # pragma: no cover - optional scientific stack
    AT = None


@unittest.skipIf(AT is None, 'SciPy and libSBML required')
class InjectionTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        path = Path(self.directory.name) / 'model.xml'
        write_sbml(path)
        self.model = ME.load_sbml_model(path)

    def tearDown(self):
        self.directory.cleanup()

    def test_rate_zero_changes_nothing_and_draws_are_reproducible(self):
        mapping, receipt = FA.inject(self.model, 0.0, 1)
        self.assertEqual(receipt['reactions_perturbed'], 0)
        self.assertEqual(FA.apply_gpr(self.model, mapping).gpr, self.model.gpr)
        again, _ = FA.inject(self.model, 0.5, 7)
        once, _ = FA.inject(self.model, 0.5, 7)
        self.assertEqual(again, once)
        self.assertNotEqual(once, FA.inject(self.model, 0.5, 8)[0])

    def test_perturbed_rule_is_a_single_other_gene_and_moves_the_knockout_set(self):
        mapping, receipt = FA.inject(self.model, 1.0, 3)
        for change in receipt['changes']:
            if 'new_gene' in change:
                self.assertEqual(mapping[change['reaction']], [[change['new_gene']]])
                own = {g for clause in change['old_rule'] for g in clause}
                self.assertNotIn(change['new_gene'], own)
        perturbed = FA.apply_gpr(self.model, mapping)
        genes = sorted({g for dnf in self.model.gpr.values() for clause in dnf for g in clause})
        moved = [g for g in genes
                 if self.model.blocked_by_genes({g}) != perturbed.blocked_by_genes({g})]
        self.assertTrue(moved, 'reassigning every rule must change some knockout set')

    def test_targeted_mode_only_touches_focal_rules_and_spares_their_genes(self):
        genes = sorted({g for dnf in self.model.gpr.values() for clause in dnf for g in clause})
        focal = genes[:1]
        mapping, receipt = FA.inject(self.model, 1.0, 1, 'targeted', focal)
        self.assertEqual(receipt['mode'], 'targeted')
        self.assertTrue(receipt['eligible_reactions'] >= 1)
        touched = {c['reaction'] for c in receipt['changes']}
        for rid, dnf in self.model.gpr.items():
            mentions = any(g in focal for clause in dnf for g in clause)
            self.assertEqual(rid in touched, mentions)
        for change in receipt['changes']:
            self.assertNotIn(change['new_gene'], focal)
        with self.assertRaises(ValueError):
            FA.inject(self.model, 1.0, 1, 'targeted', [])


@unittest.skipIf(AT is None, 'SciPy and libSBML required')
class ScoringTests(unittest.TestCase):
    def test_growth_class_mapping_and_measurement_threshold(self):
        self.assertEqual(AT.comparable('infeasible_under_constraints'), 'no_growth')
        self.assertIsNone(AT.comparable('unknown'))
        self.assertEqual(AT.observed_class({'mean_od600': 0.5}), 'growth')
        self.assertEqual(AT.observed_class({'mean_od600': 0.0}), 'no_growth')
        self.assertIsNone(AT.observed_class(None))
        self.assertIsNone(AT.observed_class({'mean_od600': None}))
        with self.assertRaises(ValueError):
            AT.observed_class({'mean_od600': float('nan')})

    def test_wilson_interval_brackets_the_point_estimate(self):
        low, high = AT.wilson(8, 10)
        self.assertLess(low, 0.8)
        self.assertGreater(high, 0.8)
        self.assertEqual(AT.wilson(0, 0), [None, None])

    def test_supplement_must_be_a_single_exchange(self):
        unsupplemented = {'manifest': {'baseline': {'medium': {'R_EX_a_e': 10.0}}}}
        supplemented = {'manifest': {'baseline': {'medium': {'R_EX_a_e': 10.0, 'R_EX_b_e': 1.0}}}}
        self.assertEqual(AT.medium_difference(unsupplemented, supplemented), ('R_EX_b_e', 1.0))
        two = {'manifest': {'baseline': {'medium': {'R_EX_a_e': 10.0, 'R_EX_b_e': 1.0, 'R_EX_c_e': 1.0}}}}
        with self.assertRaises(ValueError):
            AT.medium_difference(unsupplemented, two)


if __name__ == '__main__':
    unittest.main()
