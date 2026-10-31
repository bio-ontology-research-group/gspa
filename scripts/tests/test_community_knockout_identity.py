"""Community interventions must preserve all deletions and reject unknown IDs."""
import copy
import unittest
from unittest.mock import patch

import libsbml
from test_community_metabolism import Base, CM
from metabolic_experiments import ManifestError


class TestCommunityKnockoutIdentity(Base):
    def test_two_deletion_records_do_not_overwrite_the_first(self):
        m = self.manifest(experiments=[{'id': 'combined', 'knockouts': [
            {'organism': 'P', 'reactions': ['VITS']},
            {'organism': 'P', 'reactions': ['DM_pyr']}]}])
        out = self.run_m(m)
        self.assertEqual(self.exp(out, 'community')['outcome'], 'common_growth')
        # Removing both the vitamin branch and overflow leaves no pyr disposal
        # compatible with ATP maintenance at zero coupled biomass growth.
        self.assertEqual(self.exp(out, 'user:combined')['outcome'], 'infeasible_with_maintenance')

    def test_label_deletions_of_both_isozymes_are_one_intervention(self):
        doc = libsbml.readSBMLFromFile(str(self.dir / 'P.xml'))
        model = doc.getModel(); fbc = model.getPlugin('fbc')
        fbc.getGeneProduct('gP').setLabel('producer_one')
        gp = fbc.createGeneProduct(); gp.setId('gP2'); gp.setLabel('producer_two')
        association = model.getReaction('VITS').getPlugin('fbc').getGeneProductAssociation()
        association.unsetAssociation()
        alternatives = association.createOr()
        for gene in ('gP', 'gP2'):
            alternatives.createGeneProductRef().setGeneProduct(gene)
        libsbml.writeSBMLToFile(doc, str(self.dir / 'P_isozymes.xml'))
        m = self.manifest(p_file='P_isozymes.xml', experiments=[
            {'id': 'one', 'knockouts': [{'organism': 'P', 'genes': ['producer_one']}]},
            {'id': 'both', 'knockouts': [{'organism': 'P', 'genes': ['producer_one']},
                                       {'organism': 'P', 'genes': ['producer_two']}]}])
        out = self.run_m(m)
        self.assertEqual(self.exp(out, 'user:one')['outcome'], 'common_growth')
        self.assertEqual(self.exp(out, 'user:both')['outcome'], 'no_common_growth')

    def test_invalid_knockouts_fail_before_any_solver_call(self):
        for knockout in ({'organism': 'P', 'genes': ['typo_gP']},
                         {'organism': 'P', 'reactions': ['MISSING_REACTION']},
                         {'organism': 'missing', 'genes': ['gP']},
                         {'organism': 'P', 'genes': 'gP'}):
            with self.subTest(knockout=knockout):
                m = self.manifest(experiments=[{'id': 'invalid', 'knockouts': [knockout]}])
                with patch.object(CM.BudgetedSolver, 'optimize') as solve:
                    with self.assertRaises(ManifestError):
                        self.run_m(m)
                solve.assert_not_called()

    def test_knockout_of_excluded_organism_is_not_silently_ignored(self):
        m = self.manifest(experiments=[{'id': 'excluded', 'organisms': ['A'],
            'knockouts': [{'organism': 'P', 'genes': ['gP']}]}])
        with self.assertRaisesRegex(ManifestError, 'excluded knockout organism'):
            self.run_m(m)


if __name__ == '__main__':
    unittest.main()
