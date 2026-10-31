"""Independent LP oracle must not carry perturbations between experiment cells."""
import io
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
try:
    import cobra
    import libsbml
    import check_masked_assignment_cobra as O
except ImportError:
    cobra = None


@unittest.skipIf(cobra is None, 'COBRApy/libSBML required')
class OracleIsolation(unittest.TestCase):
    def setUp(self):
        model = cobra.Model('oracle_isolation')
        a, b = (cobra.Metabolite(x, compartment='c') for x in ('a', 'b'))
        uptake = cobra.Reaction('EX_a', lower_bound=-10, upper_bound=1000)
        uptake.add_metabolites({a: -1})
        synth = cobra.Reaction('SYNTH')
        synth.add_metabolites({a: -1, b: 1})
        synth.gene_reaction_rule = 'isozyme'
        biomass = cobra.Reaction('BIO')
        biomass.add_metabolites({b: -1})
        model.add_reactions([uptake, synth, biomass])
        model.objective = biomass
        stream = io.StringIO()
        cobra.io.write_sbml_model(model, stream)
        document = libsbml.readSBMLFromString(stream.getvalue())
        gene = document.getModel().getPlugin('fbc').createGeneProduct()
        gene.setId('G_partner')
        gene.setLabel('partner')
        # Actual SBML can declare genes unused by the baseline GPR. A later
        # independently supported complex hypothesis can introduce their role.
        self.model = cobra.io.read_sbml_model(io.StringIO(libsbml.writeSBMLToString(document)))
        self.model.solver = 'glpk'
        self.public = {'baseline': {'medium': {'EX_a': 10}},
                       'objectives': {'biomass': 'BIO', 'maintenance': 'none'}}
        self.cfg = {'uptake_rate': 10, 'growth_threshold': 1e-6,
                    'growth_fraction': .95, 'product_requirement': 'growth'}
        self.complex = {'assignments': {'SYNTH': [['MASKED_TARGET', 'partner']]}}
        self.null = {'assignments': {}}

    def simulate(self, hypothesis, genes=(), knockout=False):
        condition = {'mods': [{'op': 'knockout', 'genes': list(genes), 'reactions': []}]}
        return O.simulate(self.model, self.public, hypothesis, condition, self.cfg, knockout)

    def test_unused_gene_knockout_does_not_mutate_base_model(self):
        self.simulate(self.null, ['partner'])
        self.assertTrue(self.model.genes.get_by_id('partner').functional)
        self.assertEqual(self.model.reactions.SYNTH.gene_reaction_rule, 'isozyme')
        self.assertNotIn('MASKED_TARGET', self.model.genes)

    def test_complex_growth_and_fva_are_invariant_to_prior_null_partner_knockout(self):
        isolated = O.simulate(self.model.copy(), self.public, self.complex,
            {'mods': [{'op': 'knockout', 'genes': ['isozyme'], 'reactions': []}]}, self.cfg, False)
        self.assertEqual(isolated['growth'], 10)
        self.assertEqual(isolated['flux_intervals']['SYNTH']['minimum'], 9.5)
        self.simulate(self.null, ['partner'])
        sequential = self.simulate(self.complex, ['isozyme'])
        self.assertEqual(sequential, isolated)

    def test_condition_and_target_knockouts_disable_both_isozyme_alternatives(self):
        self.assertEqual(self.simulate(self.complex, ['isozyme'])['growth'], 10)
        self.assertEqual(self.simulate(self.complex, ['isozyme'], knockout=True)['growth'], 0)
        self.assertEqual(self.simulate(self.complex, ['partner', 'isozyme'])['growth'], 0)
        self.assertEqual(self.simulate(self.complex, ['isozyme'])['growth'], 10)


if __name__ == '__main__':
    unittest.main()
