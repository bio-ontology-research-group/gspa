"""Biological identity controls on the full retained primary tables and real SBML."""
from pathlib import Path
import copy,json,sys,unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import audit_germicidin_identity as G
ROOT=Path(__file__).resolve().parents[3]

class GermicidinIdentityTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report=G.audit(ROOT/'benchmark/v3/results/germicidin-identity-v1/sources/cihak2017.html',
            ROOT/'benchmark/v3/results/secondary-source-intake-v1/sources/repository/ModelFiles/xml/Sco-GEM.xml')
    def test_complete_primary_tables_are_retained(self):
        self.assertEqual({k:len(rows)-1 for k,rows in self.report['primary_tables'].items()},
                         {'1':15,'2':15,'3':3})
        self.assertEqual(len(self.report['model_products']),4)
        self.assertEqual(len(self.report['comparisons']),3)
        self.assertTrue(all(row['accepted_identity'] is None for row in self.report['comparisons']))
    def test_same_germicidin_a_name_does_not_establish_same_compound(self):
        rows=[r for r in self.report['comparisons'] if r['observation']['Secondary metabolite']=='Germicidin A']
        self.assertEqual(len(rows),2)
        for r in rows:
            self.assertEqual(r['observation']['Molecular formula'],'C11H16O3')
            self.assertEqual(round(float(r['observation']['Acquired mass (m/z)'])),197)
            self.assertEqual(len(r['name_matches']),1)
            self.assertEqual(r['disposition'],'name_formula_conflict')
            self.assertEqual(r['name_matches'][0]['model_formula'],'C10H14O3')
            self.assertEqual(r['name_matches'][0]['model_nominal_protonated_mass'],183)
            self.assertFalse(r['name_matches'][0]['formula_agrees'])
            self.assertEqual(r['formula_compatible_candidates'],['M_germicidinC_c'])
            self.assertIsNone(r['accepted_identity'])
    def test_same_formula_does_not_disambiguate_germicidin_b(self):
        r=self.report['comparisons'][1]
        self.assertEqual(r['observation']['Secondary metabolite'],'Germicidin B')
        self.assertTrue(r['name_matches'][0]['formula_agrees'])
        self.assertEqual(r['formula_compatible_candidates'],['M_germicidinA_c','M_germicidinB_c'])
        a,b=self.report['model_products'][:2]
        self.assertEqual(a['composition'],b['composition'])
        self.assertEqual({v['metabolite_id'] for v in a['reactants']},{'M_3oiC6ACP_c','M_emcoa__S_c','M_h_c'})
        self.assertEqual({v['metabolite_id'] for v in b['reactants']},{'M_3oaiC7ACP_c','M_mmcoa__S_c','M_h_c'})
        self.assertEqual(r['disposition'],'formula_ambiguous')
        self.assertIsNone(r['accepted_identity'])

    def test_source_phase_and_mass_disagreement_remain_explicit(self):
        a,b,c=self.report['comparisons']
        self.assertEqual(a['source_context']['assay_context'],'stationary_phase')
        self.assertEqual(c['source_context']['assay_context'],'germination')
        self.assertIn('original standard',c['source_context']['identity_evidence'])
        mass=self.report['source_consistency']
        self.assertTrue(mass['table3_caption_labels_disagree'])
        self.assertEqual(mass['table3_theoretical_mz'],mass['figure2_panel_b_acquired_mz'])
        self.assertEqual(mass['table3_acquired_mz'],mass['figure2_panel_b_theoretical_mz'])
        self.assertEqual(c['observation']['Acquired mass (m/z)'],'197.1178')
    def test_stored_report_replays_exactly_except_input_path_location(self):
        stored=json.loads((ROOT/'benchmark/v3/results/germicidin-identity-v1/audit.json').read_text())
        actual=copy.deepcopy(self.report)
        self.assertEqual(sorted(stored.pop('inputs').values()),sorted(actual.pop('inputs').values()))
        self.assertEqual(stored,actual)

if __name__=='__main__':unittest.main(verbosity=2)
