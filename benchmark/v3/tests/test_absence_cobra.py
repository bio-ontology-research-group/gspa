"""Independent LP comparison and isolation checks for the absence analysis."""
import copy
import json
import sys
from pathlib import Path
import tempfile
import unittest

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / 'benchmark/v3'))
sys.path.insert(0, str(REPO / 'scripts/tests'))
try:
    import check_absence_cobra as O
    from test_absence_experiments import problem
    from test_metabolic_experiments import write_sbml, ME, REACTIONS
except ImportError:
    O = None


@unittest.skipIf(O is None, 'COBRApy, SciPy and libSBML required')
class AbsenceOracleTests(unittest.TestCase):
    def test_independent_solver_matches_complete_fragmented_and_supply_controls(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_sbml(root / 'model.xml')
            model = ME.load_sbml_model(root / 'model.xml')
            for variation in ('complete', 'fragmented', 'no_supply', 'supplied_medium', 'infeasible_maintenance'):
                m = problem()
                m['model'] = 'model.xml'
                if variation == 'fragmented':
                    m['absence_comparison']['candidate_function_key'] = None
                    m['absence_comparison']['assembly']['status'] = 'DRAFT'
                elif variation == 'no_supply':
                    m['absence_comparison']['external_supply'] = None
                elif variation == 'supplied_medium':
                    m['baseline']['medium']['EX_vit_e'] = 2
                elif variation == 'infeasible_maintenance':
                    m['absence_comparison']['conditions'][0]['mods'] = [{'op': 'knockout', 'reactions': ['ATPM']}]
                manifest, result = root / 'manifest.json', root / 'result.json'
                manifest.write_text(json.dumps(m))
                result.write_text(json.dumps(ME.run_experiments(m, model)))
                report = O.compare(manifest, result)
                self.assertTrue(report['agrees'], (variation, report))
                self.assertGreater(len(report['comparisons']), 0)

    def test_oracle_restores_model_across_transport_knockout_and_reversed_signs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reactions = copy.deepcopy(REACTIONS)
            reactions['T_vit_alternative'] = reactions['T_vit']
            for rid in ('EX_vit_e', 'T_vit', 'T_vit_alternative'):
                a, b, lo, hi, gene = reactions[rid]
                reactions[rid] = b, a, -hi, -lo, gene
            path = root / 'model.xml'
            write_sbml(path, reactions=reactions)
            model = O.cobra.io.read_sbml_model(path, f_replace={})
            m = problem()
            before = {r.id: r.bounds for r in model.reactions}
            ko = m['absence_comparison']['conditions'][2]
            control = m['absence_comparison']['conditions'][0]
            self.assertEqual(O.simulate(model, m, 'external_supply', ko)['growth'], 0)
            self.assertGreater(O.simulate(model, m, 'external_supply', control)['growth'], 0)
            self.assertEqual({r.id: r.bounds for r in model.reactions}, before)

    def test_checker_rejects_changed_model_or_assumptions(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / 'model.xml'
            write_sbml(path)
            m = problem()
            m['model'] = 'model.xml'
            result = ME.run_experiments(m, ME.load_sbml_model(path))
            manifest, output = root / 'manifest.json', root / 'result.json'
            output.write_text(json.dumps(result))
            m['absence_comparison']['external_supply']['uptake_limit'] = 99
            manifest.write_text(json.dumps(m))
            with self.assertRaisesRegex(ValueError, 'assumptions changed'):
                O.compare(manifest, output)
