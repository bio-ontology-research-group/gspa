"""Independent metric/population controls; no models or LP solvers."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('context_audit', Path(__file__).resolve().parents[1] / 'audit_context_comparison.py')
A = importlib.util.module_from_spec(spec); spec.loader.exec_module(A)


class ContextAuditTests(unittest.TestCase):
    def rows(self):
        return [dict(genome_id='small' if i == 0 else 'large', split='calibration', label=y, scores={'a': p})
                for i, (p, y) in enumerate([(.9, 1), (.8, 0), (.6, 0), (.5, 1)])]

    def test_threshold_is_calibration_only_and_prefers_conservative_ties(self):
        rows = self.rows()
        self.assertEqual(A.threshold(rows, 'a')['value'], .9)
        rows[0]['split'] = 'test'
        with self.assertRaises(ValueError): A.threshold(rows, 'a')

    def test_all_confusion_cells_missing_predictions_and_unknowns_remain_distinct(self):
        m = A.metrics(self.rows(), 'a', .7)
        self.assertEqual([m[k] for k in ['tp', 'fp', 'fn', 'tn']], [1, 1, 1, 1])
        self.assertAlmostEqual(m['brier'], (.01 + .64 + .36 + .25) / 4)
        for p in (True, float('nan'), float('inf'), -.01, 1.01):
            with self.assertRaises(ValueError): A.probability(p)

    def test_pairing_preserves_harm_and_weights_genomes_equally(self):
        rows = [dict(genome_id=g, label=y, scores={'a': a, 'b': b}) for g, y, a, b in
                [('small', 1, 1., 0.), ('large', 0, 1., 0.), ('large', 0, 1., 0.)]]
        p = A.paired(rows, 'a', 'b', {'a': .5, 'b': .5})
        self.assertEqual((p['reference_corrected'], p['reference_harmed']), (1, 2))
        self.assertEqual(p['mean_genome_brier_delta']['estimate'], 0.)
        self.assertEqual(p['mean_genome_brier_delta'], A.paired(rows, 'a', 'b', {'a': .5, 'b': .5})['mean_genome_brier_delta'])

    def test_full_ec_panel_keeps_unemitted_references_and_reports_ties(self):
        refs = {('g1', 'p'): {'1.1.1.1', '2.2.2.2'}, ('g2', 'missing'): {'3.3.3.3'}}
        scores = {('g1', 'p'): {'1.1.1.1': .8, '0.0.0.0': .8, '2.2.2.2': .4}}
        result = A.ec_panel(refs, scores, .5)
        self.assertEqual((result['proteins'], result['reference_pairs'], result['emitted_reference_pairs']), (2, 3, 2))
        self.assertEqual((result['reference_matches'], result['unconfirmed_calls'], result['missed_reference']), (1, 1, 2))
        self.assertEqual((result['top1_matches'], result['top3_matches'], result['top1_tied']), (0, 1, 1))
        self.assertEqual(result['records'][1]['candidates'], 0)

    def test_delta_resolution_prefers_final_then_base_then_immutable_package(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); package, base, final = (root / x for x in ('package', 'base', 'final'))
            rel = Path('example/file.json'); original = A.PREFIX / rel
            for r in (package / str(A.PREFIX).lstrip('/'), base, final):
                p = r / rel; p.parent.mkdir(parents=True); p.write_text(str(r))
            locator = A.Inputs(package, base, final)
            self.assertEqual(locator.locate(original), final / rel)
            (final / rel).unlink(); self.assertEqual(locator.locate(original), base / rel)
            (base / rel).unlink(); self.assertEqual(locator.locate(original), package / str(original).lstrip('/'))
            with self.assertRaises(ValueError): locator.locate('../other')


if __name__ == '__main__': unittest.main()
