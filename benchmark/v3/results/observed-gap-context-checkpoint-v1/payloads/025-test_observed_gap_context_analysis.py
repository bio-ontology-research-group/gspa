import copy
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import analyze_observed_gap_context_models as C


def record(i, label, split='test', genome='g'):
    return dict(example_id=str(i), protein_id='p|' + str(i), genome_id=genome,
        function_type='EC', function_id='EC:1.1.1.1', split=split, label=label,
        label_basis='positive_truth' if label else 'closed_world_unconfirmed')


class ContextAnalysisTest(unittest.TestCase):
    def test_complete_population_and_normalized_identities(self):
        rows = [record(0, 1), record(1, 0)]
        scores = {a: {'g': {C.A.function_key(r): .8 if r['label'] else .2 for r in rows}} for a in C.ARMS}
        attached = C.attach_scores(rows, scores)
        cuts = {a: {'value': .5} for a in C.ARMS}
        result = C.evaluate_rows(attached, cuts, 10)
        self.assertEqual(result['label_basis'], {'positive_truth': 1, 'closed_world_unconfirmed': 1})
        for policy in result['policies'].values():
            self.assertEqual(len(policy['arms']), 9)
            self.assertEqual(len(policy['comparisons']), 12)
            for arm in policy['arms'].values():
                self.assertEqual((arm['n'], arm['tp'], arm['tn']), (2, 1, 1))
        for broken in ('missing_row', 'missing_arm', 'nan', 'bool'):
            s = copy.deepcopy(scores)
            arm = C.ARMS[-1]; key = C.A.function_key(rows[0])
            if broken == 'missing_row': del s[arm]['g'][key]
            elif broken == 'missing_arm': del s[arm]
            else: s[arm]['g'][key] = float('nan') if broken == 'nan' else True
            with self.subTest(broken=broken), self.assertRaises(ValueError): C.attach_scores(rows, s)
        with self.assertRaisesRegex(ValueError, 'duplicate'): C.attach_scores(rows + rows[:1], scores)

    def test_thresholds_only_use_calibration_and_conservative_ties(self):
        rows = [dict(record(i, y, 'calibration'), scores={a: p for a in C.ARMS})
            for i, (y, p) in enumerate([(1, .9), (0, .8), (0, .6), (1, .5)])]
        cuts = C.select_thresholds(rows)
        self.assertTrue(all(c['value'] == .9 for c in cuts.values()))
        with self.assertRaisesRegex(ValueError, 'test rows only'): C.evaluate_rows(rows, cuts, 10)
        rows[0]['split'] = 'test'
        with self.assertRaisesRegex(ValueError, 'calibration rows only'): C.select_thresholds(rows)

    def test_context_contrasts_keep_scorers_matched_and_include_harms(self):
        rows = [dict(record(0, 1, genome='small'), scores={a: .2 for a in C.ARMS}),
                dict(record(1, 0, genome='large'), scores={a: .1 for a in C.ARMS}),
                dict(record(2, 0, genome='large'), scores={a: .1 for a in C.ARMS})]
        for r in rows: r['scores']['PATHWAY/logistic'] = .9
        result = C.evaluate_rows(rows, {a: {'value': .5} for a in C.ARMS}, 100)
        comparisons = result['policies']['fixed_0.5']['comparisons']
        pair = comparisons['PATHWAY/logistic_vs_UNIVERSAL/logistic']
        self.assertEqual((pair['reference_corrected'], pair['reference_harmed']), (1, 2))
        self.assertAlmostEqual(pair['mean_genome_brier_delta']['estimate'], .085)
        self.assertNotIn('PATHWAY/logistic_vs_UNIVERSAL/gbdt', comparisons)
        self.assertIn('PATHWAY/logistic_vs_PATHWAY/rule', comparisons)

    def test_all_thresholds_are_saved_before_first_test_metric(self):
        with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR', str(Path.home() / 'tmp'))) as d:
            out = Path(d); plans = {}; scores = {}
            for transfer in C.R.G.MODES:
                table = out / (transfer + '.jsonl')
                rows = [record(0, 1, 'calibration'), record(1, 0, 'calibration'), record(2, 1)]
                table.write_text(''.join(json.dumps(r) + '\n' for r in rows))
                plans[transfer] = dict(table=str(table), paths={a: {'g': 'unused'} for a in C.ARMS})
                scores.update({C.A.function_key(r): .9 if r['label'] else .1 for r in rows})
            C.A.save(out / 'protocol.json', dict(C.POLICY, plans=plans))
            C.A.save(out / 'registration.json', dict(inputs={}, protocol_sha256=C.F.T.sha256_file(out / 'protocol.json')))
            original = C.evaluate_rows; calls = []
            def evaluate(rows, cuts):
                frozen = json.loads((out / 'thresholds.json').read_text())
                self.assertEqual(set(frozen), set(C.R.G.MODES))
                self.assertTrue(all(set(v) == set(C.ARMS) for v in frozen.values()))
                calls.append(1)
                return original(rows, cuts, 10)
            with patch.object(C.A, 'read_scores', return_value=scores), patch.object(C, 'evaluate_rows', side_effect=evaluate):
                C.run(out)
                with self.assertRaises(FileExistsError): C.run(out)
            self.assertEqual(len(calls), 2)
            self.assertEqual(json.loads((out / 'verification.json').read_text())['status'], 'verified')


if __name__ == '__main__': unittest.main()
