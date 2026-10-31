"""Fixture tests for the source x context evaluation summary; synthetic JSON only."""
import copy
import csv
import json
from pathlib import Path
import tempfile
import unittest

import summarize_source_context_evaluation as S


def cell(supported, total, auroc=.75, reasons=None):
    reasons = reasons if reasons is not None else ({'profile_lacks_both_classes_in_train_and_calibration': total - supported}
                                                   if total > supported else {})
    cov = dict(total_rows=total, supported_rows=supported, unsupported_reasons=reasons,
               locus_groups_with_unsupported_members=0)
    if supported == 0:
        return dict(status=S.NO_ROWS, n=0, input_coverage=cov, auroc=None, logloss=None, brier=None, ece_10bin=None,
                    evaluated_on='test after freezing', locus=dict(trained=False, reason='no supported test rows'))
    return dict(auroc=auroc, logloss=.5, brier=.2, ece_10bin=.1, n=supported, positive_rate=.25,
                at_display_threshold=dict(threshold=.4, precision=.6, recall=.3, recall_all_test_positives=.3,
                                          n_predicted_positive=5),
                input_coverage=cov, locus=dict(trained=False, reason='few groups'), evaluated_on='test',
                metric_scope='supported rows only', probability_target='t', model_sha256='m', source_file_sha256='f')


def evaluation():
    results = {s: {c: {k: cell(40 + i, 50 + i) for k in S.KINDS} for c in S.CONTEXTS}
               for i, s in enumerate(S.SOURCES)}
    results['clean']['SEQONLY']['gbdt'] = cell(7, 7, auroc=.123)
    results['clean']['ORGANISM']['logistic'] = cell(0, 3)
    return dict(scope='Conditional metrics only.', limitations='Reused cohort.',
                pending_stages={'production_replay': 'not_run'}, results=results)


class SummaryTests(unittest.TestCase):
    def run_summary(self, data):
        d = Path(self.tmp.name)
        (d / 'evaluation.json').write_text(json.dumps(data))
        rows = S.summarize(d / 'evaluation.json', d / 'out')
        return rows, (d / 'out' / S.TSV_NAME), (d / 'out' / S.MD_NAME).read_text()

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def test_tsv_has_one_row_per_cell_with_counts(self):
        _, tsv, _ = self.run_summary(evaluation())
        with open(tsv) as h:
            rows = list(csv.DictReader(h, delimiter='\t'))
        self.assertEqual(len(rows), len(S.SOURCES) * len(S.CONTEXTS) * len(S.KINDS))
        self.assertEqual(list(rows[0]), S.COLUMNS)
        r = next(x for x in rows if (x['source'], x['context'], x['kind']) == ('diamond_homology', 'PATHWAY', 'gbdt'))
        self.assertEqual((r['total_rows'], r['supported_rows'], r['unsupported_rows'], r['n_positive_supported']),
                         ('50', '40', '10', '10'))
        self.assertEqual(r['unsupported_reasons'], 'profile_lacks_both_classes_in_train_and_calibration=10')
        empty = next(x for x in rows if (x['source'], x['context'], x['kind']) == ('clean', 'ORGANISM', 'logistic'))
        self.assertEqual((empty['status'], empty['auroc'], empty['precision']), (S.NO_ROWS, '', ''))

    def test_markdown_tables_have_source_rows_context_columns_and_counts(self):
        _, _, md = self.run_summary(evaluation())
        for title in ('AUROC', 'Log-loss', 'Precision / recall at the display threshold'):
            for kind in S.KINDS:
                lines = md.split(f'## {title}: {kind}\n\n')[1].split('\n\n')[0].splitlines()
                self.assertEqual(len(lines), 2 + len(S.SOURCES))
                self.assertEqual(lines[0], '| source | SEQONLY | UNIVERSAL | PATHWAY | ORGANISM |')
                self.assertEqual([l.split(' | ')[0][2:] for l in lines[2:]], list(S.SOURCES))
                self.assertTrue(all(l.count('|') == 2 + len(S.CONTEXTS) for l in lines))
        self.assertIn('| clean | 0.750 (42/52) | 0.750 (42/52) | 0.750 (42/52) | n/a (0/3) |', md)
        self.assertIn('| clean | 0.123 (7/7) |', md)
        self.assertIn('0.600 / 0.300 (5 pred.+; 40/50)', md)
        self.assertIn('| diamond_homology | PATHWAY | gbdt | 50 | 40 | 10 |', md)
        self.assertIn('- scope: "Conditional metrics only."', md)
        self.assertIn('- limitations: "Reused cohort."', md)
        self.assertIn('`production_replay`: not_run', md)

    def test_output_is_deterministic(self):
        _, tsv, md = self.run_summary(evaluation())
        first = tsv.read_bytes()
        _, tsv, md2 = self.run_summary(evaluation())
        self.assertEqual((first, md), (tsv.read_bytes(), md2))

    def test_missing_cell_fails(self):
        for path in (('clean',), ('clean', 'PATHWAY'), ('clean', 'PATHWAY', 'gbdt')):
            data = evaluation()
            node = data['results']
            for key in path[:-1]:
                node = node[key]
            del node[path[-1]]
            with self.assertRaisesRegex(S.SummaryError, 'missing'):
                self.run_summary(data)

    def test_unexpected_cell_and_incomplete_cell_fail(self):
        data = evaluation()
        data['results']['clean']['EXTRA'] = copy.deepcopy(data['results']['clean']['SEQONLY'])
        with self.assertRaisesRegex(S.SummaryError, 'unexpected'):
            self.run_summary(data)
        data = evaluation()
        del data['results']['clean']['SEQONLY']['logistic']['at_display_threshold']
        with self.assertRaisesRegex(S.SummaryError, 'at_display_threshold'):
            self.run_summary(data)
        data = evaluation()
        data['results']['clean']['SEQONLY']['logistic']['input_coverage']['unsupported_reasons'] = {'x': 99}
        with self.assertRaisesRegex(S.SummaryError, 'do not sum'):
            self.run_summary(data)


if __name__ == '__main__':
    unittest.main()
