import math
import os
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import ec_reference_panel as E


class FullEcReferenceTests(unittest.TestCase):
    def test_complete_population_missing_functions_and_unconfirmed_calls(self):
        refs = {('g', 'multi'): {'1.1.1.1', '2.2.2.2'}, ('g', 'missing'): {'3.3.3.3'}}
        scores = {'g': {('multi', 'EC', '1.1.1.1'): .8, ('multi', 'EC', '9.9.9.9'): .9,
                        ('multi', 'GO', 'GO:1'): 1.0, ('unlabelled', 'EC', '4.4.4.4'): 1.0}, 'empty_reference_genome': {}}
        result = E.evaluate(refs, scores, {'g', 'empty_reference_genome'}, .5)
        self.assertEqual((result['proteins'], result['reference_pairs'], result['emitted_reference_pairs']), (2, 3, 1))
        self.assertEqual((result['reference_matches'], result['unconfirmed_calls'], result['missed_reference']), (1, 1, 2))
        self.assertEqual(result['reference_pairs_by_genome'], {'g': 3, 'empty_reference_genome': 0})
        self.assertEqual((result['reference_recall'], result['reference_micro_f1']), (1/3, .4))
        self.assertEqual((result['top1_matches'], result['top3_matches']), (0, 1))

    def test_lexical_ties_and_unchanged_cutoff_do_not_turn_missing_into_zero_scores(self):
        refs = {('g', 'p'): {'2.2.2.2'}, ('g', 'absent'): {'3.3.3.3'}}
        scores = {'g': {('p', 'EC', '2.2.2.2'): .2, ('p', 'EC', '1.1.1.1'): .2}}
        result = E.evaluate(refs, scores, {'g'}, .5)
        self.assertEqual((result['reference_matches'], result['unconfirmed_calls']), (0, 0))
        self.assertEqual((result['top1_matches'], result['top3_matches'], result['top1_tied']), (0, 1, 1))
        zero = E.evaluate(refs, scores, {'g'}, 0)
        self.assertEqual(zero['reference_matches'], 1)
        self.assertEqual(zero['missed_reference'], 1)

    def test_no_reference_is_unavailable_not_perfect_or_wrong(self):
        result = E.evaluate({}, {'g': {('p', 'EC', '1.1.1.1'): 1}}, {'g'}, .5)
        self.assertEqual(result['reference_pairs'], 0)
        self.assertIsNone(result['reference_recall'])
        self.assertIsNone(result['reference_micro_f1'])

    def test_missing_genome_and_invalid_probabilities_fail(self):
        refs = {('g', 'p'): {'1.1.1.1'}}
        with self.assertRaisesRegex(ValueError, 'missing test-genome'):
            E.evaluate(refs, {}, {'g'}, .5)
        for score in (True, math.nan, math.inf, -.1, 1.1):
            with self.subTest(score=score), self.assertRaisesRegex(ValueError, 'invalid probability'):
                E.evaluate(refs, {'g': {('p', 'EC', '1.1.1.1'): score}}, {'g'}, .5)

    def test_reference_population_uses_evidence_and_split_not_candidate_presence(self):
        with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR', str(Path.home()/'tmp'))) as folder:
            p = Path(folder); labels=p/'labels.tsv'; splits=p/'splits.tsv'
            splits.write_text('genome_id\tsplit\ng\ttest\nt\ttrain\ne\ttest\n')
            header='genome_id\tprotein_id\tfunction_id\tevidence\tlabel\n'
            body='g\tp\tEC:1.1.1.1\tECO:0000269\t1\n'
            labels.write_text(header+body+body+'g\tp\tEC:2.2.2.2\tECO:0000305\t1\n'
                              't\tq\tEC:3.3.3.3\tECO:0000269\t1\n')
            refs, genomes=E.load_references(labels, splits, ['g', 't', 'e'])
            self.assertEqual(refs, {('g','p'): {'1.1.1.1'}})
            self.assertEqual(genomes, {'g','e'})
            with self.assertRaisesRegex(ValueError, 'population differs'):
                E.load_references(labels, splits, ['g'])
            labels.write_text(header+body+body.replace('\t1\n', '\t0\n'))
            with self.assertRaisesRegex(ValueError, 'contradictory'):
                E.load_references(labels, splits, ['g','t','e'])


if __name__ == '__main__':
    unittest.main()
