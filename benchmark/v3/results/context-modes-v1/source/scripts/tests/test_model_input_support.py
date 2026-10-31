"""Input coverage is not an accuracy claim; fixtures only test exclusion mechanics."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np

from test_context_model_training import tcm as T, row, tiny_contract


class InputSupportTests(unittest.TestCase):
    def setUp(self):
        self.c = tiny_contract(('fn.type_ec', 'role.catalytic', 'seq.src.other.reported', 'ctx.layout.available'))
        self.rows = [dict(row(i, split, split, split, label, [1, 1, 1, 0]),
                          sequence_sources=['swissprot_diamond_baseline'])
                     for i, (split, label) in enumerate([('train', 0), ('train', 1), ('calibration', 0), ('calibration', 1)])]
        self.art = {'intended_use': 'production', 'contract': {'features': self.c['names']},
                    'input_support': T.build_input_support(self.rows[:2], self.rows[2:], self.c['names'])}

    def test_unknown_tools_pooled_into_same_column_remain_distinct(self):
        x = self.rows[0]['x']
        self.assertIsNone(T.unsupported_input_reason(self.art, x, ['swissprot_diamond_baseline']))
        self.assertIn('unseen', T.unsupported_input_reason(self.art, x, ['unrelated_tool']))
        self.assertIn('unseen', T.unsupported_input_reason(self.art, x, ['clean', 'swissprot_diamond_baseline']))

    def test_legacy_support_preserves_universal_and_refuses_unseen_modes(self):
        x, sources = self.rows[0]['x'], self.rows[0]['sequence_sources']
        explicit = [dict(r, cohort_context_mode='UNIVERSAL') for r in self.rows]
        self.assertEqual(T.build_input_support(explicit[:2], explicit[2:], self.c['names']), self.art['input_support'])
        for mode in ('PATHWAY', 'ORGANISM'):
            self.assertEqual(T.unsupported_input_reason(self.art, x, sources, cohort_context_mode=mode),
                             'unseen_reaction_context_mode')
            with self.assertRaisesRegex(T.ProvenanceError, 'unseen_reaction_context_mode'):
                T.predict_proba(self.art, x, sources, cohort_context_mode=mode)
            with self.assertRaisesRegex(T.ProvenanceError, 'unseen_reaction_context_mode'):
                T.assign_locus(self.art, [x, x], [sources, sources], cohort_context_modes=['UNIVERSAL', mode])

    def test_mode_profiles_cannot_borrow_other_modes_class_coverage(self):
        rows = [dict(r, cohort_context_mode='PATHWAY') for r in self.rows]
        rows[-1]['label'] = 0  # Positive calibration evidence exists only for UNIVERSAL.
        art = copy.deepcopy(self.art)
        art['input_support'] = T.build_input_support(self.rows[:2] + rows[:2], self.rows[2:] + rows[2:], self.c['names'])
        self.assertEqual(art['input_support']['schema'], 'gspa.model-input-support/2')
        x, sources = rows[0]['x'], rows[0]['sequence_sources']
        self.assertIsNone(T.unsupported_input_reason(art, x, sources))
        self.assertEqual(T.unsupported_input_reason(art, x, sources, cohort_context_mode='PATHWAY'),
                         'profile_lacks_both_classes_in_train_and_calibration')
        self.assertEqual(T.unsupported_input_reason(art, x, sources, cohort_context_mode='ORGANISM'),
                         'unseen_reaction_context_mode')

    def test_modes_survive_serialization_and_enter_training_digest(self):
        rows = [dict(r, cohort_context_mode='ORGANISM') for r in self.rows]
        self.assertNotEqual(T.training_data_digest(self.rows[:2], self.rows[2:]),
                            T.training_data_digest(rows[:2], rows[2:]))
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'rows.jsonl'
            T.write_rows_jsonl(rows, self.c['names'], p)
            self.assertEqual([r['cohort_context_mode'] for r in T.read_table(p, self.c)], ['ORGANISM'] * 4)
        for invalid in (None, '', 'pathway', 'TYPO'):
            with self.subTest(mode=invalid), self.assertRaises(T.ProvenanceError):
                T.unsupported_input_reason(self.art, rows[0]['x'], rows[0]['sequence_sources'], cohort_context_mode=invalid)

    def test_evaluator_withholds_unseen_mode_despite_identical_numeric_vectors(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'art.json'
            digest = T.write_frozen(dict(self.art, training={'synthetic': False}), p)
            rows = [dict(r, cohort_context_mode='PATHWAY', split='test') for r in self.rows[:2]]
            got = T.evaluate_frozen(p, self.c, rows, digest)
            self.assertEqual(got['n'], 0)
            self.assertIsNone(got['auroc'])

    def test_quality_guard_matches_jvm_and_blocks_partial_locus_scoring(self):
        x = self.rows[0]['x']
        sources = ['swissprot_diamond_baseline']
        quality = [{'model_id': 'ANF00021'}]
        self.assertEqual(T.unsupported_input_reason(self.art, x, sources, quality),
                         'gene_model_quality_not_supported_by_model')
        with self.assertRaisesRegex(T.ProvenanceError, 'gene_model_quality'):
            T.predict_proba(self.art, x, sources, quality)
        with self.assertRaisesRegex(T.ProvenanceError, 'gene_model_quality'):
            T.assign_locus(self.art, [x, x], [sources, sources], [[], quality])
        self.assertIsNone(T.unsupported_input_reason(self.art, x, sources, []))

    def test_serialization_retains_quality_observations_for_training_rejection(self):
        rows = copy.deepcopy(self.rows)
        rows[0]['quality_evidence'] = [{'model_id': 'ANF00021'}]
        with self.assertRaisesRegex(T.ProvenanceError, 'gene-model quality'):
            T.train(rows[:2], rows[2:], self.c, 'logistic', identity='quality-test')
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'rows.jsonl'
            T.write_rows_jsonl(rows, self.c['names'], p)
            self.assertEqual(json.loads(p.read_text().splitlines()[0])['quality_evidence'], rows[0]['quality_evidence'])
            with self.assertRaisesRegex(T.ProvenanceError, 'gene-model quality'):
                T.read_table(p, self.c)

    def test_role_function_and_availability_changes_are_not_supported(self):
        for pos in (0, 1, 3):
            x = self.rows[0]['x'].copy(); x[pos] = 1-x[pos]
            self.assertIn('unseen', T.unsupported_input_reason(self.art, x, ['swissprot_diamond_baseline']))

    def test_both_classes_are_required_in_both_splits(self):
        for split in ('train', 'calibration'):
            a = copy.deepcopy(self.art)
            a['input_support']['profiles'][0][split + '_labels'] = [0, 4]
            self.assertIn('both_classes', T.unsupported_input_reason(a, self.rows[0]['x'], ['swissprot_diamond_baseline']))

    def test_missing_identity_is_unknown_not_an_empty_source_set(self):
        self.assertIn('identity_unavailable', T.unsupported_input_reason(self.art, self.rows[0]['x']))
        rows = copy.deepcopy(self.rows)
        del rows[0]['sequence_sources']
        self.assertFalse(T.build_input_support(rows[:2], rows[2:], self.c['names'])['available'])
        a = dict(self.art); del a['input_support']
        with self.assertRaisesRegex(T.ProvenanceError, 'model_input_support_unavailable'):
            T.predict_proba(a, self.rows[0]['x'])

    def test_sources_roundtrip_and_enter_fitting_digest(self):
        before = T.training_data_digest(self.rows[:2], self.rows[2:])
        changed = copy.deepcopy(self.rows); changed[0]['sequence_sources'] = ['another_tool']
        self.assertNotEqual(before, T.training_data_digest(changed[:2], changed[2:]))
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'rows.jsonl'
            T.write_rows_jsonl(self.rows, self.c['names'], p)
            got = T.read_table(p, self.c)
            self.assertEqual(got[0]['sequence_sources'], ['swissprot_diamond_baseline'])

    def test_malformed_source_identity_and_indicator_values_fail(self):
        for sources in ('tool', ['X'], ['z','a'], ['x','x']):
            with self.subTest(sources=sources), self.assertRaises(T.ProvenanceError):
                T.validate_sequence_sources(sources)
        bad = copy.deepcopy(self.rows); bad[0]['x'][0] = .5
        with self.assertRaisesRegex(T.ProvenanceError, 'binary'):
            T.build_input_support(bad[:2], bad[2:], self.c['names'])

    def test_unsupported_test_rows_are_counted_without_manufactured_metrics(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'art.json'
            a = dict(self.art, training={'synthetic': False})
            digest = T.write_frozen(a, p)
            rows = [dict(r, sequence_sources=['unknown'], split='test') for r in self.rows[:2]]
            got = T.evaluate_frozen(p, self.c, rows, digest)
            self.assertEqual(got['input_coverage']['total_rows'], 2)
            self.assertEqual(got['n'], 0)
            self.assertIsNone(got['auroc'])


if __name__ == '__main__':
    unittest.main()
