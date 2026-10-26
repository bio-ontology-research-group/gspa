import unittest
import numpy as np
from benchmark.v3.audit_price_matrices import check_matrix


class MatrixJoinTest(unittest.TestCase):
    def test_versionless_unique_join_and_zero_scores_not_candidates(self):
        r = check_matrix(['WP_1', 'WP_2.2'], ['EC:1.1.1.1', '2.7.-.-', '3.1.1.1'],
                         np.array([[.7, .9, 0], [.2, 0, 1e-8]]), {'WP_1.1', 'WP_2.2'}, {'WP_2.2'})
        self.assertEqual(r['versionless_unique_matches'], {'WP_1': 'WP_1.1'})
        self.assertEqual(r['nonzero_complete_ec_scores'], 3)
        self.assertEqual(r['nonzero_at_least_0_01'], 2)
        self.assertEqual(r['matched_proteins_with_single_locus'], 1)

    def test_ambiguous_versions_and_wrong_explicit_version_never_guessed(self):
        r = check_matrix(['WP_1', 'WP_2.1'], ['1.1.1.1'], np.ones((2, 1)),
                         {'WP_1.1', 'WP_1.2', 'WP_2.2'}, set())
        self.assertEqual(r['ambiguous_versionless_ids'], ['WP_1'])
        self.assertEqual(r['absent_protein_ids'], ['WP_2.1'])

    def test_duplicate_alias_rows_and_bad_scores_refused(self):
        with self.assertRaisesRegex(ValueError, 'same acquired protein'):
            check_matrix(['p', 'p.1'], ['1.1.1.1'], np.ones((2, 1)), {'p.1'}, set())
        for score in [float('nan'), float('inf'), -1, 2]:
            with self.assertRaisesRegex(ValueError, 'score'):
                check_matrix(['p.1'], ['1.1.1.1'], np.array([[score]]), {'p.1'}, set())

    def test_retention_diagnostic_does_not_silently_count_partial_ecs(self):
        r = check_matrix(['p.1'], [f'1.1.1.{i}' for i in range(1, 151)] + ['1.1'],
                         np.array([[.001] * 150 + [.9]]), {'p.1'}, {'p.1'})
        self.assertEqual(r['nonzero_complete_ec_scores'], 150)
        self.assertEqual(r['threshold_0_01_plus_top100_count'], 100)
