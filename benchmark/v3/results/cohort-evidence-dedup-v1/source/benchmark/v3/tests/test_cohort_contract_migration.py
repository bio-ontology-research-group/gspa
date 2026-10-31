import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import migrate_cohort_feature_contract as M


class CohortContractMigration(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def write(self, name, rows):
        path = self.root / name
        path.write_text(''.join(json.dumps(row) + '\n' for row in rows))
        return path

    def pair(self):
        old = dict(candidate_key='p|EC|1.2.3.4|CATALYTIC', genome_id='g', sequence_sources=['diamond'],
                   features={**dict.fromkeys(M.REMOVED, 1.), 'seq.max_raw_any': .4, 'cohort.w.perm_p': .01})
        new = copy.deepcopy(old)
        new['features'] = {k: v for k, v in old['features'].items() if k not in M.REMOVED}
        new['features']['cohort.w.perm_p'] = .1
        return old, new

    def test_only_declared_columns_are_removed_and_noncohort_vectors_are_preserved(self):
        old, new = self.pair()
        count, changed = M.compare_vectors(self.write('old', [old]), self.write('new', [new]))
        self.assertEqual(count, 1)
        self.assertEqual(changed, {'cohort.w.perm_p': 1})
        new['features']['seq.max_raw_any'] = .5
        with self.assertRaisesRegex(ValueError, 'non-cohort'):
            M.compare_vectors(self.root / 'old', self.write('new', [new]))

    def test_candidate_identity_source_and_column_coverage_changes_fail(self):
        old, original = self.pair()
        previous = self.write('old', [old])
        for key, value in [('candidate_key', 'different'), ('genome_id', 'other'), ('sequence_sources', ['clean'])]:
            new = copy.deepcopy(original); new[key] = value
            with self.assertRaisesRegex(ValueError, 'identity or source'):
                M.compare_vectors(previous, self.write('new', [new]))
        with self.assertRaisesRegex(ValueError, 'universe'):
            M.compare_vectors(previous, self.write('new', []))
        new = copy.deepcopy(original); new['features']['unexpected'] = 0
        with self.assertRaisesRegex(ValueError, 'column'):
            M.compare_vectors(previous, self.write('new', [new]))

    def test_supervision_comparison_preserves_polarity_groups_splits_and_sources(self):
        old = dict(candidate_key='p', features={'old': 1}, label=0, label_basis='explicit_negative_truth',
                   group='family', locus_group_id='locus', split='train', sequence_sources=['diamond'])
        new = {**old, 'features': {'new': 2}}
        previous = self.write('old', [old])
        self.assertEqual(M.compare_supervision(previous, self.write('new', [new])), 1)
        for key, value in [('label', 1), ('label_basis', 'closed_world_unconfirmed'), ('group', 'other'),
                           ('locus_group_id', None), ('split', 'test'), ('sequence_sources', ['clean'])]:
            with self.assertRaisesRegex(ValueError, 'supervision'):
                M.compare_supervision(previous, self.write('new', [{**new, key: value}]))
        with self.assertRaisesRegex(ValueError, 'coverage'):
            M.compare_supervision(previous, self.write('new', []))

    def test_changed_sources_are_rejected_before_reuse(self):
        source = self.write('source', [{'value': 1}]); identity = M.hashes([source])
        M.check_hashes(identity)
        self.write('source', [{'value': 2}])
        with self.assertRaisesRegex(ValueError, 'changed'):
            M.check_hashes(identity)

    def test_same_contract_replay_changes_only_cohort_values_without_dropping_columns(self):
        _, old = self.pair()
        new = copy.deepcopy(old)
        new['features']['cohort.w.perm_p'] = .5
        previous = self.write('old', [old])
        self.assertEqual(M.compare_vectors(previous, self.write('new', [new]), set()),
                         (1, {'cohort.w.perm_p': 1}))
        new['features'].pop('cohort.w.perm_p')
        with self.assertRaisesRegex(ValueError, 'column'):
            M.compare_vectors(previous, self.write('new', [new]), set())


if __name__ == '__main__':
    unittest.main()
