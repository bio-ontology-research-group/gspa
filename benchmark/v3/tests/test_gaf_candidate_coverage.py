import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import assess_gaf_candidate_coverage as C
from test_gaf_negatives import row


class CoverageTests(unittest.TestCase):
    def fixture(self, root):
        r = C.A.negative_record(row(), 1)
        r['candidate_joins'] = [{'genome_id': 'g1', 'protein_id': 'p1',
            'status': 'needs_snapshot_and_reference_review', 'exclusion_reasons': []}]
        candidates = root / 'cohort-candidates.jsonl'
        candidates.write_text(json.dumps(r) + '\n')
        inventory = root / 'inventory.json'
        inventory.write_text(json.dumps({'schema': 'gspa.gaf-negative-inventory/1', 'status': 'verified',
            'output_sha256': {'cohort-candidates.jsonl': C.A.sha(candidates)}}))
        exports = root / 'features' / 'g1'
        exports.mkdir(parents=True)
        rows = exports / 'candidate_features.jsonl'
        rows.write_text(''.join(json.dumps({'genome_id': 'g1', 'protein_id': p,
            'function_type': 'GO', 'function_id': 'GO:0000001', 'role': 'UNSPECIFIED',
            'candidate_key': p + '|GO|GO:0000001|UNSPECIFIED'}) + '\n' for p in ('p1', 'p2')))
        (exports / 'export_manifest.json').write_text(json.dumps({'schema': 'gspa.context-feature-export/1',
            'rows_file': rows.name, 'rows_sha256': C.A.sha(rows), 'rows': 2}))
        return inventory, root / 'features', rows

    def test_exact_match_retains_competitors_without_creating_truth(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            inv, features, _ = self.fixture(d)
            report = C.assess(inv, features, d / 'out')
            self.assertEqual(report['candidate_join_counts'], {'exact_candidate_found_needs_review': 1})
            self.assertEqual(report['accepted_training_labels'], 0)
            q = json.loads((d / 'out/source-review-queue.jsonl').read_text())
            self.assertEqual(q['candidates'][0]['competing_proteins'], 2)
            self.assertFalse(q['accepted_training_label'])

    def test_changed_features_or_candidate_inventory_fail(self):
        for target in ('features', 'inventory'):
            with self.subTest(target=target), tempfile.TemporaryDirectory() as d:
                d = Path(d)
                inv, features, rows = self.fixture(d)
                path = rows if target == 'features' else d / 'cohort-candidates.jsonl'
                with path.open('a') as fh:
                    fh.write('\n')
                with self.assertRaisesRegex(ValueError, 'hash mismatch'):
                    C.assess(inv, features, d / 'out')
                self.assertFalse((d / 'out/coverage.json').exists())

    def test_excluded_context_cannot_become_a_candidate_match(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            inv, features, _ = self.fixture(d)
            path = d / 'cohort-candidates.jsonl'
            r = json.loads(path.read_text())
            r['candidate_joins'][0].update(status='excluded', exclusion_reasons=['interaction_context'])
            path.write_text(json.dumps(r) + '\n')
            m = json.loads(inv.read_text())
            m['output_sha256'][path.name] = C.A.sha(path)
            inv.write_text(json.dumps(m))
            report = C.assess(inv, features, d / 'out')
            self.assertEqual(report['candidate_join_counts'], {'excluded_by_inventory': 1})
            self.assertEqual(report['matched_relation_counts'], {})

    def test_other_terms_and_unexported_genomes_are_not_exact_matches(self):
        for change, expected in [('term', 'no_exact_candidate'), ('genome', 'genome_not_exported')]:
            with self.subTest(change=change), tempfile.TemporaryDirectory() as d:
                d = Path(d)
                inv, features, _ = self.fixture(d)
                path = d / 'cohort-candidates.jsonl'
                r = json.loads(path.read_text())
                if change == 'term':
                    r['go_id'] = 'GO:0000002'
                else:
                    r['candidate_joins'][0]['genome_id'] = 'g2'
                path.write_text(json.dumps(r) + '\n')
                m = json.loads(inv.read_text())
                m['output_sha256'][path.name] = C.A.sha(path)
                inv.write_text(json.dumps(m))
                report = C.assess(inv, features, d / 'out')
                self.assertEqual(report['candidate_join_counts'], {expected: 1})

    def test_incomplete_inventory_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            inv, features, _ = self.fixture(d)
            m = json.loads(inv.read_text())
            m['status'] = 'partial'
            inv.write_text(json.dumps(m))
            with self.assertRaisesRegex(ValueError, 'verified complete inventory'):
                C.assess(inv, features, d / 'out')


if __name__ == '__main__':
    unittest.main()
