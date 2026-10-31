import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import verify_price_export_collection as V


class CollectionIntegrity(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.plan_path = self.root / 'plan.json'
        self.matrix = self.root / 'matrix.bin'
        self.matrix.write_bytes(b'original matrix')
        self.job = dict(predictor='p', assembly='g', matrix=str(self.matrix),
                        matrix_sha256=V.digest(self.matrix), expected_claims=2, matrix_rows=1)
        self.plan = dict(out=str(self.root), jobs=[self.job], inputs={}, expected_claims=2)
        self.save(self.plan_path, self.plan)
        self.part = self.root / 'parts/p/g'
        self.part.mkdir(parents=True)
        for name in ('claims.jsonl', 'join_map.tsv', 'part.json'):
            (self.part / name).write_text(name)
        self.verification = dict(schema='gspa.price-part-matrix-verification/1', status='verified',
            predictor='p', assembly='g', claims_compared=2, matrix_rows=1,
            inputs={str(self.matrix): V.digest(self.matrix), **{str(self.root / 'scratch' / n): V.digest(self.part / n)
                    for n in ('claims.jsonl', 'join_map.tsv', 'part.json')}})
        self.save(self.part / 'matrix-verification.json', self.verification)
        self.receipt = dict(schema='gspa.price-full-export-job/1', status='verified', index=0,
            job=self.job, claims=2, matrix_rows=1,
            inputs={str(self.plan_path): V.digest(self.plan_path), str(self.matrix): V.digest(self.matrix)},
            artifacts={p.name: V.digest(p) for p in self.part.iterdir()})
        self.save(self.part / 'job-receipt.json', self.receipt)

    @staticmethod
    def save(path, value):
        path.write_text(json.dumps(value))

    def test_complete_collection_rehashes_bulk_outputs_and_sources(self):
        result = V.verify(self.plan_path)
        self.assertEqual((result['n_parts'], result['claims']), (1, 2))
        self.assertIn(str(self.part / 'claims.jsonl'), result['hashes'])
        self.assertIn(str(self.matrix), result['hashes'])

    def test_changed_bulk_output_fails_even_with_success_receipts(self):
        (self.part / 'claims.jsonl').write_text('changed')
        with self.assertRaisesRegex(ValueError, 'hash changed'):
            V.verify(self.plan_path)

    def test_changed_original_matrix_fails(self):
        self.matrix.write_bytes(b'changed')
        with self.assertRaisesRegex(ValueError, 'hash changed'):
            V.verify(self.plan_path)

    def test_missing_and_unplanned_parts_fail(self):
        self.part.rename(self.part.with_name('unexpected'))
        with self.assertRaisesRegex(ValueError, 'coverage'):
            V.verify(self.plan_path)

    def test_duplicate_planned_parts_fail(self):
        self.plan['jobs'].append(self.job)
        self.save(self.plan_path, self.plan)
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            V.verify(self.plan_path)

    def test_wrong_success_receipt_and_count_fail(self):
        for field, value in [('status', 'failed'), ('index', 1), ('claims', 1)]:
            receipt = {**self.receipt, field: value}
            self.save(self.part / 'job-receipt.json', receipt)
            with self.assertRaises(ValueError):
                V.verify(self.plan_path)

    def test_matrix_verification_must_bind_the_published_claims(self):
        self.verification['inputs'][str(self.root / 'scratch/claims.jsonl')] = '0' * 64
        self.save(self.part / 'matrix-verification.json', self.verification)
        self.receipt['artifacts']['matrix-verification.json'] = V.digest(self.part / 'matrix-verification.json')
        self.save(self.part / 'job-receipt.json', self.receipt)
        with self.assertRaisesRegex(ValueError, 'published artifact differ'):
            V.verify(self.plan_path)


if __name__ == '__main__':
    unittest.main()
