import json
from pathlib import Path
import unittest
from unittest.mock import patch

from benchmark.v3 import verify_price_part as V
from benchmark.v3 import export_price_matrix_job as W
from benchmark.v3.tests.test_prepare_price_context import Base, ASM


class VerifyWholeMatrix(Base):
    def prepare(self):
        self.targets = self.targets_file()
        self.rec = self.part(targets=self.targets)
        self.directory = Path(self.rec['part_dir'])
        return self.directory / 'claims.jsonl'

    def verify(self):
        return V.verify(self.directory, self.layout, self.targets)

    def refresh_receipt(self, claims):
        p = self.directory/'part.json'; rec = json.loads(p.read_text())
        rec['claims'].update(sha256=V.digest(claims), bytes=claims.stat().st_size,
                             claims=len(claims.read_text().splitlines()))
        p.write_text(json.dumps(rec))

    def test_full_matrix_includes_weak_preliminary_unplaced_and_unresolved(self):
        self.prepare()
        report = self.verify()
        self.assertEqual(report['status'], 'verified')
        self.assertEqual(report['matrix_rows'], 9)
        self.assertEqual(report['claims_compared'], 18)
        self.assertEqual(report['proteins_with_claims'], 8)

    def test_rehashed_score_change_still_fails_matrix_comparison(self):
        claims = self.prepare(); rows = claims.read_text().splitlines()
        value = json.loads(rows[0]); value['raw_score'] = .123
        rows[0] = json.dumps(value); claims.write_text('\n'.join(rows)+'\n')
        self.refresh_receipt(claims)
        with self.assertRaisesRegex(ValueError, 'differs from matrix'): self.verify()

    def test_rehashed_loss_and_extra_claim_fail(self):
        claims = self.prepare(); original = claims.read_text()
        for changed, message in [(original.rsplit('\n', 2)[0]+'\n', 'missing claim'),
                                 (original+original.splitlines()[0]+'\n', 'extra claims')]:
            claims.write_text(changed); self.refresh_receipt(claims)
            with self.assertRaisesRegex(ValueError, message): self.verify()

    def test_rehashed_wrong_protein_version_fails(self):
        claims = self.prepare(); claims.write_text(claims.read_text().replace('W1.1', 'W1.2'))
        self.refresh_receipt(claims)
        with self.assertRaisesRegex(ValueError, 'differs from matrix'): self.verify()

    def test_pruned_context_cannot_be_certified_as_whole_matrix(self):
        self.prepare(); p = self.directory/'part.json'; rec = json.loads(p.read_text())
        rec['retention']['top_k'] = 10; p.write_text(json.dumps(rec))
        with self.assertRaisesRegex(ValueError, 'all-nonzero'): self.verify()

    def test_input_mutation_during_validation_is_detected(self):
        self.prepare(); original = V.C.open_matrix
        def changing(*args):
            result = original(*args)
            p = Path(self.targets); p.write_text(p.read_text()+'\n')
            return result
        with patch.object(V.C, 'open_matrix', side_effect=changing):
            with self.assertRaisesRegex(ValueError, 'changed during'): self.verify()


class PublishVerifiedMatrix(Base):
    def setUp(self):
        super().setUp()
        matrix, sha = self.matrix()
        self.plan = self.tmp/'plan.json'
        targets = self.targets_file()
        self.out = self.tmp/'published'
        self.plan.write_text(json.dumps(dict(jobs=[dict(predictor='clean', assembly=ASM,
            matrix=matrix, matrix_sha256=sha)], out=str(self.out), layout=self.layout,
            targets=targets, inputs={targets: V.digest(targets)})))

    def test_publish_receipt_and_refuse_overwrite(self):
        W.run(self.plan, 0, self.tmp/'scratch')
        part = self.out/'parts/clean'/ASM
        receipt = json.loads((part/'job-receipt.json').read_text())
        self.assertEqual(receipt['claims'], 18)
        for name, digest in receipt['artifacts'].items(): self.assertEqual(V.digest(part/name), digest)
        original = V.digest(part/'claims.jsonl')
        with self.assertRaisesRegex(ValueError, 'already exists'): W.run(self.plan, 0, self.tmp/'scratch2')
        self.assertEqual(V.digest(part/'claims.jsonl'), original)

    def test_changed_frozen_input_has_failure_receipt_and_no_part(self):
        plan = json.loads(self.plan.read_text()); Path(plan['targets']).write_text('{}')
        with self.assertRaisesRegex(ValueError, 'frozen input changed'): W.run(self.plan, 0, self.tmp/'scratch')
        self.assertFalse((self.out/'parts/clean'/ASM).exists())
        failure = json.loads(next((self.out/'failures').glob('*.json')).read_text())
        self.assertEqual(failure['status'], 'failed')

    def test_corrupt_copy_never_publishes_completed_part(self):
        original = W.shutil.copytree
        def corrupt(src, dst):
            result = original(src, dst)
            (Path(dst)/'claims.jsonl').write_text('corrupted')
            return result
        with patch.object(W.shutil, 'copytree', side_effect=corrupt):
            with self.assertRaisesRegex(ValueError, 'published copy'): W.run(self.plan, 0, self.tmp/'scratch')
        self.assertFalse((self.out/'parts/clean'/ASM).exists())


if __name__ == '__main__': unittest.main()
