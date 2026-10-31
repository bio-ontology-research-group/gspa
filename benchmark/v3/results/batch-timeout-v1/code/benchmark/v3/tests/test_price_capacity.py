"""Exact output-retention checks, independent of biological scoring/accuracy."""
import csv
import json
import os
from pathlib import Path
import sys
import tempfile
import tracemalloc
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import run_price_matrix_capacity as C
import verify_candidate_universe as U


class CandidateUniverseTest(unittest.TestCase):
    def setUp(self):
        scratch = Path(os.environ.get('TMPDIR') or Path.home()/'tmp')
        scratch.mkdir(parents=True, exist_ok=True)
        tmp = tempfile.TemporaryDirectory(dir=scratch)
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.claims = self.root/'claims.jsonl'
        self.output = self.root/'output.tsv'
        self.keys = [('weak.1', 'EC', '1.1.1.n11'), ('unplaced.2', 'EC', '2.4.1.129'),
                     ('unicode-μ\tquoted"', 'GO', 'GO:0000001')]

    def prepare(self, keys=None):
        with self.claims.open('w') as stream:
            for p, t, f in self.keys if keys is None else keys:
                stream.write(json.dumps(dict(protein_id=p, function_type=t, function_id=f,
                                             raw_score=1e-45))+'\n')
        return U.expected_candidates(self.claims, self.root)

    def predictions(self, keys=None):
        with self.output.open('w') as stream:
            writer = csv.writer(stream, delimiter='\t')
            writer.writerow(['protein_id', 'type', 'function_id', 'posterior_prob'])
            for p, t, f in reversed(self.keys) if keys is None else keys:
                writer.writerow([p, t, ('EC:'+f) if t == 'EC' else f, .5])

    def test_exact_unordered_comparison_preserves_versions_and_preliminary_ids(self):
        expected, count = self.prepare()
        self.predictions()
        self.assertEqual(count, 3)
        self.assertEqual(U.compare_candidates(expected, self.output, self.root), 3)

    def test_missing_extra_and_same_count_wrong_identity_are_rejected(self):
        expected, _ = self.prepare()
        for keys in (self.keys[:-1], self.keys+[('extra', 'EC', '1.1.1.1')],
                     [('weak.2', 'EC', '1.1.1.n11')]+self.keys[1:]):
            self.predictions(keys)
            with self.assertRaisesRegex(ValueError, 'universe differs'):
                U.compare_candidates(expected, self.output, self.root)

    def test_duplicate_input_and_duplicate_output_are_rejected(self):
        with self.assertRaisesRegex(ValueError, 'duplicate candidate'):
            self.prepare(self.keys+self.keys[:1])
        expected, _ = self.prepare()
        self.predictions(self.keys+self.keys[:1])
        with self.assertRaisesRegex(ValueError, 'duplicate candidate|universe differs'):
            U.compare_candidates(expected, self.output, self.root)

    def test_only_ec_prefix_is_normalized(self):
        expected, _ = self.prepare([('p', 'GO', 'EC:GO:0000001')])
        self.predictions([('p', 'GO', 'GO:0000001')])
        with self.assertRaisesRegex(ValueError, 'universe differs'):
            U.compare_candidates(expected, self.output, self.root)

    def test_malformed_identity_and_header_fail(self):
        with self.assertRaisesRegex(ValueError, 'nonempty strings'):
            self.prepare([('', 'EC', '1.1.1.1')])
        expected, _ = self.prepare()
        self.output.write_text('garbage\n')
        with self.assertRaisesRegex(ValueError, 'header'):
            U.compare_candidates(expected, self.output, self.root)

    def test_large_input_and_output_are_not_materialized_in_python(self):
        n = 30000
        tracemalloc.start()
        try:
            expected, count = self.prepare((f'p{i}', 'EC', '1.1.1.1') for i in range(n))
            self.predictions((f'p{i}', 'EC', '1.1.1.1') for i in reversed(range(n)))
            self.assertEqual(U.compare_candidates(expected, self.output, self.root), count)
            self.assertEqual(count, n)
            _, peak = tracemalloc.get_traced_memory()
        finally:
            tracemalloc.stop()
        self.assertLess(peak, 3*1024*1024)

    def test_capacity_failure_retains_completed_arm_and_cannot_overwrite(self):
        out = self.root/'capacity'
        def failing(part, layout, jar, reference, output, scratch, reports, heap, arms, timeout):
            reports['stage0'] = {'candidate_universe_preserved': True}
            raise RuntimeError('stage1 failed')
        with patch.object(C, '_run', side_effect=failing):
            with self.assertRaisesRegex(RuntimeError, 'stage1 failed'):
                C.run(self.root, self.root, self.root, self.root, out, self.root)
        record = json.loads((out/'verification.json').read_text())
        self.assertEqual(record['status'], 'failed')
        self.assertEqual(record['runs']['stage0']['candidate_universe_preserved'], True)
        self.assertEqual(record['timeout_seconds_per_arm'], 7200)
        self.assertGreaterEqual(record['total_elapsed_seconds'], 0)
        with self.assertRaises(FileExistsError):
            C.run(self.root, self.root, self.root, self.root, out, self.root)
        self.assertFalse(list(self.root.glob('gspa-capacity-*')))

    def test_requested_arms_and_heap_are_validated_before_output(self):
        out = self.root/'capacity'
        for args in ({'arms': ['stage0', 'stage0']}, {'arms': []}, {'java_heap': '0g'},
                     {'timeout':0}, {'timeout':True}, {'timeout':1.5}):
            with self.assertRaises(ValueError):
                C.run(self.root, self.root, self.root, self.root, out, self.root, **args)
        self.assertFalse(out.exists())

    def test_capacity_deadline_and_total_time_are_recorded(self):
        out = self.root/'capacity'
        with patch.object(C, '_run', return_value={'status':'verified'}) as worker:
            C.run(self.root, self.root, self.root, self.root, out, self.root, timeout=14400)
        self.assertEqual(worker.call_args.args[-1], 14400)
        record = json.loads((out/'verification.json').read_text())
        self.assertEqual(record['timeout_seconds_per_arm'], 14400)
        self.assertGreaterEqual(record['total_elapsed_seconds'], 0)


if __name__ == '__main__':
    unittest.main()
