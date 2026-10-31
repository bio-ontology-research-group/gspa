"""Sharded COBRA replays merge into exactly the report of an unsharded replay (LP mocked)."""
import json
import os
import sys
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import test_wp5_cobra as T  # noqa: E402
import prepare_wp5_cases as P5  # noqa: E402


@unittest.skipUnless(T.HAVE_COBRA, 'optional COBRApy independent solver')
class CobraMerge(T.CobraReplay):
    def shards(self, n):
        import merge_wp5_cobra as M
        self.M = M
        d = self.root / ('shards-%d' % n)
        for i in range(n):
            P5.write_json(d / ('shard-%d-of-%d.json' % (i, n)),
                          self.compare(arm_stages=('deterministic',), shard='%d/%d' % (i, n)))
        return d

    def merged(self, n):
        d = self.shards(n)
        out = self.root / ('merged-%d.json' % n)
        code = self.M.main(['--runs', str(self.out), '--shard-dir', str(d), '--shards', str(n), '--out', str(out)])
        return code, out

    def test_merged_shards_are_byte_identical_to_an_unsharded_replay(self):
        whole = self.root / 'whole.json'
        P5.write_json(whole, self.compare(arm_stages=('deterministic',)))
        for n in (1, 2, 3, 7, 26):
            code, out = self.merged(n)
            self.assertEqual(code, 0)
            self.assertEqual(out.read_bytes(), whole.read_bytes(), n)
            receipt = json.loads(Path(str(out) + '.merge.json').read_text())
            self.assertEqual(receipt['shards'], n)
            self.assertEqual(len(receipt['shard_sha256']), n)

    def test_missing_swapped_or_tampered_shards_are_refused(self):
        d = self.shards(3)
        (d / 'shard-2-of-3.json').unlink()
        with self.assertRaisesRegex(self.M.MergeError, 'missing shard'):
            self.M.main(['--runs', str(self.out), '--shard-dir', str(d), '--shards', '3', '--out', str(self.root / 'x.json')])
        d = self.shards(3)
        a, b = d / 'shard-0-of-3.json', d / 'shard-1-of-3.json'
        ta, tb = a.read_bytes(), b.read_bytes()
        a.write_bytes(tb)
        b.write_bytes(ta)
        with self.assertRaisesRegex(self.M.MergeError, 'shard 0 holds'):
            self.M.main(['--runs', str(self.out), '--shard-dir', str(d), '--shards', '3', '--out', str(self.root / 'x.json')])
        a.write_bytes(ta)
        b.write_bytes(tb)
        rep = json.loads(a.read_text())
        rep['comparisons'].pop()
        rep['n_comparisons'] -= 1
        P5.write_json(a, rep)
        with self.assertRaisesRegex(self.M.MergeError, 'endpoint comparisons|do not add up'):
            self.M.main(['--runs', str(self.out), '--shard-dir', str(d), '--shards', '3', '--out', str(self.root / 'x.json')])

    def test_a_disagreeing_shard_makes_the_merge_fail(self):
        d = self.shards(2)
        with mock.patch.object(self.C.O, 'simulate', T.oracle_from_engine(0.5)):
            rep = self.C.compare(self.public, self.out, loader=lambda path: T.FakeModel(self.copies),
                                 arm_stages=('deterministic',), shard='1/2')
        P5.write_json(d / 'shard-1-of-2.json', rep)
        code = self.M.main(['--runs', str(self.out), '--shard-dir', str(d), '--shards', '2', '--out', str(self.root / 'y.json')])
        self.assertEqual(code, 1)


# Do not re-run the inherited CobraReplay tests under this module.
for _name in [n for n in dir(T.CobraReplay) if n.startswith('test_')]:
    setattr(CobraMerge, _name, None)

if __name__ == '__main__':
    unittest.main()
