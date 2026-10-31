"""WP5 step 6 (design §8.6): COBRA replay bookkeeping. COBRApy's LP is mocked; skipped without cobra."""
import importlib.util
import os
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wp5_fixtures as F  # noqa: E402
import prepare_wp5_cases as P5  # noqa: E402
import run_wp5_recovery as R  # noqa: E402

HAVE_COBRA = importlib.util.find_spec('cobra') is not None and importlib.util.find_spec('libsbml') is not None
PROBLEM = 'W001-manuscript_regulatory'


class FakeModel:
    """Stands in for the COBRA model; counts copies (one per replayed cell is required)."""

    def __init__(self, counter):
        self.counter = counter

    def copy(self):
        self.counter.append(1)
        return FakeModel(self.counter)


def oracle_from_engine(drift=0.0):
    """An LP-free replacement for ``check_masked_assignment_cobra.simulate`` mirroring the fake engine."""
    def simulate(model, public, hypothesis, condition, cfg, knockout):
        assert isinstance(model, FakeModel)
        _, phenotypes = F.fake_phenotypes(hypothesis, condition['id'])
        record = dict(phenotypes['target_knockout' if knockout else 'wild_type'])
        record.pop('knocked_out_reactions')
        if drift and record['growth'] is not None and not knockout:
            record['growth'] += drift
        return record
    return simulate


@unittest.skipUnless(HAVE_COBRA, 'optional COBRApy independent solver')
class CobraReplay(unittest.TestCase):
    def setUp(self):
        import check_wp5_cobra as C
        self.C = C
        self.root = Path(tempfile.mkdtemp(prefix='wp5-cobra-', dir=os.environ.get('TMPDIR')))
        self.addCleanup(shutil.rmtree, self.root, True)
        self.world = F.runner_world(self.root)
        self.public, self.out = self.world['public'], self.root / 'runs'
        replies = ['{"experiment_id": "%s"}' % c for c in ('C2', 'C1', 'C0', 'C2', 'C1', 'C0', 'C3')]
        with F.engine_patch():
            R.make_protocol(self.public, self.out)
            for stage in ('exhaustive', 'deterministic', 'qwen'):
                chat = F.reply_chat(replies)
                R.execute(self.public, self.out, stage, chat_factory=lambda: chat)
        self.copies = []

    def compare(self, drift=0.0, **kwargs):
        with mock.patch.object(self.C.O, 'simulate', oracle_from_engine(drift)):
            return self.C.compare(self.public, self.out, loader=lambda path: FakeModel(self.copies), **kwargs)

    def test_every_visited_cell_is_replayed_from_its_own_model_copy_and_agrees(self):
        report = self.compare()
        runs = [r for r in report['runs']]
        self.assertEqual(len(runs), 2 * (3 + 10) + 2)            # every deterministic and Qwen run
        self.assertEqual(report['n_cell_comparisons'], sum(r['cells'] for r in runs))
        self.assertEqual(len(self.copies), report['n_cell_comparisons'])
        self.assertEqual(report['n_comparisons'], 2 * report['n_cell_comparisons'])
        self.assertEqual(report['n_agree'], report['n_comparisons'])
        self.assertEqual(report['n_cell_agree'], report['n_cell_comparisons'])
        self.assertEqual(report['masked_model_sha256'][PROBLEM], self.world['record']['masked_model_sha256'])

    def test_a_disagreeing_oracle_is_reported_not_hidden(self):
        report = self.compare(drift=0.5)
        self.assertLess(report['n_agree'], report['n_comparisons'])
        with mock.patch.object(self.C.O, 'simulate', oracle_from_engine(0.5)), \
                mock.patch.object(self.C, '_load_model', lambda path: FakeModel(self.copies)):
            code = self.C.main(['--public', str(self.public), '--runs', str(self.out), '--out', str(self.root / 'o.json')])
        self.assertEqual(code, 1)

    def test_rest_sample_is_seeded_and_excludes_visited_cells(self):
        a = self.compare(rest_sample=1.0)
        b = self.compare(rest_sample=1.0)
        self.assertEqual(a['cell_comparisons'], b['cell_comparisons'])
        exhaustive = [r for r in a['runs'] if '/exhaustive-' in r['run']]
        self.assertEqual(len(exhaustive), 1)
        self.assertEqual(exhaustive[0]['cells'], 0)              # every exhaustive cell was visited by some arm
        with self.assertRaises(P5.WP5Error):
            self.compare(rest_sample=1.5)
        with self.assertRaises(P5.WP5Error):
            self.compare(include_exhaustive=True, rest_sample=0.1)

    def test_shaheen_mode_replays_deterministic_arms_without_qwen_runs(self):
        for path in (self.out / 'runs').rglob('qwen-*.json'):
            path.unlink()
        with self.assertRaisesRegex(P5.WP5Error, 'missing run'):
            self.compare()
        report = self.compare(arm_stages=('deterministic',))
        self.assertEqual(len(report['runs']), 2 * (3 + 10))
        self.assertEqual(report['arm_stages'], ['deterministic'])

    def test_changed_masked_model_or_mismatched_run_is_refused(self):
        path = self.out / 'runs' / R.run_name(PROBLEM, 'evidence', 3, 0)
        run = P5.read_json(path)
        P5.write_json(path, dict(run, status='cell_mismatch'))
        with self.assertRaisesRegex(P5.WP5Error, 'invalid run provenance'):
            self.compare()
        P5.write_json(path, run)
        model = self.public / 'models/W001.xml'
        model.write_text(model.read_text() + '\n')
        with self.assertRaises(P5.WP5Error):
            self.compare()

    def test_no_private_file_is_opened(self):
        with F.refuse_paths(F.no_private) as opened:
            self.compare()
        self.assertFalse([p for p in opened if F.no_private(p)])


if __name__ == '__main__':
    unittest.main()
