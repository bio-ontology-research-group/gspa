"""Regression tests for reuse of actual process outputs (stub engine, no biological claims)."""
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import verified_batch as V
import run_training_cohort_ablations as A

STUB = '''import csv,json,sys
from pathlib import Path
args=sys.argv[1:]
rows=list(csv.DictReader(open(args[args.index('--batch-manifest')+1]), delimiter='\\t'))
for r in rows:
 p=Path(r['out_path'])
 p.write_text('protein_id\\ttype\\tfunction_id\\tgo_aspect\\tposterior_prob\\np\\tEC\\t1.1.1.1\\t\\t0.3\\n')
 stages=[{'name':s,'status':'executed'} for s in ['reference_preparation','evidence_integration']]
 if '--rlc-suggester' in args: stages.append({'name':'reaction_context','status':'skipped','reason':'no gaps'})
 if '--broken-ledger' in args: stages[-1]['status']='failed'
 Path(str(p)+'.stages.json').write_text(json.dumps({'schema_version':2,'genome_id':r['genome_id'],'stages':stages}))
if '--fail-after-write' in args: sys.exit(2)
'''


class VerifiedBatchTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.out = self.root / 'out'
        self.out.mkdir()
        self.engine = self.root / 'engine.py'
        self.engine.write_text(STUB)
        self.claims = self.root / 'claims.jsonl'
        self.claims.write_text('{}\n')
        self.layout = self.root / 'layout.tsv'
        self.layout.write_text('layout\n')
        self.model = self.root / 'theta.json'
        self.model.write_text('{}')
        self.manifest = self.root / 'batch.tsv'
        self.prediction = self.out / 'g.integrated.tsv'
        self.manifest.write_text('genome_id\tclaims_path\tout_path\tgenome_layout_path\n'
                                 f'g\t{self.claims}\t{self.prediction}\t{self.layout}\n')
        self.cmd = [sys.executable, str(self.engine), '--batch-manifest', str(self.manifest),
                    '--theta', str(self.model)]

    def run_batch(self, extra=()):
        return V.run_batch(self.cmd + list(extra), self.manifest, self.out, sources=[self.engine])

    def spec(self, extra=()):
        return V.describe_run(self.cmd + list(extra), self.manifest, [self.engine])

    def test_matching_success_is_reused_without_reexecution(self):
        first = self.run_batch()
        time = self.prediction.stat().st_mtime_ns
        self.assertEqual(self.run_batch(), first)
        self.assertEqual(self.prediction.stat().st_mtime_ns, time)
        self.assertFalse((self.out / 'history').exists())

    def test_each_inference_dependency_invalidates_reuse(self):
        self.run_batch()
        for p in [self.claims, self.layout, self.model, self.engine]:
            with self.subTest(path=p.name):
                old = p.read_text()
                p.write_text(old + '\n')
                self.assertFalse(V.reusable(self.out, self.spec()))
                p.write_text(old)
                self.assertTrue(V.reusable(self.out, self.spec()))
        self.assertFalse(V.reusable(self.out, self.spec(['--changed-option'])))

    def test_changed_outputs_and_ledger_invalidate_reuse(self):
        self.run_batch()
        for p in [self.prediction, Path(str(self.prediction) + '.stages.json')]:
            old = p.read_text()
            p.write_text(old + '\n')
            self.assertFalse(V.reusable(self.out, self.spec()))
            p.write_text(old)

    def test_existing_unverified_outputs_are_preserved_and_rerun(self):
        self.prediction.write_text('unverified old predictions')
        self.run_batch()
        archived = list((self.out / 'history').glob('*/g.integrated.tsv'))
        self.assertEqual(len(archived), 1)
        self.assertEqual(archived[0].read_text(), 'unverified old predictions')
        self.assertTrue(V.reusable(self.out, self.spec()))

    def test_failure_after_writing_cannot_be_reused(self):
        self.run_batch()
        with self.assertRaisesRegex(RuntimeError, 'batch failed'):
            self.run_batch(['--fail-after-write'])
        self.assertFalse((self.out / V.RECEIPT).exists())
        self.assertFalse(V.reusable(self.out, self.spec()))

    def test_zero_exit_with_failed_stage_is_not_a_success(self):
        with self.assertRaisesRegex(ValueError, 'did not execute'):
            self.run_batch(['--broken-ledger'])
        self.assertFalse((self.out / V.RECEIPT).exists())

    def test_requested_ablation_stage_must_have_executed(self):
        with self.assertRaisesRegex(ValueError, 'reaction_context did not execute'):
            self.run_batch(['--rlc-suggester'])

    def test_missing_input_is_not_silently_ignored(self):
        self.run_batch()
        self.layout.unlink()
        with self.assertRaises(FileNotFoundError):
            self.run_batch()

    def test_equals_options_and_multiple_module_files_are_hashed(self):
        self.cmd = self.cmd[:-2] + [f'--theta={self.model}', f'--modules={self.layout},{self.claims}']
        self.run_batch()
        for path in (self.model, self.layout, self.claims):
            with self.subTest(path=path.name):
                original = path.read_text()
                path.write_text(original + '\n')
                self.assertFalse(V.reusable(self.out, self.spec()))
                path.write_text(original)

    def test_directory_dependency_changes_invalidate_reuse(self):
        data = self.root / 'data'
        data.mkdir()
        ref = data / 'reference.dat'
        ref.write_text('first')
        self.cmd += ['--fba-data-dir', str(data)]
        self.run_batch()
        ref.write_text('second')
        self.assertFalse(V.reusable(self.out, self.spec()))
        ref.unlink()
        self.assertFalse(V.reusable(self.out, self.spec()))

    def test_production_jar_follows_source_version_and_explicit_override(self):
        (self.root / 'build.gradle.kts').write_text('version = "3.0.0"\n')
        with patch.dict('os.environ', {}, clear=True):
            self.assertEqual(A.production_jar(self.root).name, 'gspa-3.0.0.jar')
        with patch.dict('os.environ', {'GSPA_JAR': str(self.root / 'chosen.jar')}):
            self.assertEqual(A.production_jar(self.root), self.root / 'chosen.jar')

    def test_mutation_during_execution_is_not_certified(self):
        self.engine.write_text(STUB + f'\nPath({str(self.claims)!r}).write_text("changed during inference")\n')
        with self.assertRaisesRegex(RuntimeError, 'inputs changed during execution'):
            self.run_batch()
        self.assertFalse((self.out / V.RECEIPT).exists())

    def test_outputs_outside_run_directory_are_rejected_before_mutation(self):
        self.manifest.write_text(self.manifest.read_text().replace(str(self.prediction), str(self.root / 'other.tsv')))
        with self.assertRaisesRegex(ValueError, 'inside out_dir'):
            self.run_batch()

    def test_timeout_stops_wrapped_inference_and_signal_ignoring_child(self):
        pidfile = self.root/'child.pid'
        child = ('import os,signal,sys,time; from pathlib import Path; '
                 'signal.signal(signal.SIGTERM,signal.SIG_IGN); '
                 'Path(sys.argv[1]).write_text(str(os.getpid())); time.sleep(60)')
        self.engine.write_text('import subprocess,sys\n'
                              f'subprocess.run([sys.executable,"-c",{child!r},sys.argv[1]])\n')
        wrapper = self.root/'wrapper.py'
        wrapper.write_text('import subprocess,sys\nsubprocess.run(sys.argv[1:])\n')
        command = [sys.executable, str(wrapper), sys.executable, str(self.engine), str(pidfile),
                   '--batch-manifest',str(self.manifest)]
        pid = None
        try:
            with self.assertRaises(subprocess.TimeoutExpired):
                V.run_batch(command, self.manifest, self.out,
                            sources=[wrapper, self.engine, Path(sys.executable)], timeout=1)
            pid = int(pidfile.read_text())
            def running():
                status = Path(f'/proc/{pid}/status')
                try: lines = status.read_text().splitlines()
                except FileNotFoundError: return False
                return not any(line.startswith('State:') and 'Z' in line.split()[1]
                               for line in lines)
            deadline = time.monotonic()+2
            while running() and time.monotonic() < deadline:
                time.sleep(.02)
            self.assertFalse(running(), 'timed-out inference left its child running')
            self.assertFalse((self.out/V.RECEIPT).exists())
        finally:
            # Also cleans up the deliberately failing pre-fix regression.
            if pidfile.exists():
                try: os.kill(int(pidfile.read_text()), signal.SIGKILL)
                except ProcessLookupError: pass


if __name__ == '__main__':
    unittest.main()
