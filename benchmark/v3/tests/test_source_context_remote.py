"""A successful process exit must not certify an incomplete experiment."""
import copy
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import run_source_context_remote as W


class RemoteCompletionTest(unittest.TestCase):
    def setUp(self):
        self.genomes = ['g' + str(i) for i in range(18)]
        exported = dict(status='verified', summaries={a: {**{c: {} for c in W.R.CONTEXTS}, 'comparison': {}}
                                                    for a in W.R.ARMS})
        frozen = dict(status='verified', counts={a: {c: {} for c in W.R.CONTEXTS} for a in W.R.ARMS})
        replayed = dict(status='verified', summaries={a: {c: {k: dict.fromkeys(self.genomes, {})
            for k in W.R.F.KINDS} for c in W.R.CONTEXTS} for a in W.R.ARMS})
        panel = dict(reference_pairs=39, proteins=37, test_genomes=['g0','g1','g2','g3'])
        results = dict(results={p: dict(full_ec_reference={a: copy.deepcopy(panel) for a in W.A.ARMS})
            for p in ('ec_calibration_selected','calibration_selected','fixed_0.5')})
        results['results']['native_emitted_calls'] = dict(full_ec_reference=dict(panel, threshold=0))
        self.args = [exported, frozen, replayed, dict.fromkeys(W.A.ARMS, {}), results, self.genomes]

    def test_complete_experiment(self):
        self.assertEqual(W.completion_counts(*self.args),
            dict(exports=20, fitted_models=40, replayed_genomes=720, analysis_arms=61))

    def test_empty_genome_cannot_disappear_from_replay(self):
        del self.args[2]['summaries']['clean']['SEQONLY']['logistic']['g17']
        with self.assertRaisesRegex(ValueError, 'incomplete replay genomes'):
            W.completion_counts(*self.args)

    def test_incomplete_analysis_and_changed_population_rejected(self):
        for change in ('missing_arm', 'missing_reference', 'native_cut'):
            args = copy.deepcopy(self.args)
            panels = args[4]['results']['ec_calibration_selected']['full_ec_reference']
            if change == 'missing_arm': del panels['clean_native']
            elif change == 'missing_reference': panels['clean_native']['reference_pairs'] = 38
            else: args[4]['results']['native_emitted_calls']['full_ec_reference']['threshold'] = .5
            with self.subTest(change=change), self.assertRaises(ValueError):
                W.completion_counts(*args)

    def test_prerequisite_must_be_bound_and_unchanged(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); prerequisite = root / 'audit.json'; manifest = root / 'packages.sha256'
            prerequisite.write_text('{}'); manifest.write_text('')
            with self.assertRaisesRegex(ValueError, 'unbound prerequisite'):
                W.bindings(manifest, [prerequisite])
            digest = W.R.S.sha(prerequisite)
            manifest.write_text(digest + '  ' + str(prerequisite) + '\n')
            self.assertEqual(W.bindings(manifest, [prerequisite]), {str(prerequisite): digest})
            prerequisite.write_text('{"changed":true}')
            with self.assertRaisesRegex(ValueError, 'unbound prerequisite'):
                W.bindings(manifest, [prerequisite])

    def test_missing_probe_sources_do_not_pass_vacuously(self):
        with self.assertRaisesRegex(ValueError, 'missing probed pipeline sources'):
            W.verify_probe_sources(dict(inputs={},source_root='/scratch/probe/source'))
        prefix = '/scratch/10389/source-probe/source/'
        names = ['benchmark/v3/' + name for name in ('run_source_arm_context_comparison.py',
                 'verified_batch.py', 'cohort_sequence_inputs.py')]
        inputs = {prefix + name: W.R.S.sha(W.ROOT / name) for name in names}
        W.verify_probe_sources(dict(inputs=inputs,source_root=prefix.rstrip('/')))
        inputs[prefix + names[0]] = '0' * 64
        with self.assertRaisesRegex(ValueError, 'probed source changed'):
            W.verify_probe_sources(dict(inputs=inputs,source_root=prefix.rstrip('/')))


class RemoteReadabilityTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR', str(Path.home() / 'tmp')))
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def test_readable_files_pass_without_retries(self):
        first = self.root / 'first.bin'
        second = self.root / 'second.bin'
        first.write_bytes(b'first')
        second.write_bytes(b'second')
        paths = {str(path): W.R.S.sha(path) for path in (first, second)}
        with patch.object(W.time, 'sleep') as sleep:
            self.assertEqual(W.check_readable(paths), [])
        sleep.assert_not_called()

    def test_transient_read_failure_is_retried_and_recorded(self):
        path = self.root / 'flaky.bin'
        path.write_bytes(b'restored')
        registered = W.R.S.sha(path)
        original_open = Path.open
        failed = False

        def flaky_open(candidate, *args, **kwargs):
            nonlocal failed
            if candidate == path and not failed:
                failed = True
                raise OSError(5, 'Input/output error')
            return original_open(candidate, *args, **kwargs)

        with patch.object(Path, 'open', flaky_open), patch.object(W.time, 'sleep') as sleep:
            self.assertEqual(W.check_readable({str(path): registered}), [str(path)])
        sleep.assert_called_once_with(15)

    def test_persistent_read_failure_names_path_and_errno(self):
        path = self.root / 'unreadable.bin'
        path.write_bytes(b'restored')
        registered = W.R.S.sha(path)
        with patch.object(Path, 'open', side_effect=OSError(5, 'Input/output error')), \
             patch.object(W.time, 'sleep') as sleep:
            with self.assertRaises(RuntimeError) as caught:
                W.check_readable({str(path): registered})
        message = str(caught.exception)
        self.assertIn('unreadable restored input after 3 attempts', message)
        self.assertIn(str(path), message)
        self.assertIn('[Errno 5] Input/output error', message)
        self.assertEqual(sleep.call_count, 2)

    def test_changed_restored_content_names_path(self):
        path = self.root / 'changed.bin'
        path.write_bytes(b'truncated')
        with self.assertRaises(RuntimeError) as caught:
            W.check_readable({str(path): '0' * 64})
        self.assertEqual(str(caught.exception),
                         f'restored input differs from registration: {path}')


class RemoteRestoreTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR', str(Path.home() / 'tmp')))
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.out = self.root / 'out'
        self.logs = self.root / 'logs'
        self.out.mkdir()
        self.logs.mkdir()
        self.training = self.root / 'sealed/training'
        self.replay = self.root / 'sealed/replay'
        self.training.mkdir(parents=True)
        self.replay.mkdir(parents=True)
        self.args = SimpleNamespace(out=self.out, logs=self.logs, resume_training=self.training,
                                    resume_replay=self.replay, restore_wait_seconds=7200)
        self.learned = self.out / 'training'
        self.restored_replay = self.out / 'replay'

    def write_marker(self):
        marker = dict(sources=dict(training=str(self.training.resolve()),
                                   replay=str(self.replay.resolve())))
        (self.out / 'restored.json').write_text(json.dumps(marker))

    def test_restore_marker_skips_rsync_and_readability_but_reverifies_training(self):
        self.write_marker()
        self.learned.mkdir()
        self.restored_replay.mkdir()
        with patch.object(W, 'verify_frozen_training') as verify, \
             patch.object(W, 'check_readable') as readable, \
             patch.object(W.subprocess, 'run') as run:
            self.assertIsNone(W.restore_resume_tree(self.args, self.learned, self.restored_replay))
        verify.assert_called_once_with(self.learned)
        readable.assert_not_called()
        run.assert_not_called()

    def test_concurrent_restorer_waits_for_marker_and_never_removes_lock(self):
        lock = self.out / 'restore.lock'
        descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        os.close(descriptor)
        self.learned.mkdir()
        self.restored_replay.mkdir()
        with patch.object(W.time, 'sleep', side_effect=lambda seconds: self.write_marker()) as sleep, \
             patch.object(W, 'verify_frozen_training') as verify, \
             patch.object(W, 'check_readable') as readable, \
             patch.object(W.subprocess, 'run') as run:
            self.assertIsNone(W.restore_resume_tree(self.args, self.learned, self.restored_replay))
        sleep.assert_called_once_with(30)
        verify.assert_called_once_with(self.learned)
        readable.assert_not_called()
        run.assert_not_called()
        self.assertTrue(lock.is_file())

    def test_stale_restore_lock_fails_without_being_deleted(self):
        lock = self.out / 'restore.lock'
        lock.touch()
        os.utime(lock, (0, 0))
        self.args.restore_wait_seconds = 10
        with patch.object(W.time, 'sleep') as sleep, \
             patch.object(W, 'verify_frozen_training') as verify, \
             self.assertRaisesRegex(RuntimeError, 'stale restore lock'):
            W.restore_resume_tree(self.args, self.learned, self.restored_replay)
        sleep.assert_not_called()
        verify.assert_not_called()
        self.assertTrue(lock.is_file())


class RemoteResumeTest(unittest.TestCase):
    changed_sources = ('benchmark/v3/check_source_arm_context_models.py',
                       'benchmark/v3/run_source_context_remote.py')
    required_sources = tuple('benchmark/v3/' + name for name in
        ('run_source_arm_context_comparison.py', 'verified_batch.py', 'cohort_sequence_inputs.py'))

    def setUp(self):
        RemoteCompletionTest.setUp(self)
        self.tmp = tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR', str(Path.home() / 'tmp')))
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.sealed_training = self.root / 'jobs/10406/checkpoint/training'
        self.sealed_replay = self.root / 'jobs/10406/results/replay'
        self.sealed_training.mkdir(parents=True)
        self.sealed_replay.mkdir(parents=True)
        self.out = self.root / 'resumed'
        self.logs = self.root / 'logs'

    def arguments(self):
        values = {}
        for name in ('jar', 'clean-root', 'clean-audit', 'bindings'):
            path = self.root / name
            path.write_text('{}')
            values[name] = path
        results_root = self.root / 'probe-results'
        probe = dict(schema='gspa.source-context-probe/2', status='verified',
                     summaries=dict.fromkeys(('diamond_homology', 'clean', 'clean_diamond_homology'), {}),
                     historical_model_audit=dict.fromkeys(W.R.F.KINDS, {}),
                     inputs={str(values['jar']): W.R.S.sha(values['jar'])},
                     results_root=str(results_root), outputs={str(results_root / 'result.json'): '0' * 64})
        source_root = '/scratch/10403/source'
        probe['source_root'] = source_root
        self.source_drift = {}
        for name in self.required_sources + self.changed_sources:
            current_sha = W.R.S.sha(W.ROOT / name)
            probe_sha = '0' * 64 if name in self.changed_sources else current_sha
            probe['inputs'][source_root + '/' + name] = probe_sha
            if name in self.changed_sources:
                self.source_drift[name] = dict(probe_sha256=probe_sha, current_sha256=current_sha)
        values['probe'] = self.root / 'probe.json'
        values['probe'].write_text(json.dumps(probe))
        args = ['remote.py', 'resume']
        for name, path in values.items():
            args += ['--' + name, str(path)]
        args += ['--out', str(self.out), '--logs', str(self.logs),
                 '--resume-training', str(self.sealed_training),
                 '--resume-replay', str(self.sealed_replay)]
        for name in self.changed_sources:
            args += ['--resume-changed-source', name]
        return args

    def prepare_sealed_training(self):
        frozen_outputs = {}
        for arm in W.R.ARMS:
            for context in W.R.CONTEXTS:
                for kind in W.R.F.KINDS:
                    relative = Path(arm) / context / 'models' / ('model_' + kind + '.json')
                    model = self.sealed_training / relative
                    model.parent.mkdir(parents=True, exist_ok=True)
                    model.write_text(json.dumps(dict(arm=arm, context=context, kind=kind)))
                    frozen_outputs[str(self.out / 'training' / relative)] = W.R.S.sha(model)
        frozen = dict(self.args[1], outputs=frozen_outputs)
        (self.sealed_training / 'training.json').write_text(json.dumps(frozen))
        (self.sealed_training / 'export.json').write_text(json.dumps(self.args[0]))
        for name in ('protocol.json', 'registration.json', 'evaluation.json'):
            content = '{"inputs": {}}' if name == 'registration.json' else '{}'
            (self.sealed_training / name).write_text(content)
        cache_training = self.root / 'cache/observed-gap-models-v1'
        cache_training.mkdir(parents=True)
        (cache_training / 'registration.json').write_text('{"inputs": {}}')

    def run_main(self, execute, extra=()):
        with patch.object(W.sys, 'argv', self.arguments() + list(extra)), \
             patch.object(W, 'CACHE', self.root / 'cache'), \
             patch.object(W, 'bindings', return_value={}), \
             patch.object(W.R, 'base_inputs', return_value=({'recipes': {}}, {'genomes': self.genomes}, {})), \
             patch.object(W.R.C, 'supervision_inputs', return_value={}), \
             patch.object(W.R, 'supervision_agreement'), \
             patch.object(W.R.C, 'sources', return_value=[]), \
             patch.object(W.R.F.V, 'streaming_parser_sources', return_value=[]), \
             patch.dict(W.sys.modules, {'ijson': SimpleNamespace(__version__='3.4.0.post0', backend='python')}), \
             patch.object(W.R.F, 'hashes', return_value={}), \
             patch.object(W.R.F, 'check_hashes'), \
             patch.object(W.R, 'verify_clean'), \
             patch.object(W.R, 'verify'), \
             patch.object(W.subprocess, 'run', side_effect=execute):
            W.main()

    def test_resume_runs_only_replay_and_analysis_and_records_parent(self):
        self.prepare_sealed_training()
        replay = self.out / 'replay'
        analysis = self.out / 'analysis'
        stage_calls = []

        def execute(command, **kwargs):
            command = [str(x) for x in command]
            if command[0] == 'rsync':
                source = Path(command[-2].rstrip('/'))
                destination = Path(command[-1].rstrip('/'))
                if source in (self.sealed_training, self.sealed_replay):
                    shutil.copytree(source, destination)
                return SimpleNamespace(returncode=0)
            stage_calls.append(command)
            script = Path(command[1])
            if script == Path(W.M.__file__):
                W.R.F.save(replay / 'verification.json', self.args[2])
            elif script == Path(W.A.__file__) and command[2] == 'register':
                analysis.mkdir(parents=True)
                for name in ('protocol.json', 'registration.json'):
                    W.R.F.save(analysis / name, {})
            elif script == Path(W.A.__file__) and command[2] == 'run':
                W.R.F.save(analysis / 'thresholds.json', self.args[3])
                W.R.F.save(analysis / 'results.json', self.args[4])
                W.R.F.save(analysis / 'verification.json', dict(status='verified'))
            return SimpleNamespace(returncode=0)

        self.run_main(execute)

        self.assertEqual([Path(call[1]).name for call in stage_calls],
                         ['check_source_arm_context_models.py', 'analyze_source_arm_context_models.py',
                          'analyze_source_arm_context_models.py'])
        self.assertIn('--resume', stage_calls[0])
        self.assertEqual([call[2] for call in stage_calls[1:]], ['register', 'run'])
        state = json.loads((self.logs / 'workflow.json').read_text())
        self.assertEqual([s['name'] for s in state['stages']],
                         ['production-replay', 'analysis-registration', 'analysis'])
        self.assertEqual(state['parent_job'], '10406')
        self.assertEqual(state['resume'], dict(training=str(self.sealed_training), replay=str(self.sealed_replay),
                                               source_drift=self.source_drift,
                                               readability=dict(files=40, retried=[])))
        self.assertTrue((self.logs / 'restore-training.log').is_file())
        self.assertTrue((self.logs / 'restore-replay.log').is_file())
        self.assertTrue((self.out / 'restored.json').is_file())
        self.assertFalse((self.out / 'restore.lock').exists())
        receipt = json.loads((self.out / 'final-receipt.json').read_text())
        self.assertEqual(receipt['mode'], 'resume')
        self.assertEqual(receipt['parent_job'], '10406')
        self.assertEqual(receipt['counts']['replayed_genomes'], 720)
        self.assertEqual(receipt['source_drift'], self.source_drift)

    def test_replay_arm_resume_has_no_analysis_checkpoint_or_final_verification(self):
        self.prepare_sealed_training()
        stage_calls, rsync_calls = [], []

        def execute(command, **kwargs):
            command = [str(x) for x in command]
            if command[0] == 'rsync':
                rsync_calls.append(command)
                source = Path(command[-2].rstrip('/'))
                destination = Path(command[-1].rstrip('/'))
                if source in (self.sealed_training, self.sealed_replay):
                    shutil.copytree(source, destination)
                return SimpleNamespace(returncode=0)
            stage_calls.append(command)
            receipt = Path(command[command.index('--partial-receipt') + 1])
            W.R.F.save(receipt, dict(schema='gspa.source-context-production-partial/1', status='verified'))
            return SimpleNamespace(returncode=0)

        arms = ['clean/UNIVERSAL/logistic', 'diamond_homology/SEQONLY/gbdt']
        extra = ['--replay-arm', arms[0], '--replay-arm', arms[1]]
        self.run_main(execute, extra)
        self.assertEqual(len(stage_calls), 1)
        command = stage_calls[0]
        self.assertEqual(Path(command[1]), Path(W.M.__file__))
        self.assertIn('--resume', command)
        self.assertEqual([command[i + 1] for i, value in enumerate(command) if value == '--arm'], arms)
        self.assertEqual(command[command.index('--partial-receipt') + 1],
                         str(self.logs / 'replay-partial.json'))
        self.assertEqual(len(rsync_calls), 2)
        self.assertFalse((self.out / 'analysis').exists())
        self.assertFalse((self.out / 'replay/verification.json').exists())
        self.assertFalse((self.out / 'final-receipt.json').exists())
        self.assertFalse((self.root / 'checkpoint').exists())
        state = json.loads((self.logs / 'workflow.json').read_text())
        self.assertEqual(state['mode'], 'resume-arms')
        self.assertEqual(state['status'], 'verified')
        self.assertEqual([stage['name'] for stage in state['stages']], ['production-replay'])
        restored = json.loads((self.out / 'restored.json').read_text())
        self.assertEqual(restored['sources'], dict(training=str(self.sealed_training.resolve()),
                                                   replay=str(self.sealed_replay.resolve())))
        self.assertEqual(restored['file_counts'], dict(training=45, replay=0))
        self.assertIn('restored_at', restored)
        self.assertEqual(restored['training_receipts'], {})
        self.assertFalse((self.out / 'restore.lock').exists())

    def test_resume_rejects_undeclared_changed_probe_source(self):
        source_root = '/scratch/10403/source'
        inputs = {source_root + '/' + name: W.R.S.sha(W.ROOT / name)
                  for name in self.required_sources + self.changed_sources}
        changed = self.required_sources[0]
        inputs[source_root + '/' + changed] = '0' * 64
        with self.assertRaisesRegex(ValueError, changed):
            W.verify_probe_sources(dict(inputs=inputs, source_root=source_root), resume_changed_sources=[])

    def test_resume_rejects_declaration_absent_from_probe(self):
        source_root = '/scratch/10403/source'
        inputs = {source_root + '/' + name: W.R.S.sha(W.ROOT / name)
                  for name in self.required_sources + self.changed_sources}
        missing = 'benchmark/v3/not_probed.py'
        with self.assertRaisesRegex(ValueError, missing):
            W.verify_probe_sources(dict(inputs=inputs, source_root=source_root),
                                   resume_changed_sources=[missing])

    def test_changed_source_declaration_rejected_in_full_mode(self):
        args = self.arguments()
        args[1] = 'full'
        with patch.object(W.sys, 'argv', args), \
             self.assertRaisesRegex(ValueError, '--resume-changed-source is only valid in resume mode'):
            W.main()

    def test_resume_refuses_an_existing_training_destination(self):
        (self.out / 'training').mkdir(parents=True)
        with patch.object(W.sys, 'argv', self.arguments()), \
             self.assertRaisesRegex(ValueError, 'resume destinations already exist'):
            W.main()
        self.assertFalse(self.logs.exists())


if __name__ == '__main__': unittest.main()
