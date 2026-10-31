"""Execute the frozen source/context comparison in its restored runtime.

Called by source_context_comparison.sbatch. Preflight verifies the parent data
and numerical runtime without fitting. Full execution also requires the complete
native CLEAN audit and successful real production probe. All stages fail closed;
the scheduler wrapper preserves both failed and completed attempts.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import run_source_arm_context_comparison as R
import check_source_arm_context_models as M
import analyze_source_arm_context_models as A

CACHE = Path('/home/leechuck/.cache/gspa-v3-data')
ROOT = Path(__file__).resolve().parents[2]


def completion_counts(exported, frozen, replayed, thresholds, results, genomes):
    """Reject a zero-exit partial workflow, including missing empty-genome runs."""
    for receipt in (exported, frozen, replayed):
        R.require(receipt.get('status') == 'verified', 'unverified stage receipt')
    for tree in (exported['summaries'], frozen['counts'], replayed['summaries']):
        R.require(set(tree) == set(R.ARMS), 'incomplete source arms')
        for arm in R.ARMS:
            required = set(R.CONTEXTS) | ({'comparison'} if tree is exported['summaries'] else set())
            R.require(set(tree[arm]) == required, 'incomplete context arms')
    for contexts in replayed['summaries'].values():
        for kinds in contexts.values():
            R.require(set(kinds) == set(R.F.KINDS), 'incomplete model kinds')
            for records in kinds.values():
                R.require(set(records) == set(genomes), 'incomplete replay genomes')
    R.require(set(thresholds) == set(A.ARMS), 'incomplete frozen thresholds')
    for policy in ('ec_calibration_selected', 'calibration_selected', 'fixed_0.5'):
        panels = results['results'][policy]['full_ec_reference']
        R.require(set(panels) == set(A.ARMS), 'incomplete EC analysis arms')
        for panel in panels.values():
            R.require(panel['reference_pairs'] == 39 and panel['proteins'] == 37 and
                      len(panel['test_genomes']) == 4, 'full EC reference population changed')
    native = results['results']['native_emitted_calls']['full_ec_reference']
    R.require(native['reference_pairs'] == 39 and native['proteins'] == 37 and
              native['threshold'] == 0, 'native emitted-call endpoint changed')
    return dict(exports=20, fitted_models=40, replayed_genomes=40 * len(genomes), analysis_arms=61)


def bindings(path, required):
    bound = {}
    for line in path.read_text().splitlines():
        digest, name = line.split('  ', 1)
        R.require(name not in bound, 'duplicate bundle input')
        bound[name] = digest
    for item in required:
        R.require(bound.get(str(item)) == R.S.sha(item), 'unbound prerequisite: ' + str(item))
    return {str(item): bound[str(item)] for item in required}


def resume_parent_job(path):
    parts = Path(path).parts
    return next((parts[i + 1] for i, part in enumerate(parts[:-1]) if part == 'jobs'), None)


def check_readable(paths, attempts=3, delay=15):
    """Read and digest restored files, retrying transient storage errors."""
    retried = []
    for path, registered_digest in paths.items():
        path = Path(path)
        for attempt in range(attempts):
            try:
                digest = hashlib.sha256()
                with path.open('rb') as stream:
                    while chunk := stream.read(1024 * 1024):
                        digest.update(chunk)
                break
            except OSError as error:
                if attempt == attempts - 1:
                    raise RuntimeError(
                        f'unreadable restored input after {attempts} attempts: {path}: {error}') from error
                if attempt == 0:
                    retried.append(str(path))
                time.sleep(delay)
        if digest.hexdigest() != registered_digest:
            raise RuntimeError(f'restored input differs from registration: {path}')
    return retried


def verify_frozen_training(training):
    R.verify(training)
    frozen = json.loads((training / 'training.json').read_text())
    R.require(frozen.get('status') == 'verified', 'unverified frozen training')
    for arm in R.ARMS:
        for context_name in R.CONTEXTS:
            for kind in R.F.KINDS:
                model = training / arm / context_name / 'models' / ('model_' + kind + '.json')
                R.require(frozen['outputs'].get(str(model)) == R.S.sha(model),
                          'missing or changed frozen model')
    return frozen


def check_restored_sources(args):
    restored = json.loads((args.out / 'restored.json').read_text())
    expected = dict(training=str(args.resume_training.resolve()), replay=str(args.resume_replay.resolve()))
    R.require(restored.get('sources') == expected, 'restored sources differ from resume sources')
    return restored


def restore_resume_tree(args, learned, replay):
    """Restore a shared tree once; concurrent callers wait for its marker."""
    marker = args.out / 'restored.json'
    if marker.is_file():
        check_restored_sources(args)
        verify_frozen_training(learned)
        return None
    lock = args.out / 'restore.lock'
    try:
        descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError:
        deadline = time.monotonic() + args.restore_wait_seconds
        while not marker.is_file():
            try:
                modified = lock.stat().st_mtime
            except FileNotFoundError:
                if marker.is_file():
                    break
                raise RuntimeError('restore lock disappeared before restored.json was written')
            if time.time() - modified > args.restore_wait_seconds:
                raise RuntimeError('stale restore lock: ' + str(lock))
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                raise RuntimeError('timed out waiting for restored.json: ' + str(marker))
            time.sleep(min(30, remaining))
        check_restored_sources(args)
        verify_frozen_training(learned)
        return None
    try:
        if marker.is_file():
            check_restored_sources(args)
            verify_frozen_training(learned)
            return None
        R.require(not learned.exists() and not replay.exists(), 'resume destinations already exist')
        restores = ((args.resume_training, learned, 'restore-training.log'),
                    (args.resume_replay, replay, 'restore-replay.log'))
        for source, destination, log_name in restores:
            with (args.logs / log_name).open('x') as log:
                result = subprocess.run(['rsync', '-a', str(source) + '/', str(destination) + '/'],
                                        stdout=log, stderr=subprocess.STDOUT)
            if result.returncode:
                raise RuntimeError(f'restore rsync failed for {source} (exit {result.returncode})')
        registration = json.loads((learned / 'registration.json').read_text())
        frozen = json.loads((learned / 'training.json').read_text())
        restored_inputs = dict(registration['inputs'])
        restored_inputs.update(frozen['outputs'])
        retried = check_readable(restored_inputs)
        verify_frozen_training(learned)
        receipt_paths = [learned / name for name in
                         ('protocol.json', 'registration.json', 'export.json', 'training.json', 'evaluation.json')]
        restored = dict(schema='gspa.source-context-restored/1', status='verified',
                        sources=dict(training=str(args.resume_training.resolve()),
                                     replay=str(args.resume_replay.resolve())),
                        file_counts={tree.name: sum(p.is_file() for p in tree.rglob('*'))
                                     for tree in (learned, replay)},
                        restored_at=datetime.now(timezone.utc).isoformat(),
                        training_receipts=R.F.hashes(receipt_paths))
        temporary = marker.with_suffix('.partial')
        temporary.write_text(json.dumps(restored, indent=2, sort_keys=True) + '\n')
        temporary.replace(marker)
        return dict(files=len(restored_inputs), retried=retried)
    finally:
        os.close(descriptor)
        try:
            lock.unlink()
        except FileNotFoundError:
            pass


def verify_probe_sources(probe, root=ROOT, resume_changed_sources=None):
    source_root = probe.get('source_root')
    R.require(isinstance(source_root,str) and Path(source_root).is_absolute() and source_root != '/',
              'missing probe source root')
    prefix = source_root.rstrip('/') + '/'
    found = {old[len(prefix):]: digest for old, digest in probe['inputs'].items() if old.startswith(prefix)}
    R.require({'benchmark/v3/' + name for name in ('run_source_arm_context_comparison.py',
        'verified_batch.py', 'cohort_sequence_inputs.py')} <= set(found), 'missing probed pipeline sources')
    current = {name: R.S.sha(root / name) for name in found}
    if resume_changed_sources is None:
        for name, digest in found.items():
            R.require(current[name] == digest, 'probed source changed: ' + name)
        return
    declared = set(resume_changed_sources)
    missing = sorted(declared - set(found))
    R.require(not missing, 'declared resume sources absent from probe: ' + ', '.join(missing))
    undeclared = sorted(name for name, digest in found.items()
                        if current[name] != digest and name not in declared)
    R.require(not undeclared, 'undeclared changed probed sources: ' + ', '.join(undeclared))
    return {name: dict(probe_sha256=found[name], current_sha256=current[name])
            for name in sorted(declared)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('preflight', 'full', 'resume'))
    for key in ('out', 'logs', 'jar', 'clean-root', 'clean-audit', 'probe', 'bindings'):
        parser.add_argument('--' + key, type=Path, required=True)
    for key in ('resume-training', 'resume-replay'):
        parser.add_argument('--' + key, type=Path)
    parser.add_argument('--resume-changed-source', action='append', default=[], metavar='RELPATH',
                        help='declare a probed source allowed to differ during resume (repeatable)')
    parser.add_argument('--replay-arm', action='append', default=[], metavar='SOURCE/CONTEXT/KIND',
                        help='audit only this production replay arm (repeatable; resume mode only)')
    parser.add_argument('--restore-wait-seconds', type=float, default=7200,
                        help='maximum wait for another job restoring this results tree')
    args = parser.parse_args()
    R.require(args.mode == 'resume' or not args.resume_changed_source,
              '--resume-changed-source is only valid in resume mode')
    R.require(args.mode == 'resume' or not args.replay_arm,
              '--replay-arm is only valid in resume mode')
    R.require(args.restore_wait_seconds >= 0, '--restore-wait-seconds must be nonnegative')
    parent_job = None
    if args.mode == 'resume':
        R.require(args.resume_training is not None and args.resume_replay is not None,
                  'resume requires sealed training and replay directories')
        R.require(args.resume_training.is_dir() and args.resume_replay.is_dir(),
                  'resume sources must be directories')
        restore_present = (args.out / 'restored.json').exists() or (args.out / 'restore.lock').exists()
        R.require(restore_present or
                  (not (args.out / 'training').exists() and not (args.out / 'replay').exists()),
                  'resume destinations already exist')
        parent_job = resume_parent_job(args.resume_training)
    args.out.mkdir(parents=True, exist_ok=args.mode == 'resume')
    args.logs.mkdir(parents=True, exist_ok=False)
    state = dict(mode='resume-arms' if args.replay_arm else args.mode,
                 job_id=os.environ.get('SLURM_JOB_ID'),
                 hostname=os.uname().nodename, stages=[], status='running')
    if args.mode == 'resume':
        state.update(parent_job=parent_job, resume=dict(training=str(args.resume_training),
                                                          replay=str(args.resume_replay)))
    completion = None
    source_drift = None

    def save():
        target = args.logs / 'workflow.json'
        temporary = target.with_suffix('.partial')
        temporary.write_text(json.dumps(state, indent=2) + '\n')
        temporary.replace(target)

    save()
    try:
        required = [args.jar] + ([args.probe, args.clean_audit] if args.mode in ('full', 'resume') else [])
        prerequisite_inputs = bindings(args.bindings, required)
        training = CACHE / 'observed-gap-models-v1'
        context = CACHE / 'cohort-context-modes-v2-inputs'
        fitting, base, inputs = R.base_inputs(training, context)
        inputs.update(prerequisite_inputs)
        inputs.update(R.C.supervision_inputs(fitting['recipes'],
            json.loads((training / 'registration.json').read_text())['inputs']))
        R.supervision_agreement(fitting['recipes'])
        R.require(len(base['genomes']) == 18, 'complete 18-genome parent required')
        inputs.update(R.F.hashes([args.jar, Path(__file__), *R.C.sources(),
                                 *R.F.V.streaming_parser_sources()]))
        import ijson
        R.require(ijson.__version__ == '3.4.0.post0' and ijson.backend == 'python',
                  'unexpected streaming parser runtime')
        receipt = dict(status='verified', genomes=base['genomes'], inputs=inputs,
                       python=sys.version, numpy=R.F.T.np.__version__,
                       parser=R.F.V.parser_identity(), scope='Parent data and runtime only; no inference or accuracy claim.')
        R.F.save(args.out / 'parent-preflight.json', receipt)
        if args.mode in ('full', 'resume'):
            probe = json.loads(args.probe.read_text())
            R.require(probe.get('schema') == 'gspa.source-context-probe/2' and probe.get('status') == 'verified' and
                      set(probe['summaries']) == {'diamond_homology', 'clean', 'clean_diamond_homology'} and
                      set(probe['historical_model_audit']) == set(R.F.KINDS),
                      'complete three-source production probe required')
            R.require(probe['inputs'].get(str(args.jar)) == R.S.sha(args.jar), 'probe JAR differs')
            if args.mode == 'resume':
                source_drift = verify_probe_sources(probe, resume_changed_sources=args.resume_changed_source)
                state['resume']['source_drift'] = source_drift
                save()
            else:
                verify_probe_sources(probe)
            probe_results = Path(probe['results_root'])
            relocated = {str(args.probe.parent / Path(old).relative_to(probe_results)): digest
                         for old, digest in probe['outputs'].items()}
            R.require(relocated, 'empty probe outputs')
            R.F.check_hashes(relocated)
            inputs.update(relocated)
            R.verify_clean(args.clean_root, args.clean_audit, base['genomes'])
            R.F.save(args.out / 'prerequisites.json', dict(status='verified',
                inputs=R.F.hashes([args.probe, args.clean_audit])))
            learned, replay, analysis = (args.out / x for x in ('training', 'replay', 'analysis'))
            runner = str(ROOT / 'benchmark/v3/run_source_arm_context_comparison.py')
            replay_command = [sys.executable, str(Path(M.__file__)),
                '--training', str(learned), '--out', str(replay)]
            analysis_commands = [
                ('analysis-registration', [sys.executable, str(Path(A.__file__)), 'register',
                 '--training', str(learned), '--replay', str(replay), '--out', str(analysis)]),
                ('analysis', [sys.executable, str(Path(A.__file__)), 'run', '--out', str(analysis)])]
            if args.mode == 'full':
                commands = [('prepare', [sys.executable, runner, 'prepare', '--out', str(learned),
                    '--training', str(training), '--context', str(context), '--jar', str(args.jar),
                    '--clean-root', str(args.clean_root), '--clean-audit', str(args.clean_audit)])]
                commands += [(action, [sys.executable, runner, action, '--out', str(learned)])
                             for action in ('export', 'train', 'evaluate')]
                commands += [('production-replay', replay_command)] + analysis_commands
            else:
                readability = restore_resume_tree(args, learned, replay)
                if readability is not None:
                    state['resume']['readability'] = readability
                    save()
                replay_command.append('--resume')
                if args.replay_arm:
                    for arm in args.replay_arm:
                        replay_command.extend(('--arm', arm))
                    replay_command.extend(('--partial-receipt', str(args.logs / 'replay-partial.json')))
                    commands = [('production-replay', replay_command)]
                else:
                    commands = [('production-replay', replay_command)] + analysis_commands
            state['stages'] = [dict(name=name, command=command, status='queued') for name, command in commands]
            save()
            for stage in state['stages']:
                stage.update(status='running', started_at=datetime.now(timezone.utc).isoformat()); save()
                with (args.logs / (stage['name'] + '.log')).open('x') as log:
                    result = subprocess.run(stage['command'], cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
                stage.update(status='completed' if result.returncode == 0 else 'failed',
                             returncode=result.returncode, finished_at=datetime.now(timezone.utc).isoformat()); save()
                print(json.dumps(stage), flush=True)
                if result.returncode:
                    raise RuntimeError('stage failed: ' + stage['name'])
                # A checkpoint is partial until final verification and publication.
                if not args.replay_arm:
                    subprocess.run(['rsync', '-a', str(args.out) + '/',
                                    str(args.logs.parent / 'checkpoint') + '/'], check=True)
            if not args.replay_arm:
                read = lambda path: json.loads(path.read_text())
                exported, frozen = (read(learned / name) for name in ('export.json', 'training.json'))
                replayed = read(replay / 'verification.json')
                analysis_verified = read(analysis / 'verification.json')
                R.require(analysis_verified.get('status') == 'verified', 'unverified analysis')
                counts = completion_counts(exported, frozen, replayed, read(analysis / 'thresholds.json'),
                                           read(analysis / 'results.json'), base['genomes'])
                verify_frozen_training(learned)
                receipt_paths = [learned / name for name in ('protocol.json', 'registration.json', 'export.json', 'training.json', 'evaluation.json')]
                receipt_paths += [replay / 'verification.json'] + [analysis / name for name in
                    ('protocol.json', 'registration.json', 'thresholds.json', 'results.json', 'verification.json')]
                completion = dict(status='verified', counts=counts, mode=args.mode, parent_job=parent_job,
                                  receipts=R.F.hashes(receipt_paths), model_promotion='none')
                if args.mode == 'resume':
                    completion['source_drift'] = source_drift
        R.F.check_hashes(inputs)
        if completion is not None:
            R.F.save(args.out / 'final-receipt.json', completion)
        state['status'] = 'verified'
    except BaseException as error:
        state.update(status='failed', error=str(error))
        raise
    finally:
        state['finished_at'] = datetime.now(timezone.utc).isoformat()
        save()


if __name__ == '__main__':
    main()
