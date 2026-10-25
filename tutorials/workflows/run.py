#!/usr/bin/env python3
"""A small, standard-library launcher; long container commands stay out of lessons."""
import argparse
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import webbrowser

ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / '.tutorial-work'
HELPERS = ROOT / 'tutorials/workflows'
CASES = ('bacteria', 'archaea', 'eukaryote', 'mag', 'metagenome', 'virus', 'prokaryote')


def execute(command, log=None):
    if log is None:
        subprocess.run(command, cwd=ROOT, check=True)
        return
    log.parent.mkdir(parents=True, exist_ok=True)
    print(f'Working… log: {log.relative_to(ROOT)}', flush=True)
    with log.open('w') as handle:
        subprocess.run(command, cwd=ROOT, stdout=handle, stderr=subprocess.STDOUT, check=True)


def docker_ready():
    if not shutil.which('docker'):
        raise SystemExit('Docker is required. See docs/tutorials/install.md, then run ./tutorial setup.')
    result = subprocess.run(['docker', 'info'], capture_output=True, text=True)
    if result.returncode:
        raise SystemExit('Docker is not available. Start Docker and check that `docker info` works.\n' + result.stderr.strip())


def container(command, image='cpu'):
    args = ['docker', 'run', '--rm', '--user', f'{os.getuid()}:{os.getgid()}',
            '--cpus', '2', '--memory', '6g']
    if image == 'neural-cuda':
        args += ['--gpus', 'all']
    return args + ['-v', f'{WORK}:/work', '-v', f'{HELPERS}:/tutorial:ro',
                   f'gspa-tutorial:{image}'] + command


def require_setup():
    docker_ready()
    image = subprocess.run(['docker', 'image', 'inspect', 'gspa-tutorial:cpu'], capture_output=True)
    needed = ('reference.dmnd', 'cpu.json', 'virus.json', 'downloads.json')
    if image.returncode or not all((WORK / 'data' / name).is_file() for name in needed):
        raise SystemExit('Prepare the teaching environment first: ./tutorial setup')


def setup():
    docker_ready()
    WORK.mkdir(exist_ok=True)
    print('Building the CPU environment. First setup needs internet and may take several minutes.', flush=True)
    execute([sys.executable, str(HELPERS / 'prepare_context.py'), str(WORK)])
    execute(['docker', 'build', '-t', 'gspa-tutorial:cpu', '-f',
             str(WORK / 'source/tutorials/workflows/Dockerfile'), str(WORK / 'source')], WORK / 'build-cpu.log')
    execute(container(['gspa', '--version']))
    print('Preparing and checking the teaching sequences…', flush=True)
    execute(container(['python3', '/tutorial/data.py']))
    execute(container(['diamond', 'makedb', '--in', '/work/data/reference.faa',
                       '--db', '/work/data/reference', '--threads', '2']), WORK / 'makedb.log')
    if not (WORK / 'data/reference.dmnd').stat().st_size:
        raise SystemExit('The teaching database is empty; inspect .tutorial-work/makedb.log.')
    print('Ready. Start the first lesson: ./tutorial run bacteria')


def run_neural(device):
    require_setup()
    if device == 'cuda':
        required = ('centroids.npz', 'config.json', 'query.faa')
        if (not all((WORK / 'neural' / name).is_file() for name in required)
                or not list((WORK / 'neural/cpu').glob('*_annotations.tsv'))):
            raise SystemExit('Run the CPU comparison first: ./tutorial neural cpu')
    (WORK / 'neural/comparison.json').unlink(missing_ok=True)
    tag = 'neural-' + device
    command = ['docker', 'build', '-t', 'gspa-tutorial:' + tag,
               '-f', str(HELPERS / 'Dockerfile.neural')]
    if device == 'cuda':
        command += ['--build-arg', 'TORCH_INDEX=https://download.pytorch.org/whl/cu128']
    execute(command + [str(HELPERS)], WORK / f'build-{tag}.log')
    if device == 'cpu':
        execute(container(['gspa', 'download', 'esm2-35m', '--dir', '/work/resources']), WORK / 'download-esm2.log')
        # The symlink target must be the container path, not the host path.
        cache = WORK / 'torch/hub/checkpoints'
        cache.mkdir(parents=True, exist_ok=True)
        for path in (WORK / 'resources/esm2-35m').glob('*.pt'):
            link = cache / path.name
            link.unlink(missing_ok=True)
            link.symlink_to('/work/resources/esm2-35m/' + path.name)
    execute(container(['python3', '/tutorial/workflow.py', tag], image=tag))
    print('Results and logs: .tutorial-work/neural/')


def main():
    parser = argparse.ArgumentParser(description='Explore protein function with the GSPA teaching datasets.')
    sub = parser.add_subparsers(dest='action', required=True)
    sub.add_parser('setup', help='build the CPU environment and prepare the small teaching reference')
    for name, help_text in [('run', 'run a biological lesson and generate its browser'),
                            ('inspect', 'print the biological observations from existing results'),
                            ('browser', 'open the results in your web browser'),
                            ('commands', 'show the GSPA commands behind a lesson')]:
        child = sub.add_parser(name, help=help_text)
        child.add_argument('case', choices=CASES)
    child = sub.add_parser('neural', help='run the optional ESM-2 lesson (CPU first, then optional CUDA)')
    child.add_argument('device', choices=('cpu', 'cuda'))
    args = parser.parse_args()
    os.chdir(ROOT)
    try:
        if args.action == 'setup':
            setup()
        elif args.action == 'run':
            require_setup()
            execute(container(['python3', '/tutorial/workflow.py', args.case]))
            print(f'Results: .tutorial-work/{args.case}/\nOpen them: ./tutorial browser {args.case}')
        elif args.action == 'neural':
            run_neural(args.device)
        elif args.action == 'commands':
            from workflow import commands
            for stage, command in commands(args.case):
                print(f'# {stage}\n{shlex.join(command)}\n')
            print('These paths are inside the tutorial container; /work maps to .tutorial-work/.')
            print('The workflow also converts claims, summarizes results and prepares browser metadata.')
        elif args.action == 'inspect':
            if not (WORK / args.case / 'summary.json').exists():
                raise SystemExit(f'Run the lesson first: ./tutorial run {args.case}')
            import results
            results.inspect(args.case)
        elif args.action == 'browser':
            path = WORK / args.case / 'gspa_out' / f'{args.case}_browser.html'
            if not path.is_file() or path.stat().st_size == 0:
                raise SystemExit(f'Run the lesson first: ./tutorial run {args.case}')
            print(f'Results browser: {path.as_uri()}', flush=True)
            if not webbrowser.open(path.as_uri()):
                print('No desktop browser detected. Open the HTML file on your local computer.')
    except subprocess.CalledProcessError as error:
        print(f'Command failed (exit {error.returncode}): {shlex.join(error.cmd)}', file=sys.stderr)
        print('See the log path printed above for details. Fix the problem, then rerun the same command.', file=sys.stderr)
        raise SystemExit(error.returncode)
    except OSError as error:
        raise SystemExit(str(error))


if __name__ == '__main__':
    main()
