#!/usr/bin/env python3
"""Native gene-caller and actual annotation CLI before/after contract checks; not accuracy."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import signal
import sys
import time
import zipfile


def sha(path):
    with Path(path).open('rb') as f:
        h = hashlib.sha256()
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def fasta(path):
    entries = {}
    key = None
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            key = line[1:].split()[0]
            assert key not in entries
            entries[key] = ''
        elif line:
            assert key is not None
            entries[key] += line.strip()
    return entries


def native_reference(path):
    entries = {}
    for block in Path(path).read_text().split('>'):
        if not block:
            continue
        lines = block.splitlines()
        header = lines[0].split(' # ')
        key = (header[0].rsplit('_', 1)[0], int(header[1]), int(header[2]), '+' if header[3] == '1' else '-')
        assert key not in entries
        entries[key] = ''.join(lines[1:])
    gff = Path(path).with_suffix('.gff')
    rows = [line.split('\t') for line in gff.read_text().splitlines() if line and not line.startswith('#')]
    cds = [(r[0], int(r[3]), int(r[4]), r[6]) for r in rows if r[2] == 'CDS']
    assert len(cds) == len(set(cds)) == len(entries) and set(cds) == set(entries)
    return entries


def run(args):
    for name, value in vars(args).items():
        if isinstance(value, Path):
            setattr(args, name, value.resolve())
    out = args.out
    out.mkdir(parents=True, exist_ok=False)
    (out / 'tmp').mkdir()
    libs = out / 'libraries'
    libs.mkdir()
    for wheel in args.wheels.glob('*.whl'):
        with zipfile.ZipFile(wheel) as z:
            z.extractall(libs)
    assert sha(args.before_jar) == '069263f7b876222815928a637d02106b9c3b31a756edd1cf692432acc2e52303'
    assert sha(args.prodigal) == '80a3d7acb08ba9f63dde6b6ef4e95ec356128ac355f36960ea722f381150f246'
    static = [args.jar, args.before_jar, args.prodigal, args.native_json, args.java, Path(__file__), Path("/usr/bin/python3").resolve()]
    static += sorted(args.wheels.glob('*.whl'))
    for directory in (args.observed, args.numeric, args.prodigal_observed):
        static += sorted(p for p in directory.rglob('*') if p.is_file())
    bindings = {str(p): sha(p) for p in static}
    runs = []
    summaries = {}

    def executable(path, body):
        path.write_text(body)
        path.chmod(0o755)
        bindings[str(path)] = sha(path)
        return path

    def environment(name, caller, replay=None):
        directory = out / (name + '-bin')
        directory.mkdir()
        for tool in ('pyrodigal', 'prodigal'):
            path = directory / tool
            if tool != caller:
                executable(path, '#!/bin/bash\nexit 127\n')
            elif replay:
                executable(path, '#!/usr/bin/python3\nimport sys,shutil\na=sys.argv\n'
                    "if '-V' in a or '-v' in a: print('pyrodigal fixture'); sys.exit(0)\n"
                    + 'shutil.copyfile(' + repr(str(replay[0])) + ",a[a.index('-o')+1])\n"
                    + 'shutil.copyfile(' + repr(str(replay[1])) + ",a[a.index('-a')+1])\n")
            elif tool == 'pyrodigal':
                executable(path, '#!/bin/bash\nexport PYTHONPATH=' + shlex.quote(str(libs))
                    + '\nexec /usr/bin/python3 -m pyrodigal "$@"\n')
            else:
                executable(path, '#!/bin/bash\nexec ' + shlex.quote(str(args.prodigal)) + ' "$@"\n')
        return dict(os.environ, PATH=str(directory) + ':' + os.environ['PATH'],
                    TMPDIR=str(out / 'tmp'), JAVA_TOOL_OPTIONS='-Djava.io.tmpdir=' + str(out / 'tmp'))

    def execute(name, command, env):
        started = time.monotonic()
        timed_out = False
        with (out / (name + '.log')).open('w') as log:
            proc = subprocess.Popen(command, env=env, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            try:
                proc.wait(timeout=600)
            except subprocess.TimeoutExpired:
                timed_out = True
                try:
                    os.killpg(proc.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                proc.wait()
        runs.append(dict(name=name, command=command, exit_code=proc.returncode,
                         timed_out=timed_out, elapsed_seconds=time.monotonic() - started))
        (out / 'runs.json').write_text(json.dumps(runs, indent=2) + '\n')
        assert not timed_out and proc.returncode == 0, (name, proc.returncode, timed_out)

    def annotate(name, jar, source, caller, mode, replay=None):
        env = environment(name, caller, replay)
        folder = out / name
        config = out / (name + '.config.json')
        config.write_text(json.dumps(dict(input=dict(kingdom='bacteria', type='metagenome' if mode == 'meta' else 'genome'),
            predictors={key: dict(enabled=False) for key in ('similarity', 'domains', 'pathway', 'operons', 'localization')},
            integration=dict(enabled=False), quality=dict(enabled=False))) + '\n')
        bindings[str(config)] = sha(config)
        execute(name, [str(args.java), '-XX:ActiveProcessorCount=3', '-jar', str(jar), 'annotate',
            '--genome', str(source), '--config', str(config), '--output', str(folder), '--threads', '1', '--no-metrics'], env)
        prefix = source.stem
        status = json.loads((folder / (prefix + '_predictor_status.json')).read_text())['predictors']
        assert len(status) == 1, (name, status)
        with (folder / (prefix + '_proteins.tsv')).open() as stream:
            rows = list(csv.DictReader(stream, delimiter='\t'))
        proteins = fasta(folder / (prefix + '.faa'))
        assert {r['protein_id'] for r in rows} == set(proteins)
        got = {(r['contig'], int(r['start']), int(r['end']), r['strand']): proteins[r['protein_id']] for r in rows}
        assert len(got) == len(rows)
        summaries[name] = dict(proteins=len(rows), status=status[0])
        return folder, prefix, got

    # Whole-assembly parity uses each native tool's OWN complete output as reference.
    for caller, observed in [('pyrodigal', args.observed), ('prodigal', args.prodigal_observed)]:
        for mode in ('single', 'meta'):
            source = observed / 'inputs/genome.fna'
            reference = native_reference(observed / ('genome-' + mode + '.faa'))
            assert reference
            old = annotate('old-' + caller + '-' + mode, args.before_jar, source, caller, mode)
            new = annotate('new-' + caller + '-' + mode, args.jar, source, caller, mode)
            assert old[2] == new[2] == reference
            for suffix in ('.faa', '_proteins.tsv', '.gff3', '_annotations.tsv'):
                assert (old[0] / (old[1] + suffix)).read_bytes() == (new[0] / (new[1] + suffix)).read_bytes(), (caller, mode, suffix)
            assert summaries['new-' + caller + '-' + mode]['status']['outcome'] == 'SUCCESS'
            assert summaries['new-' + caller + '-' + mode]['status']['exit_code'] == 0
        for mode in ('single', 'meta'):
            for version, jar in [('old', args.before_jar), ('new', args.jar)]:
                name = version + '-' + caller + '-short-' + mode
                _, _, got = annotate(name, jar, observed / 'inputs/short.fna', caller, mode)
                assert not got
                expected = 'SUCCESS_ZERO' if version == 'new' and mode == 'meta' else 'FAILED'
                assert summaries[name]['status']['outcome'] == expected, summaries[name]
                if version == 'new':
                    native = json.loads((observed / 'observations.json').read_text())['cases']
                    native_exit = next(r['exit_code'] for r in native if r['name'] == 'short-' + mode)
                    assert summaries[name]['status']['exit_code'] == native_exit

    for mode in ('single', 'meta'):
        name = 'numeric-' + mode
        _, _, got = annotate(name, args.jar, args.numeric / 'inputs/numeric.fna', 'pyrodigal', mode)
        assert got == native_reference(args.numeric / ('numeric-' + mode + '.faa'))
        assert summaries[name]['status']['outcome'] == 'SUCCESS'

    original = args.observed / 'genome-single.faa'
    records = ['>' + r for r in original.read_text().split('>') if r]
    for name, selected, outcome in [('reordered', list(reversed(records)), 'SUCCESS'), ('missing-protein', records[1:], 'FAILED')]:
        changed = out / (name + '.faa')
        changed.write_text(''.join(selected))
        bindings[str(changed)] = sha(changed)
        folder, prefix, got = annotate(name, args.jar, args.observed / 'inputs/genome.fna', 'pyrodigal', 'single',
            (args.observed / 'genome-single.gff', changed))
        assert summaries[name]['status']['outcome'] == outcome
        if outcome == 'SUCCESS':
            assert got == native_reference(original)
        else:
            assert not got
            gff = (folder / (prefix + '.gff3')).read_text().split('##FASTA')[0]
            assert not any('\tCDS\t' in line for line in gff.splitlines())

    _, _, got = annotate('unavailable', args.jar, args.observed / 'inputs/short.fna', 'none', 'meta')
    assert not got and summaries['unavailable']['status']['outcome'] == 'UNAVAILABLE'

    # The extracted coordinate reader still yields byte-identical imported native BGC claims.
    for version, jar in [('old', args.before_jar), ('new', args.jar)]:
        env = dict(os.environ, TMPDIR=str(out / 'tmp'))
        execute('import-' + version, [str(args.java), '-jar', str(jar), 'antismash-claims',
            '--results', str(args.native_json), '--genome', str(args.observed / 'inputs/genome.fna'),
            '--gff3', str(args.observed / 'genome-single.gff'), '--proteins', str(original),
            '--output', str(out / ('import-' + version + '.jsonl'))], env)
    assert (out / 'import-old.jsonl').read_bytes() == (out / 'import-new.jsonl').read_bytes()
    claims = [json.loads(line) for line in (out / 'import-new.jsonl').read_text().splitlines()]
    assert claims and all(r['metadata']['score_available'] is False for r in claims)
    for path, digest in bindings.items():
        assert sha(path) == digest, path
    report = dict(status='verified', scope=__doc__, runs=len(runs), annotations=summaries,
                  runtime=dict(checker_python=sys.version, native_python=subprocess.check_output(['/usr/bin/python3', '--version'], text=True).strip()),
                  imported_unscored_memberships=len(claims), bindings=bindings,
                  output_hashes={str(p.relative_to(out)): sha(p) for p in out.rglob('*') if p.is_file() and 'libraries' not in p.parts})
    (out / 'verification.json').write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps(dict(status=report['status'], runs=len(runs), imported_unscored_memberships=len(claims))))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('jar', 'before-jar', 'java', 'observed', 'numeric', 'prodigal-observed', 'prodigal', 'wheels', 'native-json', 'out'):
        parser.add_argument('--' + name, type=Path, required=True)
    run(parser.parse_args())
