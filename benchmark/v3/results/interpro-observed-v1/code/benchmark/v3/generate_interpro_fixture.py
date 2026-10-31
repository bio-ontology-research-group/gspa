"""Capture genuine installed InterProScan 5 output and input/tool identities.

Run on a compute allocation. Inputs are public proteins plus a declared synthetic
control; output is parser/dispatch evidence, not biological accuracy truth.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def capture(tool, fasta, out):
    out.mkdir(parents=True, exist_ok=False)
    install = tool.resolve().parent
    shutil.copy2(fasta, out / 'input.faa')
    files = {tool.resolve(), install / 'interproscan-5.jar', install / 'interproscan.properties',
             Path(__file__).resolve(), fasta.resolve(), Path(shutil.which('java')).resolve()}
    for directory in ('data/pfam', 'data/hamap', 'data/prosite', 'bin/prosite', 'bin/hmmer', 'lib'):
        files.update(p for p in (install / directory).rglob('*') if p.is_file())
    identity = {str(p): digest(p) for p in sorted(files)}
    (out / 'inputs.json').write_text(json.dumps(identity, indent=2, sort_keys=True) + '\n')
    records = []
    for name, args in [('double-dash-version', ['--version']), ('version', ['-version']),
                       ('help', ['--help']), ('matches', ['-i', str(out / 'input.faa'), '-o',
                        str(out / 'interproscan_results.tsv'), '-f', 'TSV', '-appl',
                        'Pfam,Hamap,ProSiteProfiles,ProSitePatterns', '-goterms', '-pa', '-dp',
                        '-cpu', '2', '-T', str(Path(os.environ['TMPDIR']) / 'interpro-work')])]:
        command = [str(tool), *args]
        start = time.monotonic()
        with (out / (name + '.stdout.txt')).open('w') as stdout, (out / (name + '.stderr.txt')).open('w') as stderr:
            proc = subprocess.run(command, cwd=out, stdout=stdout, stderr=stderr, timeout=3600)
        records.append({'name': name, 'argv': command, 'cwd': str(out), 'exit_code': proc.returncode,
                        'elapsed_seconds': time.monotonic() - start})
        (out / 'commands.json').write_text(json.dumps(records, indent=2) + '\n')
        if name == 'matches' and proc.returncode:
            raise RuntimeError('actual InterProScan execution failed; see retained logs')
    for path, expected in identity.items():
        if digest(path) != expected:
            raise RuntimeError('input/tool/database changed during run: ' + path)
    report = {'schema': 'gspa.observed-tool-fixture/1', 'tool': 'InterProScan 5',
              'scope': 'Real execution on public proteins and synthetic control; not accuracy truth.',
              'identity_scope': 'launcher, main jar, libraries, properties, selected database directories and binaries, Java, generator and input',
              'artifacts': {p.name: digest(p) for p in sorted(out.iterdir()) if p.is_file()}}
    (out / 'provenance.json').write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': 'captured', 'out': str(out), 'bound_files': len(identity)}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tool', type=Path, required=True)
    parser.add_argument('--fasta', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    capture(args.tool.resolve(), args.fasta.resolve(), args.out.resolve())
