"""Replay genuine HMMER through old/new annotation CLI and shared integration.

The supplied observed fixture uses synthetic sequences and a synthetic model.
This is an execution/identity check, not a biological accuracy measurement.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def run(before_jar, after_jar, fixture, out, hmmsearch):
    out.mkdir(parents=True, exist_ok=False)
    provenance = json.loads((fixture / 'provenance.json').read_text())
    for name, expected in provenance['artifacts'].items():
        if digest(fixture / name) != expected:
            raise ValueError('changed observed fixture: ' + name)
    if digest(hmmsearch) != provenance['tools']['hmmsearch']['sha256']:
        raise ValueError('tool differs from the observed fixture; regenerate it first')
    config = out / 'config.json'
    config.write_text(json.dumps({
        'input': {'kingdom': 'bacteria'},
        'predictors': {
            'similarity': {'enabled': False},
            'domains': {'enabled': True, 'tool': 'hmmer', 'database': str(fixture / 'profile.hmm')},
            'operons': {'enabled': False}, 'pathway': {'enabled': False},
            'localization': {'enabled': False}},
        'quality': {'enabled': False}}, indent=2) + '\n')
    workflow = out / 'workflow.json'
    workflow.write_text(json.dumps({'lite': True, 'enabledPriors': 'none', 'genomeId': 'input'}) + '\n')
    java = Path(shutil.which('java')).resolve()
    inputs = [before_jar, after_jar, fixture / 'profile.hmm', fixture / 'input.faa',
              fixture / 'provenance.json', config, workflow, hmmsearch, java, Path(__file__).resolve()]
    initial = {str(p): digest(p) for p in inputs}
    env = dict(os.environ, PATH=str(hmmsearch.parent) + os.pathsep + os.environ['PATH'])
    runs = []

    def execute(name, cmd):
        with (out / (name + '.log')).open('w') as log:
            process = subprocess.run(cmd, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=120)
        runs.append({'name': name, 'command': cmd, 'exit_code': process.returncode})
        (out / 'runs.json').write_text(json.dumps(runs, indent=2) + '\n')
        if process.returncode:
            raise RuntimeError('CLI failed: ' + name)

    claims = {}
    statuses = {}
    for name, jar, integrated in [('before', before_jar, False), ('after', after_jar, False),
                                   ('annotate_integrated', after_jar, True)]:
        directory = out / name
        command = [str(java), '-jar', str(jar), 'annotate', '--proteins', str(fixture / 'input.faa'),
                   '--config', str(config), '--output', str(directory), '--claims-out',
                   str(directory / 'raw.jsonl'), '--threads', '1', '--no-metrics']
        if integrated:
            command += ['--integration-config', str(workflow)]
        execute(name, command)
        status = json.loads((directory / 'input_predictor_status.json').read_text())['predictors']
        if len(status) != 1 or status[0]['predictor'] != 'hmmer':
            raise ValueError('unexpected predictor set')
        statuses[name] = status[0]['outcome']
        claims[name] = [json.loads(line) for line in (directory / 'raw.jsonl').read_text().splitlines()]
        if name == 'before':
            if statuses[name] != 'UNAVAILABLE' or claims[name]:
                raise ValueError('old build did not reproduce the availability defect')
        else:
            if statuses[name] != 'SUCCESS' or len(claims[name]) != 1:
                raise ValueError('HMMER did not produce exactly one successful hit')
            hit = claims[name][0]
            if (hit['protein_id'], hit['function_type'], hit['function_id'], hit['source']) != (
                    'target', 'PFAM', 'GSPA_TEST_001', 'hmmer'):
                raise ValueError('wrong claim identity or negative-control hit')
            if hit['metadata']['domainEvalue'] != 2.7e-65 or hit['metadata']['domainScore'] != 204.0:
                raise ValueError('wrong domain statistics')
            if '--cpu 1 ' not in (out / (name + '.log')).read_text():
                raise ValueError('HMMER ignored the CLI thread budget')
    if claims['after'] != claims['annotate_integrated']:
        raise ValueError('integration altered the captured raw claims')
    direct = out / 'direct.tsv'
    execute('direct_integrate', [str(java), '-jar', str(after_jar), 'integrate',
            '--claims', str(out / 'after/raw.jsonl'), '--workflow-config', str(workflow), '--out', str(direct)])
    annotated = out / 'annotate_integrated/input_integrated.tsv'
    if direct.read_bytes() != annotated.read_bytes():
        raise ValueError('annotation and direct integration predictions differ')
    with direct.open() as stream:
        if len(list(csv.DictReader(stream, delimiter='\t'))) != 1:
            raise ValueError('integration lost or added candidates')
    if {str(p): digest(p) for p in inputs} != initial:
        raise ValueError('input changed during replay')
    report = {'schema': 'gspa.hmmer-route-replay/1', 'status': 'verified',
              'scope': 'Real HMMER on synthetic inputs; no biological accuracy claim.',
              'predictor_outcomes': statuses, 'raw_claim_counts': {k: len(v) for k, v in claims.items()},
              'integration_predictions_byte_identical': True, 'inputs': initial,
              'outputs': {str(p.relative_to(out)): digest(p) for p in out.rglob('*') if p.is_file()}}
    (out / 'verification.json').write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ('inputs', 'outputs')}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('before-jar', 'after-jar', 'fixture', 'out', 'hmmsearch'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    run(*(getattr(args, name).resolve() for name in ('before_jar', 'after_jar', 'fixture', 'out', 'hmmsearch')))
