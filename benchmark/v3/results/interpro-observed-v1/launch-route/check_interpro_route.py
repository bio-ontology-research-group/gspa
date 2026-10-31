"""Run installed InterProScan through annotation and compare direct integration.

Requires an observed fixture captured by generate_interpro_fixture.py and the
same installed tool/database files. Use a compute allocation for real searches.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess

from generate_interpro_fixture import digest


def run(jar, fixture, tool, java, out):
    out.mkdir(parents=True, exist_ok=False)
    provenance = json.loads((fixture / 'provenance.json').read_text())
    for name, expected in provenance['artifacts'].items():
        if digest(fixture / name) != expected:
            raise ValueError('changed observed fixture: ' + name)
    bound = json.loads((fixture / 'inputs.json').read_text())
    for p in (jar, Path(__file__).resolve(), Path(__file__).with_name('generate_interpro_fixture.py'), java):
        bound[str(p)] = digest(p)

    def verify_bound():
        for path, expected in bound.items():
            if digest(path) != expected:
                raise ValueError('changed tool/database/input: ' + path)

    verify_bound()
    config = out / 'config.json'
    config.write_text(json.dumps({'input': {'kingdom': 'bacteria'}, 'predictors': {
        'similarity': {'enabled': False},
        'domains': {'enabled': True, 'tool': 'interproscan', 'disablePrecalc': True,
                    'applications': ['Pfam', 'Hamap', 'ProSiteProfiles', 'ProSitePatterns']},
        'operons': {'enabled': False}, 'pathway': {'enabled': False},
        'localization': {'enabled': False}}, 'quality': {'enabled': False}}, indent=2) + '\n')
    workflow = out / 'workflow.json'
    workflow.write_text(json.dumps({'lite': True, 'enabledPriors': 'none', 'genomeId': 'input'}) + '\n')
    for p in (config, workflow):
        bound[str(p)] = digest(p)
    env = dict(os.environ, PATH=str(tool.parent) + os.pathsep + os.environ['PATH'])
    runs, claims, outcomes = [], {}, {}

    def execute(name, command):
        with (out / (name + '.log')).open('w') as log:
            process = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=3600)
        runs.append({'name': name, 'command': command, 'exit_code': process.returncode})
        (out / 'runs.json').write_text(json.dumps(runs, indent=2) + '\n')
        if process.returncode:
            raise RuntimeError('CLI failed: ' + name)

    for name, integrated in [('annotate', False), ('annotate_integrated', True)]:
        folder = out / name
        command = [str(java), '-XX:ActiveProcessorCount=2', '-jar', str(jar), 'annotate', '--proteins',
            str(fixture / 'input.faa'), '--config', str(config), '--output', str(folder),
            '--claims-out', str(folder / 'raw.jsonl'), '--threads', '1', '--no-metrics']
        if integrated:
            command += ['--integration-config', str(workflow)]
        execute(name, command)
        status = json.loads((folder / 'input_predictor_status.json').read_text())['predictors']
        if len(status) != 1 or status[0]['predictor'] != 'interproscan' or status[0]['outcome'] != 'SUCCESS':
            raise ValueError('InterProScan did not succeed: ' + str(status))
        outcomes[name] = status[0]
        claims[name] = [json.loads(s) for s in (folder / 'raw.jsonl').read_text().splitlines()]
        if not claims[name] or any(c['protein_id'] == 'synthetic_negative' for c in claims[name]):
            raise ValueError('missing matches or unexpected negative-control match')
        kinds = {c['metadata']['statisticType'] for c in claims[name]}
        if not {'evalue', 'profile_score'}.issubset(kinds):
            raise ValueError('fixture did not exercise both statistic kinds')
        for c in claims[name]:
            if not 0 <= c['raw_score'] <= 1 or c['metadata']['scoreCalibrated'] is not False:
                raise ValueError('invalid confidence/provenance')
            if c['function_type'] == 'EC':
                raise ValueError('TSV produced an invented EC annotation')
        log = (out / (name + '.log')).read_text()
        if '-cpu 1 ' not in log or '-dp' not in log:
            raise ValueError('thread or local database setting missing from actual command')
    canonical = lambda rows: sorted(json.dumps(r, sort_keys=True) for r in rows)
    if canonical(claims['annotate']) != canonical(claims['annotate_integrated']):
        raise ValueError('raw claims differ between annotation routes')
    direct = out / 'direct.tsv'
    execute('direct_integrate', [str(java), '-jar', str(jar), 'integrate', '--claims',
        str(out / 'annotate/raw.jsonl'), '--workflow-config', str(workflow), '--out', str(direct)])
    if direct.read_bytes() != (out / 'annotate_integrated/input_integrated.tsv').read_bytes():
        raise ValueError('direct/annotation integration differs')
    verify_bound()
    report = {'schema': 'gspa.interpro-route-replay/1', 'status': 'verified',
              'scope': 'Installed InterProScan on public proteins/control; no accuracy or calibration claim.',
              'predictor_outcomes': outcomes, 'raw_claim_counts': {k: len(v) for k, v in claims.items()},
              'integration_predictions_byte_identical': True, 'inputs': bound,
              'outputs': {str(p.relative_to(out)): digest(p) for p in out.rglob('*') if p.is_file()}}
    (out / 'verification.json').write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ('inputs', 'outputs')}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('jar', 'fixture', 'tool', 'java', 'out'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    run(*(getattr(args, n).resolve() for n in ('jar', 'fixture', 'tool', 'java', 'out')))
