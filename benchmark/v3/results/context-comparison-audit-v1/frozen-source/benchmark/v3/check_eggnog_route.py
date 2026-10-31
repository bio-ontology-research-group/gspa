"""Replay installed eggNOG 2/3 through annotation and direct integration.

Uses actual upstream miniature databases and independently captured tool tables.
Checks execution, provenance and route parity; does not establish biological accuracy.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def expected_claims(table):
    """Named columns from the separately executed reference tool, before GSPA."""
    expected, header = [], None
    for line in table.read_text().splitlines():
        if line.startswith('#query'):
            header = line[1:].split('\t')
        elif line and not line.startswith('#'):
            if header is None:
                raise ValueError('reference header missing')
            row = dict(zip(header, line.split('\t')))
            e = float(row['evalue'])
            score = 1. if e == 0 else max(0., min(1., -math.log10(e) / 30.))
            for field, kind, prefix in [('GOs', 'GO', ''), ('EC', 'EC', 'EC:'), ('KEGG_ko', 'KEGG', '')]:
                if row[field] != '-':
                    for value in row[field].split(','):
                        expected.append((row['query'], kind, prefix + value.strip(), score))
            if row['COG_category'] not in ('-', 'S'):
                expected.append((row['query'], 'COG', row['COG_category'], score))
    return sorted(expected)


def project(claims):
    return sorted((c['protein_id'], c['function_type'], c['function_id'], c['raw_score']) for c in claims)


def run(jar, fixture, tool, database, java, out, old_claims=None):
    out.mkdir(parents=True, exist_ok=False)
    provenance = json.loads((fixture / 'provenance.json').read_text())
    for name, expected in provenance['artifacts'].items():
        if sha(fixture / name) != expected:
            raise ValueError('changed reference fixture: ' + name)
    reference = expected_claims(fixture / 'eggnog_results.emapper.annotations')
    bound = {str(p): sha(p) for p in [jar, java, tool, Path(__file__),
             *sorted(fixture.iterdir()), *sorted(p for p in database.rglob('*') if p.is_file())]}
    # Bind the installed Python implementation as well as its launcher.
    for folder in tool.parent.parent.glob('lib/python*/site-packages/eggnogmapper'):
        for p in sorted(folder.rglob('*')):
            if p.is_file() and p.suffix in ('.py', '.so'):
                bound[str(p)] = sha(p)
    diamond = Path(shutil.which('diamond'))
    bound[str(diamond)] = sha(diamond)
    config = out / 'config.json'
    config.write_text(json.dumps({'input': {'kingdom': 'bacteria'}, 'predictors': {
        'similarity': {'enabled': False}, 'domains': {'enabled': False},
        'pathway': {'enabled': False}, 'operons': {'enabled': False},
        'localization': {'enabled': False}, 'enable': ['eggnog-mapper'],
        'eggnog': {'threads': 2, 'dataDir': str(database), 'executablePath': str(tool)}},
        'quality': {'enabled': False}}, indent=2) + '\n')
    workflow = out / 'workflow.json'
    workflow.write_text(json.dumps({'lite': True, 'enabledPriors': 'none', 'genomeId': 'input'}) + '\n')
    for p in (config, workflow):
        bound[str(p)] = sha(p)
    env = dict(os.environ)
    env.pop('EGGNOG_DATA_DIR', None)  # Explicit configuration must suffice.
    runs, claims, statuses = [], {}, {}

    def execute(name, command):
        with (out / (name + '.log')).open('w') as log:
            result = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=600)
        runs.append({'name': name, 'command': command, 'exit_code': result.returncode})
        (out / 'runs.json').write_text(json.dumps(runs, indent=2) + '\n')
        if result.returncode:
            raise RuntimeError('CLI failed: ' + name)

    for name, integrated in [('annotate', False), ('annotate_integrated', True)]:
        folder = out / name
        command = [str(java), '-XX:ActiveProcessorCount=4', '-jar', str(jar), 'annotate',
            '--proteins', str(fixture / 'input.faa'), '--config', str(config), '--output', str(folder),
            '--claims-out', str(folder / 'raw.jsonl'), '--threads', '1', '--no-metrics']
        if integrated:
            command += ['--integration-config', str(workflow)]
        execute(name, command)
        status = json.loads((folder / 'input_predictor_status.json').read_text())['predictors']
        if len(status) != 1 or status[0]['predictor'] != 'eggnog-mapper' or status[0]['outcome'] != 'SUCCESS':
            raise ValueError('eggNOG failed: ' + str(status))
        statuses[name] = status[0]
        rows = [json.loads(line) for line in (folder / 'raw.jsonl').read_text().splitlines()]
        claims[name] = rows
        if project(rows) != reference:
            raise ValueError('claims disagree with independent reference table')
        for row in rows:
            m = row['metadata']
            if m['scoreCalibrated'] is not False or m['statisticType'] != 'evalue' or not m['seedOrtholog']:
                raise ValueError('score/seed provenance lost')
            if not m['toolVersion'].startswith('emapper-' + provenance['version']):
                raise ValueError('wrong reported tool version')
        commands = [line for line in (out / (name + '.log')).read_text().splitlines()
                    if 'Running eggnog-mapper:' in line]
        if len(commands) != 1 or '--cpu 1 ' not in commands[0] or '--data_dir ' + str(database) not in commands[0]:
            raise ValueError('thread/database configuration did not reach tool')
        if ('--go_evidence non-electronic' in commands[0]) != provenance['version'].startswith('2.'):
            raise ValueError('wrong version-specific donor filter')
    canonical = lambda rows: sorted(json.dumps(r, sort_keys=True) for r in rows)
    if canonical(claims['annotate']) != canonical(claims['annotate_integrated']):
        raise ValueError('raw provenance differs across annotation routes')
    direct = out / 'direct.tsv'
    execute('direct_integrate', [str(java), '-jar', str(jar), 'integrate', '--claims',
        str(out / 'annotate/raw.jsonl'), '--workflow-config', str(workflow), '--out', str(direct)])
    if direct.read_bytes() != (out / 'annotate_integrated/input_integrated.tsv').read_bytes():
        raise ValueError('direct/annotation integration differs')
    if old_claims:
        bound[str(old_claims)] = sha(old_claims)
        old = [json.loads(line) for line in old_claims.read_text().splitlines()]
        if project(old) != reference:
            raise ValueError('version-2 supported claims/scores changed')
    for path, expected in bound.items():
        if sha(path) != expected:
            raise ValueError('input changed during replay: ' + path)
    report = {'schema': 'gspa.eggnog-route-replay/1', 'status': 'verified',
        'version': provenance['version'], 'scope': 'Miniature upstream database execution and route parity; no biological accuracy claim.',
        'raw_claims': len(reference), 'predictor_outcomes': statuses,
        'raw_provenance_identical': True, 'integration_predictions_byte_identical': True,
        'old_supported_claims_preserved': True if old_claims else None, 'inputs': bound,
        'outputs': {str(p.relative_to(out)): sha(p) for p in out.rglob('*') if p.is_file()}}
    (out / 'verification.json').write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ('inputs', 'outputs', 'predictor_outcomes')}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('jar', 'fixture', 'tool', 'database', 'java', 'out'):
        parser.add_argument('--' + key, type=Path, required=True)
    parser.add_argument('--old-claims', type=Path)
    args = parser.parse_args()
    run(**{k: v.resolve() if v else None for k, v in vars(args).items()})
