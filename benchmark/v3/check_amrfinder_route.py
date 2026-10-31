"""Verify real AMRFinder execution, typed outcomes, and shared integration parity.

Uses versioned upstream proteins and an independently executed reference table.
This is adapter evidence, not an experimental resistance benchmark.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def run(jar, fixture, tool, database, java, out):
    out.mkdir(parents=True, exist_ok=False)
    provenance = json.loads((fixture / 'provenance.json').read_text())
    for name, digest in provenance['artifacts'].items():
        assert sha(fixture / name) == digest, name
    with (fixture / 'amrfinder_results.tsv').open() as f:
        reference = {r['Protein id']: r for r in csv.DictReader(f, delimiter='\t')}
    files = [jar, java, Path(__file__), *sorted(fixture.iterdir()),
             *sorted(p for p in tool.parent.rglob('*') if p.is_file()),
             *sorted(p for p in database.rglob('*') if p.is_file())]
    files += [Path(shutil.which(exe)).resolve() for exe in ('blastp', 'hmmsearch')]
    bound = {str(p): sha(p) for p in files}
    config = {'input': {'kingdom': 'bacteria'}, 'predictors': {
        'similarity': {'enabled': False}, 'domains': {'enabled': False},
        'pathway': {'enabled': False}, 'operons': {'enabled': False},
        'localization': {'enabled': False}, 'enable': ['amrfinder'],
        'amrfinder': {'threads': 2, 'database': str(database), 'executablePath': str(tool),
                      'organism': 'Escherichia'}}, 'quality': {'enabled': False}}
    workflow = out / 'workflow.json'
    workflow.write_text(json.dumps({'lite': True, 'enabledPriors': 'none', 'genomeId': 'input'}) + '\n')
    bound[str(workflow)] = sha(workflow)
    env = dict(os.environ, TMPDIR=str(out / 'tmp'))
    Path(env['TMPDIR']).mkdir()
    runs, claims, statuses = [], {}, {}

    def execute(name, command):
        with (out / (name + '.log')).open('w') as log:
            result = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=600)
        runs.append({'name': name, 'command': command, 'exit_code': result.returncode})
        (out / 'runs.json').write_text(json.dumps(runs, indent=2) + '\n')
        assert result.returncode == 0, name

    for name in ('annotate', 'annotate_integrated', 'zero', 'invalid_input', 'bad_database', 'unavailable'):
        cfg = json.loads(json.dumps(config))
        fasta = fixture / 'input.faa'
        if name == 'zero':
            fasta = out / 'zero.faa'
            fasta.write_bytes((fixture / 'unreported.faa').read_bytes())
            assert len((fixture / 'unreported.tsv').read_text().splitlines()) == 1
            bound[str(fasta)] = sha(fasta)
        if name == 'invalid_input':
            fasta = out / 'invalid.faa'
            fasta.write_text('>nucleotide_like\nM' + 'A' * 100 + '\n')
            bound[str(fasta)] = sha(fasta)
        if name == 'bad_database':
            cfg['predictors']['amrfinder']['database'] = str(out / 'absent-database')
        if name == 'unavailable':
            cfg['predictors']['amrfinder']['executablePath'] = str(out / 'absent-tool')
        config_file = out / (name + '.json')
        config_file.write_text(json.dumps(cfg, indent=2) + '\n')
        bound[str(config_file)] = sha(config_file)
        folder = out / name
        command = [str(java), '-XX:ActiveProcessorCount=4', '-jar', str(jar), 'annotate',
            '--proteins', str(fasta), '--config', str(config_file), '--output', str(folder),
            '--claims-out', str(folder / 'raw.jsonl'), '--threads', '1', '--no-metrics']
        if name == 'annotate_integrated':
            command += ['--integration-config', str(workflow)]
        execute(name, command)
        status = json.loads((folder / (fasta.stem + '_predictor_status.json')).read_text())['predictors']
        expected = {'zero': 'SUCCESS_ZERO', 'invalid_input': 'FAILED', 'bad_database': 'FAILED', 'unavailable': 'UNAVAILABLE'}.get(name, 'SUCCESS')
        assert len(status) == 1 and status[0]['predictor'] == 'amrfinder' and status[0]['outcome'] == expected, status
        statuses[name] = status[0]
        if name == 'invalid_input':
            assert 'looks like a nucleotide' in status[0]['detail']
        if name == 'bad_database':
            assert 'No valid AMRFinder database' in status[0]['detail']
        rows = [json.loads(line) for line in (folder / 'raw.jsonl').read_text().splitlines()]
        claims[name] = rows
        if expected != 'SUCCESS':
            assert not rows, name
            continue
        assert len(rows) == len(reference) and {r['protein_id'] for r in rows} == set(reference)
        for row in rows:
            r, m = reference[row['protein_id']], row['metadata']
            assert row['function_type'] == 'AMR'
            assert row['function_id'] == r['Type'] + ':' + r['Subtype']
            assert row['raw_score'] == float(r['% Identity to reference']) / 100
            assert m['identity'] == float(r['% Identity to reference'])
            assert m['coverage'] == float(r['% Coverage of reference'])
            assert m['method'] == r['Method'] and m['geneName'] == r['Element symbol']
            assert m['scope'] == r['Scope'] and m['drugClass'] == r['Class']
            assert m['tool_version'] == provenance['tool_version']
            assert m['score_calibrated'] is False and m['score_available'] is True
        commands = [line for line in (out / (name + '.log')).read_text().splitlines() if 'Running amrfinder:' in line]
        assert len(commands) == 1 and '--threads 1 ' in commands[0] and '-d ' + str(database) in commands[0]
    canonical = lambda rows: sorted(json.dumps(r, sort_keys=True) for r in rows)
    assert canonical(claims['annotate']) == canonical(claims['annotate_integrated'])
    direct = out / 'direct.tsv'
    execute('direct_integrate', [str(java), '-jar', str(jar), 'integrate', '--claims',
        str(out / 'annotate/raw.jsonl'), '--workflow-config', str(workflow), '--out', str(direct)])
    assert direct.read_bytes() == (out / 'annotate_integrated/input_integrated.tsv').read_bytes()
    for path, digest in bound.items():
        assert sha(path) == digest, 'input changed: ' + path
    report = {'schema': 'gspa.amrfinder-route-replay/1', 'status': 'verified', 'raw_claims': len(reference),
        'scope': __doc__, 'predictor_outcomes': statuses, 'raw_provenance_identical': True,
        'integration_predictions_byte_identical': True, 'inputs': bound,
        'outputs': {str(p.relative_to(out)): sha(p) for p in out.rglob('*') if p.is_file()}}
    (out / 'verification.json').write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': 'verified', 'raw_claims': len(reference), 'cli_runs': len(runs)}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('jar', 'fixture', 'tool', 'database', 'java', 'out'):
        parser.add_argument('--' + key, type=Path, required=True)
    run(**{k: v.resolve() for k, v in vars(parser.parse_args()).items()})
