"""Replay AntiFam annotation availability using real hmmsearch and a synthetic HMM.

This does not establish biological spurious-ORF detection or functional rescoring.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess

from check_hmmer_route import digest


def run(before, after, fixture, out, hmmsearch):
    out.mkdir(parents=True, exist_ok=False)
    provenance = json.loads((fixture / 'provenance.json').read_text())
    for name, expected in provenance['artifacts'].items():
        if digest(fixture / name) != expected:
            raise ValueError('fixture changed: ' + name)
    if digest(hmmsearch) != provenance['tools']['hmmsearch']['sha256']:
        raise ValueError('tool differs from observed fixture')
    config = out / 'config.json'
    config.write_text(json.dumps({'input': {'kingdom': 'bacteria'}, 'quality': {'enabled': False},
        'predictors': {**{name: {'enabled': False} for name in
                         ['similarity', 'domains', 'operons', 'pathway', 'localization']},
                       'antifam': {'enabled': True, 'database': str(fixture / 'profile.hmm'),
                                   'useCutGa': False}}}, indent=2) + '\n')
    java = Path(shutil.which('java')).resolve()
    paths = [before, after, fixture / 'profile.hmm', fixture / 'input.faa', fixture / 'provenance.json',
             config, hmmsearch, java, Path(__file__).resolve(), Path(__file__).with_name('check_hmmer_route.py').resolve()]
    inputs = {str(p): digest(p) for p in paths}
    runs = []
    env = dict(os.environ, PATH=str(hmmsearch.parent) + os.pathsep + os.environ['PATH'])
    for arm, jar in [('before', before), ('after', after)]:
        directory = out / arm
        command = [str(java), '-jar', str(jar), 'annotate', '--proteins', str(fixture / 'input.faa'),
                   '--config', str(config), '--output', str(directory), '--claims-out',
                   str(directory / 'raw.jsonl'), '--threads', '1', '--no-metrics']
        with (out / (arm + '.log')).open('w') as log:
            process = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=120)
        if process.returncode:
            raise RuntimeError('CLI failed: ' + arm)
        statuses = json.loads((directory / 'input_predictor_status.json').read_text())['predictors']
        if len(statuses) != 1 or statuses[0]['predictor'] != 'antifam':
            raise ValueError('wrong predictor set')
        with (directory / 'input_annotations.tsv').open() as stream:
            annotations = list(csv.DictReader(stream, delimiter='\t'))
        if arm == 'before':
            if statuses[0]['outcome'] != 'UNAVAILABLE' or annotations:
                raise ValueError('old build did not reproduce unavailable AntiFam')
        else:
            if statuses[0]['outcome'] != 'SUCCESS' or len(annotations) != 1:
                raise ValueError('wrong repaired AntiFam outcome')
            hit = annotations[0]
            if (hit['protein_id'], hit['type'], hit['value'], hit['source']) != (
                    'target', 'PSEUDOGENE', 'GSPA_TEST_001', 'antifam'):
                raise ValueError('wrong annotation or negative-control hit')
            if '--cpu 1 ' not in (out / (arm + '.log')).read_text():
                raise ValueError('thread budget not honored')
        # Observe, but do not certify functional integration of this quality flag.
        raw_count = len((directory / 'raw.jsonl').read_text().splitlines())
        runs.append({'arm': arm, 'command': command, 'exit_code': process.returncode,
                     'predictor_outcome': statuses[0]['outcome'],
                     'annotation_count': len(annotations), 'raw_claim_count': raw_count})
    if inputs != {str(p): digest(p) for p in paths}:
        raise ValueError('inputs changed during replay')
    report = {'schema': 'gspa.antifam-route-replay/1', 'status': 'verified', 'runs': runs,
              'scope': 'Synthetic-model execution only; no real AntiFam detection or functional integration claim.',
              'inputs': inputs, 'outputs': {str(p.relative_to(out)): digest(p)
                  for p in sorted(out.rglob('*')) if p.is_file()}}
    (out / 'verification.json').write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps([{k: v for k, v in r.items() if k != 'command'} for r in runs]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('before', 'after', 'fixture', 'out', 'hmmsearch'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    run(*(getattr(args, name).resolve() for name in ('before', 'after', 'fixture', 'out', 'hmmsearch')))
