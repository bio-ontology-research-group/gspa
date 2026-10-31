"""Reproduce a pipe-containing protein identifier crash through production CLI builds.

Synthetic essential GO claim; no biological annotation-accuracy assertion.
"""
import argparse
import csv
import json
from pathlib import Path
import shutil
import subprocess

from check_hmmer_route import digest


def run(before, after, out):
    out.mkdir(parents=True, exist_ok=False)
    claims = {}
    for label, protein in [('plain', 'protein'), ('pipe', 'sp|P00956|SYI_ECOLI')]:
        path = out / (label + '.jsonl')
        path.write_text(json.dumps({'protein_id': protein, 'function_type': 'GO',
            'function_id': 'GO:0006412', 'source': 'diamond', 'evidence_type': 'SEQUENCE_SIMILARITY',
            'raw_score': .4, 'calibrated_prob': .4}) + '\n')
        claims[label] = path
    java = Path(shutil.which('java')).resolve()
    paths = [before, after, java, *claims.values(), Path(__file__).resolve(),
             Path(__file__).with_name('check_hmmer_route.py').resolve()]
    inputs = {str(p): digest(p) for p in paths}
    runs = []
    predictions = {}
    for arm, jar in [('before', before), ('after', after)]:
        for label, claim in claims.items():
            name = arm + '-' + label
            output = out / (name + '.tsv')
            command = [str(java), '-jar', str(jar), 'integrate', '--claims', str(claim),
                '--out', str(output), '--lite', '--essential-profile', 'bacteria',
                '--enable-priors', 'essentiality']
            with (out / (name + '.log')).open('w') as log:
                process = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=60)
            runs.append({'name': name, 'command': command, 'exit_code': process.returncode})
            if name == 'before-pipe':
                if process.returncode == 0 or 'NullPointerException' not in (out / (name + '.log')).read_text():
                    raise ValueError('old build did not reproduce the targeted-rescue crash')
                continue
            if process.returncode:
                raise RuntimeError('CLI failed: ' + name)
            with output.open() as stream:
                rows = list(csv.DictReader(stream, delimiter='\t'))
            if len(rows) != 1:
                raise ValueError('changed candidate count')
            if rows[0]['protein_id'] != json.loads(claim.read_text())['protein_id']:
                raise ValueError('protein identity changed')
            predictions[name] = {k: v for k, v in rows[0].items() if k != 'protein_id'}
    if predictions['before-plain'] != predictions['after-plain']:
        raise ValueError('ordinary identifier behavior changed')
    if predictions['after-plain'] != predictions['after-pipe']:
        raise ValueError('protein spelling changes the repaired prediction')
    if inputs != {str(p): digest(p) for p in paths}:
        raise ValueError('input changed during replay')
    report = {'schema': 'gspa.claim-key-replay/1', 'status': 'verified',
        'scope': 'Synthetic identity/crash regression, no biological accuracy claim.',
        'runs': runs, 'predictions_except_protein_id': predictions, 'inputs': inputs,
        'outputs': {str(p.relative_to(out)): digest(p) for p in out.rglob('*') if p.is_file()}}
    (out / 'verification.json').write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': 'verified', 'exit_codes': {r['name']: r['exit_code'] for r in runs},
                      'ordinary_id_unchanged': True, 'repaired_prediction_id_invariant': True}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('before', 'after', 'out'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    run(args.before.resolve(), args.after.resolve(), args.out.resolve())
