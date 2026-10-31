#!/usr/bin/env python3
"""Actual production CLI and independent LP replay on explicit synthetic mechanisms."""
import argparse
import copy
import json
from pathlib import Path
import subprocess
import sys

from check_absence_cobra import compare, sha

ROOT = Path(__file__).resolve().parents[2]


def run(jar, out):
    out.mkdir(parents=True, exist_ok=False)
    helper = ROOT / 'gspa-core/src/test/resources/experiments/make_fixtures.py'
    subprocess.run([sys.executable, '-I', str(helper), str(out / 'models')], check=True, capture_output=True)
    claims = out / 'claims.jsonl'
    claims.write_text(json.dumps({'protein_id': 'pE', 'function_type': 'EC', 'function_id': 'EC:2.3.3.3',
                                 'source': 'diamond', 'evidence_type': 'SEQUENCE_SIMILARITY',
                                 'raw_score': .6, 'calibrated_prob': .6}) + '\n')
    template = json.loads((out / 'models/template.json').read_text())
    template['model'] = str(out / 'models/toy.xml')
    template['absence_comparison'] = {
        'reaction': 'VITS', 'candidate_function_key': 'pE|EC|EC:2.3.3.3',
        'assembly': {'status': 'COMPLETE', 'completeness_fraction': 1.0, 'source': 'synthetic closed genome'},
        'function_search': {'status': 'COMPLETE', 'source': 'synthetic function search'},
        'external_supply': {'exchange': 'EX_vit_e', 'uptake_limit': 2, 'transport_reactions': ['T_vit'], 'source': 'synthetic medium'},
        'conditions': [{'id': 'control', 'mods': []},
                       {'id': 'withdrawal', 'mods': [{'op': 'close', 'exchange': 'EX_vit_e'}]},
                       {'id': 'transport_ko', 'mods': [{'op': 'knockout', 'reactions': ['T_vit']}]}]}
    commands, cases = [], []

    def execute(name, args):
        cmd = ['java', '-jar', str(jar), 'integrate', '--claims', str(claims), '--out', str(out / (name + '.tsv')),
               '--lite', '--enable-priors', 'none', '--feature-export-dir', str(out / (name + '.features')), *args]
        completed = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
        (out / (name + '.log')).write_text(completed.stdout + completed.stderr)
        commands.append({'name': name, 'argv': cmd, 'returncode': completed.returncode})
        (out / 'commands.json').write_text(json.dumps(commands, indent=2) + '\n')
        if completed.returncode:
            raise ValueError('production CLI failed: ' + name)

    execute('baseline', [])
    for name in ('complete', 'fragmented', 'unknown', 'no_supply', 'supplied_medium', 'low_supply', 'budget', 'transport_disabled'):
        manifest = copy.deepcopy(template)
        request = manifest['absence_comparison']
        budget = 24
        if name == 'fragmented':
            request['candidate_function_key'] = None
            request['assembly'].update(status='DRAFT', completeness_fraction=.5)
            request['function_search']['status'] = 'INCOMPLETE'
        elif name == 'unknown':
            request['assembly'].update(status='UNKNOWN', completeness_fraction=None)
            request['function_search']['status'] = 'UNKNOWN'
        elif name == 'no_supply':
            request['external_supply'] = None
        elif name == 'supplied_medium':
            manifest['baseline']['medium']['EX_vit_e'] = 2
        elif name == 'low_supply':
            request['external_supply']['uptake_limit'] = .03
        elif name == 'budget':
            budget = 7
        elif name == 'transport_disabled':
            request['conditions'][0]['mods'] = [{'op': 'knockout', 'reactions': ['T_vit']}]
        path = out / (name + '.json')
        path.write_text(json.dumps(manifest, indent=2) + '\n')
        execute(name, ['--experiment-manifest', str(path), '--experiment-mode', 'manifest', '--experiment-strategy', 'fixed',
                       '--experiment-python-bin', sys.executable, '--experiment-sidecar', str(ROOT / 'scripts/metabolic_experiments.py'),
                       '--experiment-max-experiments', '3', '--experiment-max-solver-calls', str(budget)])
        result_path = out / (name + '.tsv.experiments/result.json')
        result = json.loads(result_path.read_text())
        baseline = (out / 'baseline.tsv').read_bytes()
        if (out / (name + '.tsv')).read_bytes() != baseline:
            raise ValueError('conditional analysis changed predictions: ' + name)
        features = (out / (name + '.features/candidate_features.jsonl')).read_bytes()
        if features != (out / 'baseline.features/candidate_features.jsonl').read_bytes():
            raise ValueError('absence analysis changed learned features: ' + name)
        oracle = compare(path, result_path)
        (out / (name + '.cobra.json')).write_text(json.dumps(oracle, indent=2) + '\n')
        if name == 'budget':
            assert result['cells'] == [] and not oracle['comparisons']
        elif not oracle['agrees']:
            raise ValueError('independent LP disagreement: ' + name)
        cases.append({'name': name, 'solver_calls': len(result['solver_calls']), 'cells': len(result['cells']),
                      'oracle_comparisons': len(oracle['comparisons']), 'oracle_agrees': oracle['agrees'] if oracle['comparisons'] else None,
                      'prediction_and_feature_bytes_unchanged': True, 'gene_attribution': result['gene_attribution']})
        print(json.dumps(cases[-1]), flush=True)
    inputs = [jar, claims, helper, ROOT / 'scripts/metabolic_experiments.py', ROOT / 'scripts/absence_experiments.py',
              ROOT / 'benchmark/v3/check_absence_cobra.py', Path(__file__).resolve()]
    report = {'schema': 'gspa.absence_mechanism_replay/1', 'status': 'verified',
              'interpretation': 'Synthetic mechanism and independent implementation checks, not biological validation or sequence-fragmentation accuracy.',
              'inputs': {str(p): sha(p) for p in inputs}, 'cases': cases,
              'outputs': {str(p.relative_to(out)): sha(p) for p in sorted(out.rglob('*')) if p.is_file()}}
    (out / 'verification.json').write_text(json.dumps(report, indent=2) + '\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--jar', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    run(args.jar.resolve(), args.out.resolve())
