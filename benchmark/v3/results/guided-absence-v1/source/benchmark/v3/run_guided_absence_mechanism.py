#!/usr/bin/env python3
"""Replay automatic gap-to-absence comparisons through the CLI and independent GLPK."""
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
    model = out / 'models/toy.xml'
    unannotated = model.read_text()
    import re
    model.write_text(re.sub(r'(<reaction\b[^>]*id="VITS"[^>]*>)', r'\1' +
        '<annotation><rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" '
        'xmlns:bqbiol="http://biomodels.net/biology-qualifiers/"><rdf:Description rdf:about="#VITS">'
        '<bqbiol:is><rdf:Bag><rdf:li rdf:resource="https://identifiers.org/ec-code/2.3.3.3"/>'
        '</rdf:Bag></bqbiol:is></rdf:Description></rdf:RDF></annotation>', unannotated))
    table_only = out / 'models/table-only.xml'
    table_only.write_text(unannotated)
    gaps = out / 'gaps.jsonl'
    gaps.write_text(json.dumps({'pathway_id': 'synthetic:P', 'reaction_id': 'source:VITS', 'ec_number': '2.3.3.3'}) + '\n')
    claims = out / 'claims.jsonl'
    claims.write_text(''.join(json.dumps({'protein_id': p, 'function_type': 'EC', 'function_id': 'EC:2.3.3.3',
        'source': 'diamond', 'evidence_type': 'SEQUENCE_SIMILARITY', 'raw_score': .6, 'calibrated_prob': .6}) + '\n'
        for p in ('pE', 'pG')))
    template = json.loads((out / 'models/template.json').read_text())
    template['model'] = str(model)
    template['absence_comparison'] = {
        'assembly': {'status': 'UNKNOWN', 'completeness_fraction': None, 'source': 'synthetic test; not measured'},
        'function_search': {'status': 'UNKNOWN', 'source': 'synthetic test; not measured'},
        'external_supply_by_reaction': {'VITS': {'exchange': 'EX_vit_e', 'uptake_limit': 2,
            'transport_reactions': ['T_vit'], 'source': 'synthetic medium'}},
        'conditions': [{'id': 'control', 'mods': []},
                       {'id': 'withdrawal', 'mods': [{'op': 'close', 'exchange': 'EX_vit_e'}]},
                       {'id': 'transport_ko', 'mods': [{'op': 'knockout', 'reactions': ['T_vit']}]}]}
    commands, cases = [], []

    def execute(name, args):
        cmd = ['java', '-jar', str(jar), 'integrate', '--claims', str(claims), '--out', str(out / (name + '.tsv')),
               '--lite', '--enable-priors', 'none', '--gaps', str(gaps),
               '--feature-export-dir', str(out / (name + '.features')), *args]
        completed = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
        (out / (name + '.log')).write_text(completed.stdout + completed.stderr)
        commands.append({'name': name, 'argv': cmd, 'returncode': completed.returncode})
        (out / 'commands.json').write_text(json.dumps(commands, indent=2) + '\n')
        if completed.returncode:
            raise ValueError('production CLI failed: ' + name)

    execute('baseline', [])
    for name, budget, max_conditions, expected_calls in (
        ('complete', 48, 6, 48), ('solver_budget', 8, 6, 8), ('condition_budget', 48, 1, 8),
        ('no_supply', 48, 6, 36), ('table_only', 48, 6, 0)):
        manifest = copy.deepcopy(template)
        if name == 'no_supply':
            manifest['absence_comparison'].pop('external_supply_by_reaction')
        elif name == 'table_only':
            manifest['model'] = str(table_only)
        path = out / (name + '.json')
        path.write_text(json.dumps(manifest, indent=2) + '\n')
        execute(name, ['--experiment-manifest', str(path), '--experiment-mode', 'evidence_guided',
            '--experiment-strategy', 'fixed', '--experiment-python-bin', sys.executable,
            '--experiment-sidecar', str(ROOT / 'scripts/metabolic_experiments.py'),
            '--experiment-max-experiments', str(max_conditions), '--experiment-max-solver-calls', str(budget)])
        result = json.loads((out / (name + '.tsv.experiments/result.json')).read_text())
        assert result['solver_calls_used'] == expected_calls
        assert result['conditions_visited'] <= max_conditions
        assert result['selected_explanation'] is None and result['gene_attribution'] == 'unverified'
        assert (out / (name + '.tsv')).read_bytes() == (out / 'baseline.tsv').read_bytes()
        assert (out / (name + '.features/candidate_features.jsonl')).read_bytes() == (out / 'baseline.features/candidate_features.jsonl').read_bytes()
        oracles = []
        for comparison in result['comparisons']:
            child_result = Path(comparison['result_file'])
            oracle = compare(child_result.parent / 'manifest.json', child_result)
            if oracle['comparisons']:
                assert oracle['agrees'], (name, comparison['id'])
            oracles.append(oracle)
        (out / (name + '.cobra.json')).write_text(json.dumps(oracles, indent=2) + '\n')
        cases.append({'name': name, 'solver_calls': expected_calls, 'conditions_visited': result['conditions_visited'],
            'oracle_comparisons': sum(len(o['comparisons']) for o in oracles),
            'verified_comparisons': len(result['comparisons']), 'unvisited_comparisons': len(result['unvisited']),
            'predictions_and_features_unchanged': True})
        print(json.dumps(cases[-1]), flush=True)
    report = {'schema': 'gspa.guided_absence_mechanism_replay/1', 'status': 'verified',
        'interpretation': 'Synthetic execution and independent LP verification; not biological absence or gene-function validation.',
        'inputs': {str(p): sha(p) for p in (jar, helper, Path(__file__).resolve(), ROOT / 'benchmark/v3/check_absence_cobra.py',
            ROOT / 'scripts/metabolic_experiments.py', ROOT / 'scripts/absence_experiments.py')},
        'cases': cases, 'outputs': {str(p.relative_to(out)): sha(p) for p in sorted(out.rglob('*')) if p.is_file()}}
    (out / 'verification.json').write_text(json.dumps(report, indent=2) + '\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--jar', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    run(args.jar.resolve(), args.out.resolve())
