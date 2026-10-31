#!/usr/bin/env python3
"""Negative sulfur-supply controls; no phenotype file is read or scored."""
import argparse
import copy
import json
from pathlib import Path

import run_sulfonate_benchmark as S
import check_nutrient_cobra as O
import cobra


def check(plan, model_path):
    frozen, plans = S.verify(plan, model_path)
    original = next(p for p in plans if p['sulfur_source'] == 'Sulfate'
                    and p['genotype_label'] == 'None' and p['arm'] == 'sulfur_uptake=1')
    model = S.ME.load_sbml_model(model_path)
    other = cobra.io.read_sbml_model(model_path, f_replace={})
    other.solver = 'glpk'
    results = []
    for arm in ('no_sulfur', 'no_sulfur_no_thiamine'):
        entry = copy.deepcopy(original)
        entry['arm'] = arm
        medium = entry['manifest']['baseline']['medium']
        medium['R_EX_so4_e'] = 0.
        if arm == 'no_sulfur_no_thiamine':
            medium['R_EX_thm_e'] = 0.
        production = S.N.simulate(entry, model)
        record = production['production']['experiments'][0]
        oracle = O.simulate(other, entry['manifest'], [])
        results.append({'arm': arm, 'plan': entry, 'production': production,
                        'oracle': oracle, 'agrees': O.agreement(record, oracle),
                        'growth_detected': production['prediction']})
    return {'schema': 'gspa.sulfonate_medium_controls/1',
            'protocol_sha256': S.N.sha(plan / 'protocol.json'),
            'implementation_sha256': S.N.sha(Path(__file__)),
            'model_sha256': frozen['model_sha256'], 'results': results,
            'interpretation': 'Model controls, not study phenotypes. Growth would invalidate the claimed sulfur restriction; unknown is not a passing negative control.',
            'passes': all(r['agrees'] and r['growth_detected'] is False for r in results)}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('plan', 'model', 'out'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error('output exists')
    report = check(args.plan, args.model)
    S.N.write_json(args.out, report)
    print(json.dumps({'passes': report['passes'], 'controls': len(report['results'])}))
    raise SystemExit(0 if report['passes'] else 1)
