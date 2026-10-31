#!/usr/bin/env python3
"""Independent SBML/GPR/GLPK replay of every registered nutrient phenotype LP."""
import argparse
import json
import math
from pathlib import Path
import cobra
from run_nutrient_phenotypes import read_rows, sha, digest, write_json


def simulate(model, manifest, genes):
    with model:
        medium = manifest['baseline']['medium']
        for reaction in model.reactions:
            if len(reaction.metabolites) == 1:
                coefficient = next(iter(reaction.metabolites.values()))
                low, high = reaction.bounds
                if coefficient < 0:
                    low = -medium[reaction.id] if reaction.id in medium else max(low, 0.)
                else:
                    high = medium[reaction.id] if reaction.id in medium else min(high, 0.)
                if low > high:
                    return {'status': 'infeasible', 'growth': None, 'blocked_reactions': None}
                reaction.bounds = low, high
        maintenance = manifest['objectives']['maintenance']
        if maintenance != 'none':
            reaction = model.reactions.get_by_id(maintenance['reaction'])
            required = max(reaction.lower_bound, maintenance['lower_bound'])
            if required > reaction.upper_bound:
                return {'status': 'infeasible', 'growth': None, 'blocked_reactions': None}
            reaction.lower_bound = required
        for gene in genes:
            model.genes.get_by_id(gene).knock_out()
        blocked = sorted(r.id for r in model.reactions if r.gene_reaction_rule and not r.functional)
        model.objective = model.reactions.get_by_id(manifest['objectives']['biomass'])
        model.objective_direction = 'max'
        growth = float(model.slim_optimize())
        return {'status': model.solver.status, 'growth': growth if math.isfinite(growth) else None,
                'blocked_reactions': blocked}


def agreement(record, oracle, tolerance=1e-6):
    if record['status'] != 'completed':
        return False
    if record['outcome'] == 'infeasible_with_maintenance':
        return oracle['status'] == 'infeasible'
    value = record.get('measurements', {}).get('max_growth')
    return (oracle['status'] == 'optimal' and value is not None and
            math.isclose(value, oracle['growth'], abs_tol=tolerance, rel_tol=tolerance))


def compare(plan_root, run_root, model_paths):
    frozen = json.loads((plan_root / 'protocol.json').read_text())
    receipt = json.loads((run_root / 'run.json').read_text())
    if (sha(plan_root / 'plans.jsonl') != frozen['plans_sha256'] or
            receipt['plans_sha256'] != frozen['plans_sha256'] or
            receipt['protocol_sha256'] != sha(plan_root / 'protocol.json') or
            receipt['predictions_sha256'] != sha(run_root / 'predictions.jsonl')):
        raise ValueError('run or plan receipt mismatch')
    plans = {(r['case_id'], r['arm']): r for r in read_rows(plan_root / 'plans.jsonl')}
    predictions = read_rows(run_root / 'predictions.jsonl')
    if len(predictions) != len(plans) or {(r['case_id'], r['arm']) for r in predictions} != set(plans):
        raise ValueError('missing or duplicated prediction records')
    models = {}
    for sp, path in model_paths.items():
        if sha(path) != frozen['models'][sp]['sha256']:
            raise ValueError('oracle model differs from frozen model')
        models[sp] = cobra.io.read_sbml_model(path, f_replace={})
        models[sp].solver = 'glpk'
    results = []
    for i, row in enumerate(predictions):
        plan = plans[row['case_id'], row['arm']]
        if digest(plan) != row['plan_sha256']:
            raise ValueError('prediction plan changed')
        result = {'case_id': row['case_id'], 'arm': row['arm'], 'comparisons': []}
        for record in row.get('production', {}).get('experiments', []):
            genes = [] if record['kind'] == 'baseline' else plan['genes']
            oracle = simulate(models[plan['species']], plan['manifest'], genes)
            blocked = record['measurements'].get('knocked_out_reactions')
            result['comparisons'].append({'kind': record['kind'], 'oracle': oracle,
                'flux_agrees': agreement(record, oracle),
                'gpr_agrees': oracle['blocked_reactions'] is not None and blocked == oracle['blocked_reactions']})
        result['agrees'] = (len(result['comparisons']) == (2 if plan['genes'] else 1) and
                            all(c['flux_agrees'] and c['gpr_agrees'] for c in result['comparisons']))
        results.append(result)
        if (i + 1) % 100 == 0:
            print(json.dumps({'oracle_completed': i + 1}), flush=True)
    return {'schema': 'gspa.nutrient_phenotype_oracle/1', 'implementation_sha256': sha(Path(__file__)),
            'run_sha256': sha(run_root / 'run.json'), 'cobra_version': cobra.__version__,
            'solver': 'GLPK; independently parsed SBML and GPRs', 'results': results,
            'comparisons': sum(len(r['comparisons']) for r in results),
            'agrees': all(r['agrees'] for r in results),
            'disagreements': sum(not r['agrees'] for r in results)}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('plan', 'run', 'ec-model', 'bs-model', 'out'):
        p.add_argument('--' + key, type=Path, required=True)
    args = p.parse_args()
    if args.out.exists():
        p.error('oracle output exists')
    report = compare(args.plan, args.run, {'Ec': args.ec_model, 'Bs': args.bs_model})
    write_json(args.out, report)
    print(json.dumps({k: report[k] for k in ('comparisons', 'agrees', 'disagreements')}))
    raise SystemExit(0 if report['agrees'] else 1)
