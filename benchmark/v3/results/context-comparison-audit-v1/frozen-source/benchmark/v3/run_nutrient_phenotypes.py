#!/usr/bin/env python3
"""Freeze, execute and score known-genotype nutrient phenotypes in separate stages.

Planning and simulation never read measurements. The production experiment
engine supplies the predictions; endpoint OD supplies phenotype observations,
not gene-function truth. All predeclared arms and all failed cases are retained.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
import metabolic_experiments as ME

PROTOCOL = {
    'schema': 'gspa.nutrient_phenotype_protocol/1',
    'scope': 'All single-organism assays; no-nitrogen controls reported separately; no outcome filtering.',
    'oxygen_rates': [0., 2., 20.],
    'medium_arms': ['printed', 'additional_inorganic'],
    'arm_selection': 'All six reported; no arm selected using outcomes.',
    'rates_mmol_gdw_h': {'named_carbon': 10., 'named_nitrogen': 10.,
                        'inorganic': 1000., 'tryptophan': .1, 'thiamine': .01,
                        'MOPS': 10., 'supplement': 1., 'serine_supplement': 5.},
    'rate_interpretation': 'Fixed modeling assumptions, not converted from concentrations or fitted to OD.',
    'additional_inorganic': ['so4', 'fe2', 'cu2', 'cobalt2', 'mobd', 'ni2'],
    'additional_inorganic_note': 'Unprinted supply sensitivity, including sulfate; not the measured recipe.',
    'buffer': 'Supply printed MOPS if the exact exchange exists; otherwise explicitly unrepresented.',
    'missing_components': 'Record every missing exact exchange; do not substitute related molecules.',
    'maintenance': 'Model-native ATPM lower bound; never tuned or disabled.',
    'genotype': 'Exact symbol to canonical gene mapping; production GPR knockout, including isozyme bypasses.',
    'growth_threshold': 1e-6,
    'observed_growth': 'Raw-reconstructed mean OD600 > 0.02; endpoint presence, not growth rate.',
    'paired_rescue': 'Supplemented growth and unsupplemented no-growth, at the same strain/environment/arm.',
    'interpretation': 'Model-conditional phenotype feasibility; no independent gene attribution or calibrated probability.',
}
CARBON = {None: [], 'D-glucose': ['glc__D'], 'glycerol': ['glyc'], 'sodium fumarate': ['fum', 'na1']}
NITROGEN = {None: [], 'L-glutamic acid': ['glu__L'], 'ammonium sulfate': ['nh4', 'so4'],
            'urea': ['urea'], 'L-aspartic acid': ['asp__L'], 'L-asparagine': ['asn__L'],
            'L-glutamine': ['gln__L'], 'N-acetylglucosamine': ['acgam'],
            'D-glucosamine hydrochloride': ['gam', 'cl'], 'L-alanine': ['ala__L'], 'L-proline': ['pro__L']}
SUPPLEMENT = {'methionine': 'met__L', 'phenylalanine': 'phe__L', 'arginine': 'arg__L',
              'cysteine': 'cys__L', 'serine': 'ser__L', 'leucine': 'leu__L'}
BASE = {'potassium phosphate': ['k', 'pi'], 'MgCl2': ['mg2', 'cl'], 'MnCl2': ['mn2', 'cl'],
        'ZnCl2': ['zn2', 'cl'], 'CaCl2': ['ca2', 'cl'], 'FeCl3': ['fe3', 'cl'],
        'L-tryptophan': ['trp__L'], 'thiamine': ['thm'], 'MOPS': ['mops']}


def sha(path):
    hasher = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            hasher.update(chunk)
    return hasher.hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def read_rows(path):
    return [json.loads(line) for line in Path(path).read_text().splitlines() if line.strip()]


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def code_hashes():
    return {str(p.relative_to(ROOT)): sha(p) for p in [Path(__file__), ROOT / 'scripts/metabolic_experiments.py']}


def medium_for(case, model, oxygen, arm):
    if oxygen not in PROTOCOL['oxygen_rates'] or arm not in PROTOCOL['medium_arms']:
        raise ValueError('unregistered sensitivity arm')
    rates, supplied, bindings = PROTOCOL['rates_mmol_gdw_h'], {}, []

    def add(component, molecules, rate):
        for molecule in molecules:
            rid = 'R_EX_' + molecule + '_e'
            present = rid in model.rxn_index
            bindings.append({'component': component, 'exchange': rid,
                             'status': 'represented' if present else 'unrepresented', 'assumed_rate': rate})
            if present:
                # Exact extracellular chemistry is checked by the production validator.
                supplied[rid] = max(supplied.get(rid, 0.), rate)

    raw = case['medium']
    if set(raw['base']) != set(BASE):
        raise ValueError('base medium differs from registered ingredient set')
    for component, molecules in BASE.items():
        rate_key = {'L-tryptophan': 'tryptophan', 'thiamine': 'thiamine', 'MOPS': 'MOPS'}.get(component, 'inorganic')
        add(component, molecules, rates[rate_key])
    add('named carbon', CARBON[raw['carbon_source']], rates['named_carbon'])
    add('named nitrogen', NITROGEN[raw['nitrogen_source']], rates['named_nitrogen'])
    add('water/proton balance', ['h2o', 'h'], rates['inorganic'])
    add('assumed oxygen', ['o2'], oxygen)
    if raw['supplement']:
        add('measured supplement', [SUPPLEMENT[raw['supplement']]],
            rates['serine_supplement' if raw['supplement'] == 'serine' else 'supplement'])
    if arm == 'additional_inorganic':
        add('unprinted inorganic sensitivity', PROTOCOL['additional_inorganic'], rates['inorganic'])
    return dict(sorted(supplied.items())), bindings


def make_plan(case, model, registration, oxygen, arm):
    if len(case['organisms']) != 1 or case['kind'] not in ('monoculture', 'supplemented_monoculture', 'wild_type'):
        raise ValueError('not a single-organism phenotype assay')
    mappings = {g['symbol']: g for g in registration['genes']}
    genes = []
    for symbol in case['organisms'][0]['deleted_genes']:
        found = mappings.get(symbol, {}).get('canonical')
        if not found or model.gene_alias.get(found) != found:
            raise ValueError('unmapped deletion: ' + symbol)
        genes.append(found)
    medium, bindings = medium_for(case, model, oxygen, arm)
    if len(model.objective_terms) != 1 or model.objective_terms[0][1] <= 0 or model.objective_sense != 'maximize':
        raise ValueError('requires a single maximized biomass objective')
    manifest = {'model': registration['path'], 'model_sha256_expected': registration['sha256'],
                'objectives': {'biomass': model.objective_terms[0][0],
                               'maintenance': {'reaction': 'R_ATPM', 'lower_bound': float(model.lb[model.rxn_index['R_ATPM']])}},
                'baseline': {'exchanges': list(medium), 'medium': medium},
                'context': {'hypotheses': [{'id': 'observed_genotype', 'genes': genes}] if genes else []},
                'config': {'strategy': 'fixed', 'max_experiments': 1, 'max_solver_calls': 2,
                           'growth_threshold': PROTOCOL['growth_threshold'],
                           'max_focus_reactions': 0, 'max_baseline_focus_reactions': 0}}
    ME.validate_manifest(manifest)
    ME.validate_boundary(model, manifest)
    return {'case_id': case['case_id'], 'species': case['species'], 'kind': case['kind'],
            'environment': case['environment'], 'strain_id': case['organisms'][0]['strain_id'],
            'primary_environment': case['primary_environment'], 'genes': genes,
            'arm': f'{arm}:oxygen={oxygen:g}', 'bindings': bindings, 'manifest': manifest}


def prepare(inputs, registration_path, out):
    registration = json.loads(registration_path.read_text())
    if sha(inputs) != registration['inputs_sha256']:
        raise ValueError('inputs do not match model registration')
    records = read_rows(inputs)
    if len({r['case_id'] for r in records}) != len(records):
        raise ValueError('duplicate input cases')
    models = {}
    for species, meta in registration['models'].items():
        if sha(meta['path']) != meta['sha256']:
            raise ValueError('registered model changed')
        models[species] = ME.load_sbml_model(Path(meta['path']))
    selected = [r for r in records if r['kind'] != 'coculture']
    out.mkdir(parents=True, exist_ok=False)
    with (out / 'plans.jsonl').open('w') as stream:
        for case in selected:
            for arm in PROTOCOL['medium_arms']:
                for oxygen in PROTOCOL['oxygen_rates']:
                    row = make_plan(case, models[case['species']], registration['models'][case['species']], oxygen, arm)
                    stream.write(json.dumps(row, sort_keys=True, allow_nan=False) + '\n')
    frozen = {'schema': 'gspa.nutrient_phenotype_plan/1', 'protocol': PROTOCOL,
              'inputs_sha256': sha(inputs), 'registration_sha256': sha(registration_path),
              'models': registration['models'], 'plans_sha256': sha(out / 'plans.jsonl'),
              'implementation': code_hashes(), 'cases': len(selected), 'plans': len(selected) * 6,
              'excluded': [{'case_id': r['case_id'], 'reason': 'separate community evaluation'} for r in records if r['kind'] == 'coculture']}
    write_json(out / 'protocol.json', frozen)
    return frozen


def prediction(record):
    if record['status'] != 'completed':
        return None
    if record['outcome'] == 'infeasible_with_maintenance':
        return False
    growth = record.get('measurements', {}).get('max_growth')
    if growth is None or not math.isfinite(growth):
        return None
    return growth >= PROTOCOL['growth_threshold']


def simulate(plan, model):
    manifest = plan['manifest']
    if model.sha256 != manifest['model_sha256_expected']:
        raise ValueError('execution model hash mismatch')
    if plan['genes']:
        result = ME.run_experiments(manifest, model)
        records = result['experiments']
        if [r['kind'] for r in records] != ['baseline', 'knockout']:
            raise ValueError('production engine did not execute the registered WT/KO pair')
        target = records[1]
        entrypoint = 'metabolic_experiments.run_experiments'
    else:
        # The public experiment selector intentionally does not solve empty menus.
        # A measured WT control needs a baseline LP, not an invented hypothesis.
        runner = ME.Runner(model, manifest, ME.resolve_config(manifest))
        target = runner.run(ME.baseline_item())
        result = {'experiments': [target], 'solver_log': runner.solver.log}
        entrypoint = 'metabolic_experiments.Runner.run(baseline_item())'
    return {'case_id': plan['case_id'], 'arm': plan['arm'], 'plan_sha256': digest(plan),
            'entrypoint': entrypoint, 'prediction': prediction(target), 'production': result}


def verify_plan(root):
    frozen = json.loads((root / 'protocol.json').read_text())
    if frozen['protocol'] != PROTOCOL or frozen['implementation'] != code_hashes():
        raise ValueError('frozen protocol or execution source changed')
    if frozen['plans_sha256'] != sha(root / 'plans.jsonl'):
        raise ValueError('frozen plans changed')
    return frozen


def run(root, out, model_paths):
    frozen = verify_plan(root)
    models = {sp: ME.load_sbml_model(path) for sp, path in model_paths.items()}
    for sp, model in models.items():
        if model.sha256 != frozen['models'][sp]['sha256']:
            raise ValueError('execution model differs from registration')
    out.mkdir(parents=True, exist_ok=False)
    start = time.time()
    counts = Counter()
    with (out / 'predictions.jsonl').open('w') as stream:
        for i, plan in enumerate(read_rows(root / 'plans.jsonl')):
            try:
                row = simulate(plan, models[plan['species']])
            except Exception as exc:
                row = {'case_id': plan['case_id'], 'arm': plan['arm'], 'plan_sha256': digest(plan),
                       'prediction': None, 'failure': f'{type(exc).__name__}: {exc}'}
            counts[str(row['prediction'])] += 1
            stream.write(json.dumps(row, sort_keys=True, allow_nan=False) + '\n')
            if (i + 1) % 100 == 0:
                stream.flush()
                print(json.dumps({'completed': i + 1, 'elapsed_s': round(time.time() - start, 1)}), flush=True)
    report = {'schema': 'gspa.nutrient_phenotype_run/1', 'protocol_sha256': sha(root / 'protocol.json'),
              'plans_sha256': sha(root / 'plans.jsonl'), 'predictions_sha256': sha(out / 'predictions.jsonl'),
              'counts': dict(counts), 'elapsed_s': time.time() - start,
              'python': sys.version, 'platform': platform.platform(), 'implementation': code_hashes()}
    write_json(out / 'run.json', report)
    return report


def confusion(pairs):
    count = Counter()
    for observed, predicted in pairs:
        if observed is None or predicted is None:
            count['unknown'] += 1
        else:
            count[('t' if observed == predicted else 'f') + ('p' if predicted else 'n')] += 1
    result = {k: count[k] for k in ('tp', 'tn', 'fp', 'fn', 'unknown')}
    n = sum(count[k] for k in ('tp', 'tn', 'fp', 'fn'))
    result.update(n=len(pairs), accuracy=(count['tp'] + count['tn']) / n if n else None)
    return result


def evaluate(root, run_root, observations, out):
    frozen = verify_plan(root)
    receipt = json.loads((run_root / 'run.json').read_text())
    if (receipt['protocol_sha256'] != sha(root / 'protocol.json') or
            receipt['plans_sha256'] != frozen['plans_sha256'] or
            receipt['predictions_sha256'] != sha(run_root / 'predictions.jsonl')):
        raise ValueError('prediction provenance mismatch')
    plans = {(r['case_id'], r['arm']): r for r in read_rows(root / 'plans.jsonl')}
    rows = read_rows(run_root / 'predictions.jsonl')
    if len(rows) != len(plans) or {(r['case_id'], r['arm']) for r in rows} != set(plans):
        raise ValueError('incomplete or duplicated prediction catalogue')
    obs_rows = read_rows(observations)
    obs = {r['case_id']: r for r in obs_rows}
    if len(obs) != len(obs_rows):
        raise ValueError('duplicate observations')
    strata, paired, scored = defaultdict(list), defaultdict(dict), []
    for row in rows:
        plan = plans[row['case_id'], row['arm']]
        if digest(plan) != row['plan_sha256']:
            raise ValueError('prediction plan hash mismatch')
        outcome = obs.get(row['case_id'], {})
        mean = outcome.get('mean_od600')
        if mean is not None and (isinstance(mean, bool) or not isinstance(mean, (int, float)) or not math.isfinite(mean)):
            raise ValueError('nonfinite or malformed observed OD')
        observed = None if mean is None else mean > .02
        predicted = row['prediction']
        if predicted is not None and not isinstance(predicted, bool):
            raise ValueError('prediction must be boolean or unknown')
        scope = 'primary' if plan['primary_environment'] else 'no_nitrogen_control'
        group = (row['arm'], scope, plan['species'])
        strata[group + ('all',)].append((observed, predicted))
        strata[group + (plan['kind'],)].append((observed, predicted))
        strata[group + ('strain:' + plan['strain_id'],)].append((observed, predicted))
        strata[group + ('environment:' + plan['environment'],)].append((observed, predicted))
        paired[group + (plan['strain_id'], plan['environment'])][plan['kind']] = (observed, predicted)
        scored.append({'case_id': row['case_id'], 'arm': row['arm'], 'observed_growth': observed,
                       'predicted_growth': predicted, 'mean_od600': mean,
                       'experiment_mean_sem': outcome.get('experiment_mean_sem')})
    rescues = defaultdict(list)
    for key, assays in paired.items():
        if 'wild_type' in assays:
            continue
        no, yes = assays.get('monoculture', (None, None)), assays.get('supplemented_monoculture', (None, None))
        outcome = [None if no[i] is None or yes[i] is None else yes[i] and not no[i] for i in (0, 1)]
        rescues[key[:3]].append(outcome)
    report = {'schema': 'gspa.nutrient_phenotype_evaluation/1',
              'protocol_sha256': sha(root / 'protocol.json'), 'run_sha256': sha(run_root / 'run.json'),
              'observations_sha256': sha(observations), 'implementation': code_hashes(),
              'strata': [{'arm': k[0], 'scope': k[1], 'species': k[2], 'stratum': k[3], **confusion(v)} for k, v in sorted(strata.items())],
              'paired_rescue': [{'arm': k[0], 'scope': k[1], 'species': k[2], **confusion(v)} for k, v in sorted(rescues.items())],
              'cases': scored, 'interpretation': PROTOCOL['interpretation']}
    if out.exists():
        raise ValueError('evaluation output exists')
    write_json(out, report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    prep = sub.add_parser('prepare')
    for key in ('inputs', 'registration', 'out'):
        prep.add_argument('--' + key, type=Path, required=True)
    execute = sub.add_parser('run')
    for key in ('plan', 'out', 'ec-model', 'bs-model'):
        execute.add_argument('--' + key, type=Path, required=True)
    score = sub.add_parser('evaluate')
    for key in ('plan', 'run', 'observations', 'out'):
        score.add_argument('--' + key, type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'prepare':
        report = prepare(args.inputs, args.registration, args.out)
        print(json.dumps({k: report[k] for k in ('cases', 'plans', 'plans_sha256')}))
    elif args.command == 'run':
        print(json.dumps(run(args.plan, args.out, {'Ec': args.ec_model, 'Bs': args.bs_model})))
    else:
        report = evaluate(args.plan, args.run, args.observations, args.out)
        print(json.dumps({'strata': len(report['strata']), 'paired_rescue': report['paired_rescue']}))


if __name__ == '__main__':
    main()
