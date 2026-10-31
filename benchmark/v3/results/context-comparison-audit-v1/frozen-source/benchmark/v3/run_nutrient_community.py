#!/usr/bin/env python3
"""Outcome-blind matched coculture execution using the production community engine."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import copy
import json
import itertools
import math
from pathlib import Path
import platform
import sys
import time

import run_nutrient_phenotypes as N

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
import community_metabolism as C

PROTOCOL = {
    'schema': 'gspa.nutrient_community_protocol/1',
    'scope': 'Every two-strain coculture, including no-added-nitrogen controls; no outcome filtering.',
    'medium_protocol': N.PROTOCOL,
    'biomass_ratio': [1., 1.],
    'ratio_interpretation': 'Equal inoculum OD is modeled as equal fixed biomass fractions; neither equal biomass nor balanced growth was measured.',
    'port_ceiling': 1000., 'pool_removal_ceiling': 1000.,
    'ports': 'Exact extracellular model species; all exchange ports predicted, not curated transporter evidence. Internal transport bounds unchanged.',
    'policies': ['permissive', 'native_unsupplied', 'no_unsupplied_uptake'],
    'permissive': 'Open each exchange uptake to at least 1000; finite pool supply and partner secretion still constrain uptake.',
    'native_unsupplied': 'Keep native uptake bounds for compounds with zero environment supply; supplied compounds retain open ports.',
    'no_unsupplied_uptake': 'Close uptake of every compound with zero environment supply; sharing supplied compounds remains possible.',
    'max_solver_calls': 5,
    'solve_order': ['monoculture_0', 'monoculture_1', 'permissive', 'native_unsupplied', 'no_unsupplied_uptake'],
    'growth_threshold': 1e-6,
    'observed_endpoint': 'Raw-reconstructed bulk mean OD600 > 0.02; does not establish growth of both members.',
    'paired_rescue_proxy': 'Coculture endpoint above threshold and both matched monoculture endpoints at or below threshold.',
    'interpretation': 'Conditional balanced-growth feasibility against a bulk endpoint proxy; no measured transfer, gene attribution, calibrated probability or promotion.',
}


def code_hashes():
    return {str(p.relative_to(ROOT)): N.sha(p) for p in [Path(__file__), Path(N.__file__),
        ROOT / 'scripts/metabolic_experiments.py', ROOT / 'scripts/community_metabolism.py']}


def context_for(model, registration):
    ext = sorted(r for r in model.rxn_ids if model.is_single_metabolite(r)
                 and model.met_compartment[model.exchange_metabolite(r)] == 'e')
    metabolites = [model.exchange_metabolite(r) for r in ext]
    if len(set(metabolites)) != len(metabolites):
        raise ValueError('multiple extracellular exchanges for one species require an explicit mapping audit')
    if len(model.objective_terms) != 1 or model.objective_terms[0][1] <= 0 or model.objective_sense != 'maximize':
        raise ValueError('requires one maximized biomass objective')
    if N.sha(registration['path']) != model.sha256 or registration['sha256'] != model.sha256:
        raise ValueError('model registration mismatch')
    return {'model_sha256': model.sha256, 'biomass': model.objective_terms[0][0],
        'maintenance': {'reaction': 'R_ATPM', 'lower_bound': float(model.lb[model.rxn_index['R_ATPM']])},
        'pool_metabolites': {m: {'formula': model.met_formula.get(m)} for m in metabolites},
        'ports': {r: {'pool': m, 'metabolite': m, 'support': 'predicted',
                      'evidence': 'supplied model exchange and native internal transport; not experimental transport evidence',
                      'open_uptake': PROTOCOL['port_ceiling']} for r, m in zip(ext, metabolites)}}


def make_plan(case, model, registration, context, oxygen, arm):
    organisms = case['organisms']
    if case['kind'] != 'coculture' or len(organisms) != 2 or len({o['strain_id'] for o in organisms}) != 2:
        raise ValueError('requires two distinct coculture strains')
    if any(o['species'] != case['species'] for o in organisms) or case['assay']['coculture_od_ratio'] != [1, 1]:
        raise ValueError('unregistered species combination or inoculum ratio')
    mappings = {r['symbol']: r['canonical'] for r in registration['genes']}
    members = []
    for o in organisms:
        genes = [mappings.get(g) for g in o['deleted_genes']]
        if not genes or any(g is None or model.gene_alias.get(g) != g for g in genes):
            raise ValueError('unmapped deletion')
        members.append({'id': o['strain_id'], 'genes': genes,
                        'blocked_reactions': sorted(model.blocked_by_genes(set(genes)))})
    medium, bindings = N.medium_for(case, model, oxygen, arm)
    if set(medium) - set(context['ports']):
        raise ValueError('medium ingredient does not map to an extracellular port')
    return {'case_id': case['case_id'], 'species': case['species'], 'pair': case['pair'],
        'environment': case['environment'], 'primary_environment': case['primary_environment'],
        'arm': f'{arm}:oxygen={oxygen:g}', 'members': members,
        'supply_exchanges': medium, 'bindings': bindings}


def prepare(inputs, registration_path, out):
    registration = json.loads(registration_path.read_text())
    if N.sha(inputs) != registration['inputs_sha256']:
        raise ValueError('inputs do not match registration')
    rows = N.read_rows(inputs)
    if len({r['case_id'] for r in rows}) != len(rows):
        raise ValueError('duplicate cases')
    models = {s: N.ME.load_sbml_model(Path(m['path'])) for s, m in registration['models'].items()}
    contexts = {s: context_for(m, registration['models'][s]) for s, m in models.items()}
    selected = [r for r in rows if r['kind'] == 'coculture']
    out.mkdir(parents=True, exist_ok=False)
    N.write_json(out / 'contexts.json', contexts)
    with (out / 'plans.jsonl').open('w') as stream:
        for case in selected:
            for arm in N.PROTOCOL['medium_arms']:
                for oxygen in N.PROTOCOL['oxygen_rates']:
                    p = make_plan(case, models[case['species']], registration['models'][case['species']],
                                  contexts[case['species']], oxygen, arm)
                    stream.write(json.dumps(p, sort_keys=True, allow_nan=False) + '\n')
    frozen = {'protocol': PROTOCOL, 'implementation': code_hashes(), 'inputs_sha256': N.sha(inputs),
        'registration_sha256': N.sha(registration_path), 'models': registration['models'],
        'contexts_sha256': N.sha(out / 'contexts.json'), 'plans_sha256': N.sha(out / 'plans.jsonl'),
        'cases': len(selected), 'plans': len(selected) * 6,
        'primary_cases': sum(r['primary_environment'] for r in selected),
        'excluded': [{'case_id': r['case_id'], 'reason': 'single-organism evaluation'} for r in rows if r['kind'] != 'coculture']}
    N.write_json(out / 'protocol.json', frozen)
    return frozen


def verify_plan(root):
    frozen = json.loads((root / 'protocol.json').read_text())
    if frozen['protocol'] != PROTOCOL or frozen['implementation'] != code_hashes():
        raise ValueError('frozen protocol or implementation changed')
    for name in ('plans', 'contexts'):
        if frozen[name + '_sha256'] != N.sha(root / (name + ('.jsonl' if name == 'plans' else '.json'))):
            raise ValueError('frozen ' + name + ' changed')
    return frozen


def make_manifest(plan, context, model_path):
    supply = {context['ports'][rid]['pool']: rate for rid, rate in plan['supply_exchanges'].items()}
    specs = [{'id': m['id'], 'model': str(model_path), 'sha256': context['model_sha256'],
        'biomass': context['biomass'], 'maintenance': context['maintenance'], 'biomass_units': C.BIOMASS_UNITS,
        'flux_units': C.FLUX_UNITS, 'ports': context['ports'], 'knockouts': m['blocked_reactions']} for m in plan['members']]
    return {'organisms': specs, 'pool': {'metabolites': context['pool_metabolites'], 'supply': supply,
        'removal': {m: PROTOCOL['pool_removal_ceiling'] for m in context['pool_metabolites']}},
        'community': {'biomass_ratio': {m['id']: 1. for m in plan['members']}},
        'config': {'max_solver_calls': 5, 'auto_port_closure': False, 'max_fva_candidates': 0,
                   'growth_threshold': PROTOCOL['growth_threshold']}}


def prediction(record):
    if record['status'] != 'completed':
        return None
    if record['outcome'] == 'infeasible_with_maintenance':
        return False
    mu = record.get('mu')
    return None if mu is None or not math.isfinite(mu) else mu >= PROTOCOL['growth_threshold']


def simulate(plan, context, model, model_path):
    if model.sha256 != context['model_sha256']:
        raise ValueError('execution model differs from frozen context')
    for member in plan['members']:
        if sorted(model.blocked_by_genes(set(member['genes']))) != member['blocked_reactions']:
            raise ValueError('frozen genotype knockout differs from execution model')
    manifest = make_manifest(plan, context, model_path)
    cfg = {**C.DEFAULT_CONFIG, **manifest['config']}
    orgs = [C.build_organism(s['id'], model, s, context['pool_metabolites'], cfg) for s in manifest['organisms']]
    analysis = C.Analysis(manifest, orgs, cfg)
    for i, o in enumerate(orgs):
        analysis.growth(f'monoculture_{i}', 'monoculture', [o.id], analysis.lp([o.id]))
    analysis.growth('permissive', 'community', analysis.order, analysis.lp(analysis.order))
    native_orgs = []
    for spec in manifest['organisms']:
        native = copy.deepcopy(spec)
        for rid, port in native['ports'].items():
            if plan['supply_exchanges'].get(rid, 0.) <= 0:
                port.pop('open_uptake')
        native_orgs.append(C.build_organism(spec['id'], model, native, context['pool_metabolites'], cfg))
    native_lp = C.build_lp(native_orgs, analysis.ratio, analysis.supply, analysis.removal)
    analysis.growth('native_unsupplied', 'community', analysis.order, native_lp)
    closures = [(o.id, rid, 'uptake') for o in orgs for rid in o.ports if plan['supply_exchanges'].get(rid, 0.) <= 0]
    analysis.growth('no_unsupplied_uptake', 'community', analysis.order, analysis.lp(analysis.order, closures))
    for rec in analysis.records:
        rec.pop('_x', None)
        rec.pop('_lp', None)
    return {'case_id': plan['case_id'], 'arm': plan['arm'], 'plan_sha256': N.digest(plan),
        'manifest_sha256': N.digest(manifest), 'entrypoint': 'community_metabolism.Analysis.growth / build_lp',
        'predictions': {r['id']: prediction(r) for r in analysis.records}, 'production': analysis.records,
        'solver_calls': analysis.solver.log, 'promotion': 'none', 'gene_attribution': 'none'}


def run(root, out, paths, shard_index=0, shard_count=1):
    if shard_count < 1 or not 0 <= shard_index < shard_count:
        raise ValueError('invalid shard index/count')
    frozen = verify_plan(root)
    contexts = json.loads((root / 'contexts.json').read_text())
    models = {sp: N.ME.load_sbml_model(p) for sp, p in paths.items()}
    if any(m.sha256 != frozen['models'][sp]['sha256'] for sp, m in models.items()):
        raise ValueError('execution model differs from registration')
    out.mkdir(parents=True, exist_ok=False)
    counts, start = Counter(), time.time()
    with (out / 'predictions.jsonl').open('w') as stream:
        for i, plan in enumerate(N.read_rows(root / 'plans.jsonl')):
            if i % shard_count != shard_index:
                continue
            try:
                row = simulate(plan, contexts[plan['species']], models[plan['species']], paths[plan['species']])
            except Exception as exc:
                row = {'case_id': plan['case_id'], 'arm': plan['arm'], 'plan_sha256': N.digest(plan),
                    'predictions': {k: None for k in PROTOCOL['solve_order']}, 'failure': f'{type(exc).__name__}: {exc}'}
            counts[str(row['predictions']['permissive'])] += 1
            stream.write(json.dumps(row, sort_keys=True, allow_nan=False) + '\n')
            if (i + 1) % 100 == 0:
                stream.flush()
                print(json.dumps({'completed': i + 1, 'elapsed_s': round(time.time() - start, 1)}), flush=True)
    receipt = {'protocol_sha256': N.sha(root / 'protocol.json'), 'plans_sha256': frozen['plans_sha256'],
        'predictions_sha256': N.sha(out / 'predictions.jsonl'), 'counts': dict(counts),
        'python': sys.version, 'platform': platform.platform(), 'elapsed_s': time.time() - start,
        'implementation': code_hashes(), 'shard_index': shard_index, 'shard_count': shard_count}
    N.write_json(out / 'run.json', receipt)
    return receipt


def merge(root, runs, out):
    frozen = verify_plan(root)
    shards = {}
    for directory in runs:
        receipt = json.loads((directory / 'run.json').read_text())
        if (receipt['protocol_sha256'] != N.sha(root / 'protocol.json') or receipt['plans_sha256'] != frozen['plans_sha256']
                or receipt['implementation'] != code_hashes() or receipt['predictions_sha256'] != N.sha(directory / 'predictions.jsonl')):
            raise ValueError('shard provenance mismatch')
        if receipt['shard_count'] != len(runs) or receipt['shard_index'] in shards:
            raise ValueError('duplicate or incomplete shards')
        shards[receipt['shard_index']] = directory
    if set(shards) != set(range(len(runs))):
        raise ValueError('incomplete shards')
    out.mkdir(parents=True, exist_ok=False)
    handles = [(shards[i] / 'predictions.jsonl').open() for i in range(len(runs))]
    counts, n = Counter(), 0
    try:
        with (root / 'plans.jsonl').open() as plans, (out / 'predictions.jsonl').open('w') as dest:
            for group in itertools.zip_longest(*handles):
                for line in group:
                    if line is None: continue
                    plan_line = next(plans, None)
                    if plan_line is None: raise ValueError('too many predictions')
                    plan, row = json.loads(plan_line), json.loads(line)
                    if (row['case_id'], row['arm'], row['plan_sha256']) != (plan['case_id'], plan['arm'], N.digest(plan)):
                        raise ValueError('shard order or prediction plan mismatch')
                    counts[str(row['predictions']['permissive'])] += 1
                    dest.write(line)
                    n += 1
            if next(plans, None) is not None: raise ValueError('incomplete predictions')
    finally:
        for handle in handles: handle.close()
    receipt = {'protocol_sha256': N.sha(root / 'protocol.json'), 'plans_sha256': frozen['plans_sha256'],
        'predictions_sha256': N.sha(out / 'predictions.jsonl'), 'counts': dict(counts), 'plans': n,
        'implementation': code_hashes(), 'shards': {str(p): N.sha(p / 'run.json') for p in runs}}
    N.write_json(out / 'run.json', receipt)
    return receipt


def endpoint(row):
    value = row.get('mean_od600')
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError('malformed observed OD')
    return value > .02


def rescue(co, a, b):
    return None if any(x is None for x in (co, a, b)) else co and not a and not b


def prediction_rows(path):
    with path.open() as stream:
        for line in stream:
            if line.strip(): yield json.loads(line)


def evaluate(root, run_root, observations, out):
    frozen = verify_plan(root)
    receipt = json.loads((run_root / 'run.json').read_text())
    if receipt['protocol_sha256'] != N.sha(root / 'protocol.json') or receipt['plans_sha256'] != frozen['plans_sha256'] or receipt['predictions_sha256'] != N.sha(run_root / 'predictions.jsonl'):
        raise ValueError('prediction provenance mismatch')
    plans = {(r['case_id'], r['arm']): r for r in N.read_rows(root / 'plans.jsonl')}
    rows = prediction_rows(run_root / 'predictions.jsonl')
    obs = N.read_rows(observations)
    by_id = {r['case_id']: r for r in obs}
    if len(by_id) != len(obs):
        raise ValueError('duplicate observations')
    mono = {(r['species'], r['environment'], r['strains'][0]): r for r in obs if r['kind'] == 'monoculture'}
    strata, scored, seen = defaultdict(list), [], set()
    for row in rows:
        identity = (row['case_id'], row['arm'])
        if identity not in plans or identity in seen: raise ValueError('extra or duplicated predictions')
        seen.add(identity)
        plan = plans[identity]
        if row['plan_sha256'] != N.digest(plan):
            raise ValueError('prediction plan mismatch')
        p = row['predictions']
        if set(p) != set(PROTOCOL['solve_order']) or any(x is not None and not isinstance(x, bool) for x in p.values()):
            raise ValueError('malformed predictions')
        observed = endpoint(by_id.get(plan['case_id'], {}))
        a, b = [endpoint(mono.get((plan['species'], plan['environment'], m['id'].split(':')[-1]), {})) for m in plan['members']]
        rescue_observed = rescue(observed, a, b)
        group = (plan['species'], plan['arm'], 'primary' if plan['primary_environment'] else 'no_nitrogen_control')
        for policy in PROTOCOL['policies']:
            strata[group + (policy, 'bulk_growth')].append((observed, p[policy]))
        strata[group + ('permissive', 'paired_rescue_proxy')].append((rescue_observed, rescue(p['permissive'], p['monoculture_0'], p['monoculture_1'])))
        for label, constant in [('always_growth', True), ('never_growth', False)]:
            strata[group + (label, 'bulk_growth')].append((observed, constant))
            strata[group + (label, 'paired_rescue_proxy')].append((rescue_observed, constant))
        scored.append({'case_id': plan['case_id'], 'arm': plan['arm'], 'observed_bulk_growth': observed,
            'observed_monocultures': [a, b], 'observed_rescue_proxy': rescue_observed, 'predictions': p})
    if seen != set(plans): raise ValueError('incomplete predictions')
    report = {'protocol_sha256': N.sha(root / 'protocol.json'), 'run_sha256': N.sha(run_root / 'run.json'),
        'observations_sha256': N.sha(observations), 'interpretation': PROTOCOL['interpretation'],
        'strata': [{'species': k[0], 'arm': k[1], 'scope': k[2], 'method': k[3], 'endpoint': k[4], **N.confusion(v)} for k, v in sorted(strata.items())],
        'cases': scored}
    if out.exists():
        raise ValueError('evaluation output exists')
    N.write_json(out, report)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    for command, keys in [('prepare', ['inputs', 'registration', 'out']), ('run', ['plan', 'out', 'ec-model', 'bs-model']),
                          ('evaluate', ['plan', 'run', 'observations', 'out']), ('merge', ['plan', 'out'])]:
        child = sub.add_parser(command)
        for key in keys: child.add_argument('--' + key, type=Path, required=True)
        if command == 'run':
            child.add_argument('--shard-index', type=int, default=0)
            child.add_argument('--shard-count', type=int, default=1)
        if command == 'merge': child.add_argument('--runs', type=Path, nargs='+', required=True)
    args = parser.parse_args()
    if args.command == 'prepare':
        r = prepare(args.inputs, args.registration, args.out)
        print(json.dumps({k: r[k] for k in ('cases', 'primary_cases', 'plans', 'plans_sha256')}))
    elif args.command == 'run':
        print(json.dumps(run(args.plan, args.out, {'Ec': args.ec_model, 'Bs': args.bs_model}, args.shard_index, args.shard_count)))
    elif args.command == 'merge':
        print(json.dumps(merge(args.plan, args.runs, args.out)))
    else:
        r = evaluate(args.plan, args.run, args.observations, args.out)
        print(json.dumps({'cases': len(r['cases']), 'strata': r['strata']}))
