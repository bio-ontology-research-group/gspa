"""Conditional, reaction-specific comparisons of annotation/sequence loss and supply.

This is an explicit scenario comparison, not a classifier of biological absence.
The original GPR is never consulted: candidate deletion removes only the proposed
assignment, not all activities of a real gene. Measured truth belongs in a separate
evaluation. All scenarios for a condition are reserved before spending LP calls.
"""
from __future__ import annotations

import copy
import itertools
import math
import time

import metabolic_experiments as ME

SCHEMA = 'gspa.absence_experiments/1'
HYPOTHESES = ('missing_annotation', 'missing_sequence', 'pathway_loss', 'external_supply')
PERTURBATIONS = ('wild_type', 'candidate_assignment_knockout')


def _fields(value, required, optional=()):
    if not isinstance(value, dict) or set(value) - set(required) - set(optional) or set(required) - set(value):
        raise ME.ManifestError('invalid absence fields; required: ' + ', '.join(sorted(required)))


def _text(value):
    return isinstance(value, str) and bool(value.strip())


def _number(value, minimum=0):
    return not isinstance(value, bool) and isinstance(value, (int, float)) and math.isfinite(value) and value >= minimum


def transport_reactions(model, exchange, compartments):
    """Model reactions connecting this external species to another compartment."""
    met = model.exchange_metabolite(exchange)
    row = model.met_ids.index(met)
    result = set()
    for col in model.S.getrow(row).indices:
        column = model.S[:, col]
        external_coeff = float(model.S[row, col])
        if any(model.met_compartment[model.met_ids[i]] not in compartments and coeff * external_coeff < 0
               for i, coeff in zip(column.indices, column.data)):
            result.add(model.rxn_ids[col])
    return result


def validate(manifest, model, cfg):
    request = manifest['absence_comparison']
    _fields(request, ('reaction', 'candidate_function_key', 'assembly', 'function_search', 'external_supply', 'conditions'))
    rid = request['reaction']
    if not isinstance(rid, str) or rid not in model.rxn_index or model.is_single_metabolite(rid):
        raise ME.ManifestError('absence comparison requires an internal model reaction')
    candidate = request['candidate_function_key']
    if candidate is not None and not _text(candidate):
        raise ME.ManifestError('candidate_function_key must be nonempty or null')
    assembly = request['assembly']
    _fields(assembly, ('status', 'completeness_fraction', 'source'))
    if assembly['status'] not in ('COMPLETE', 'DRAFT', 'UNKNOWN') or not _text(assembly['source']):
        raise ME.ManifestError('assembly requires explicit status and source')
    fraction = assembly['completeness_fraction']
    if fraction is not None and (not _number(fraction) or fraction > 1):
        raise ME.ManifestError('assembly completeness must be null or in [0,1]')
    search = request['function_search']
    _fields(search, ('status', 'source'))
    if search['status'] not in ('COMPLETE', 'INCOMPLETE', 'UNKNOWN') or not _text(search['source']):
        raise ME.ManifestError('function_search requires explicit status and source')
    # An assembly statistic or family search does not establish functional absence.
    external = ME.validate_boundary(model, manifest)
    run = ME.Runner(model, manifest, cfg)
    if rid == model.rxn_ids[run.bio] or (manifest['objectives']['maintenance'] != 'none'
            and rid == manifest['objectives']['maintenance']['reaction']):
        raise ME.ManifestError('the focal reaction cannot be the biomass or maintenance objective')
    supply = request['external_supply']
    closed_ports = []
    if supply is not None:
        _fields(supply, ('exchange', 'uptake_limit', 'transport_reactions', 'source'))
        if supply['exchange'] not in external or not _number(supply['uptake_limit']) or supply['uptake_limit'] <= 0:
            raise ME.ManifestError('supply requires a declared extracellular exchange and positive finite uptake limit')
        ports = supply['transport_reactions']
        if not isinstance(ports, list) or not ports or any(not isinstance(p, str) for p in ports) or len(set(ports)) != len(ports):
            raise ME.ManifestError('supply requires unique explicit transport reactions')
        possible = transport_reactions(model, supply['exchange'], set(manifest['baseline'].get('extracellular_compartments', ['e'])))
        if not set(ports) <= possible or rid in possible or not _text(supply['source']):
            raise ME.ManifestError('supply ports must connect its external species to an internal compartment, independently of the focal reaction')
        closed_ports = sorted(possible - set(ports))
    conditions = request['conditions']
    if not isinstance(conditions, list) or not conditions:
        raise ME.ManifestError('absence comparison requires conditions')
    ids = set()
    for condition in conditions:
        _fields(condition, ('id', 'mods'))
        if not _text(condition['id']) or condition['id'] in ids or not isinstance(condition['mods'], list):
            raise ME.ManifestError('conditions need unique IDs and a list of modifications')
        ids.add(condition['id'])
        for mod in condition['mods']:
            op = mod.get('op') if isinstance(mod, dict) else None
            if op == 'open':
                _fields(mod, ('op', 'exchange', 'rate'))
                if not _number(mod['rate']):
                    raise ME.ManifestError('condition uptake rate must be finite and nonnegative')
            elif op == 'close':
                _fields(mod, ('op', 'exchange'))
            elif op == 'knockout':
                _fields(mod, ('op', 'reactions'))
                if not isinstance(mod['reactions'], list) or not mod['reactions'] or any(r not in model.rxn_index for r in mod['reactions']):
                    raise ME.ManifestError('condition knockouts require known reactions')
            else:
                raise ME.ManifestError('absence conditions accept only explicit exchange open/close and reaction knockout')
        ME.apply_mods(model, manifest, cfg, run.base_lb.copy(), run.base_ub.copy(), condition['mods'], external)
    if cfg['strategy'] != 'fixed' or cfg.get('panel') is not None:
        raise ME.ManifestError('absence comparisons use strategy fixed with the explicit conditions in manifest order')
    for key in ('max_solver_calls', 'max_experiments'):
        if isinstance(cfg[key], bool) or not isinstance(cfg[key], int) or cfg[key] < 0:
            raise ME.ManifestError('absence budgets must be nonnegative integers')
    if not _number(cfg['growth_threshold']) or cfg['growth_threshold'] <= 0:
        raise ME.ManifestError('growth_threshold must be positive and finite')
    return run, closed_ports


def phenotype(status, growth, threshold):
    if status == 'optimal':
        return 'growth' if growth > threshold else 'no_growth'
    if status == 'infeasible':
        return 'infeasible_under_constraints'
    return 'unknown'


def close_reaction(model, lb, ub, rid):
    col = model.rxn_index[rid]
    # A stated positive requirement remains a contradiction, not relaxed rescue.
    lb[col], ub[col] = max(lb[col], 0.0), min(ub[col], 0.0)


def simulate(run, request, hypothesis, condition, cfg, closed_ports):
    model, solver = run.model, run.solver
    candidate = request['candidate_function_key']
    unavailable = ('no_retained_candidate' if hypothesis == 'missing_annotation' and candidate is None else
                   'no_explicit_supply_and_transport' if hypothesis == 'external_supply' and request['external_supply'] is None else None)
    cell = {'condition_id': condition['id'], 'hypothesis': hypothesis, 'phenotypes': {}, 'solver_calls': []}
    if unavailable:
        cell.update(status='unavailable', reason=unavailable)
        return cell
    lb, ub = run.base_lb.copy(), run.base_ub.copy()
    for rid in closed_ports:
        close_reaction(model, lb, ub, rid)
    if hypothesis == 'external_supply':
        supply = request['external_supply']
        col, coeff = model.exchange(supply['exchange'])
        # Exact total cap, not an increment or an unbounded rescue.
        if coeff > 0:
            ub[col] = supply['uptake_limit']
        else:
            lb[col] = -supply['uptake_limit']
    # Interventions apply after hypothesized supply, so withdrawal really removes it.
    ME.apply_mods(model, run.manifest, cfg, lb, ub,
                  [m for m in condition['mods'] if m['op'] != 'knockout'], run.external)
    killed = sorted({r for m in condition['mods'] if m['op'] == 'knockout' for r in m['reactions']})
    for rid in killed:
        close_reaction(model, lb, ub, rid)
    active = ME.apply_condition_constraints(model, run.condition_constraints, lb, ub)
    if hypothesis in ('pathway_loss', 'external_supply'):
        close_reaction(model, lb, ub, request['reaction'])
    cell['active_constraints'] = active
    first = solver.used
    for perturbation in PERTURBATIONS:
        if perturbation == 'candidate_assignment_knockout' and candidate is None:
            cell['phenotypes'][perturbation] = {'status': 'unavailable', 'growth': None, 'class': 'unknown'}
            continue
        plb, pub = lb.copy(), ub.copy()
        disabled = list(killed)
        if hypothesis in ('pathway_loss', 'external_supply') or (
                hypothesis == 'missing_annotation' and perturbation == 'candidate_assignment_knockout'):
            close_reaction(model, plb, pub, request['reaction'])
            disabled.append(request['reaction'])
        result, growth = run._growth(plb, pub, condition['id'] + ':' + hypothesis + ':' + perturbation)
        cell['phenotypes'][perturbation] = {
            'status': result.status, 'growth': growth if result.status == 'optimal' else None,
            'class': phenotype(result.status, growth, cfg['growth_threshold']),
            'disabled_reactions': sorted(set(disabled)),
        }
    cell['solver_calls'] = list(range(first, solver.used))
    cell['status'] = 'completed' if all(p['class'] != 'unknown' or p['status'] == 'unavailable'
                                     for p in cell['phenotypes'].values()) else 'inconclusive'
    return cell


def compare(cells):
    comparisons = []
    for left, right in itertools.combinations(HYPOTHESES, 2):
        distinguishing, comparable = [], 0
        by_key = {(c['condition_id'], c['hypothesis']): c for c in cells}
        for cid in dict.fromkeys(c['condition_id'] for c in cells):
            a, b = by_key[cid, left], by_key[cid, right]
            for perturbation in PERTURBATIONS:
                x = a['phenotypes'].get(perturbation, {}).get('class', 'unknown')
                y = b['phenotypes'].get(perturbation, {}).get('class', 'unknown')
                # Both predict failure to grow; LP feasibility alone is not an
                # independently observable phenotype that distinguishes causes.
                x = 'no_growth' if x == 'infeasible_under_constraints' else x
                y = 'no_growth' if y == 'infeasible_under_constraints' else y
                if 'unknown' in (x, y):
                    continue
                comparable += 1
                if x != y:
                    distinguishing.append({'condition_id': cid, 'perturbation': perturbation})
        comparisons.append({'hypotheses': [left, right], 'comparable_phenotypes': comparable,
                            'status': 'distinguishable_under_assumptions' if distinguishing else
                                      'same_evaluated_growth_classes' if comparable else 'unavailable',
                            'distinguishing_interventions': distinguishing})
    return comparisons


def run_absence_experiments(manifest, model, overrides=None):
    started = time.monotonic()
    cfg = ME.resolve_config(manifest, overrides)
    run, closed_ports = validate(manifest, model, cfg)
    request = manifest['absence_comparison']
    hypotheses_available = 2 + (request['candidate_function_key'] is not None) + (request['external_supply'] is not None)
    cost = hypotheses_available * (2 if request['candidate_function_key'] is not None else 1)
    cells, visited = [], []
    for condition in request['conditions']:
        if len(visited) >= cfg['max_experiments'] or cost > run.solver.remaining:
            break
        visited.append(condition['id'])
        cells.extend(simulate(run, request, h, condition, cfg, closed_ports) for h in HYPOTHESES)
    problem = {key: copy.deepcopy(manifest[key]) for key in
               ('baseline', 'objectives', 'absence_comparison', 'condition_constraints') if key in manifest}
    return {'schema': SCHEMA, 'input_model_sha256': model.sha256,
            'problem': problem, 'problem_sha256': ME.canonical_hash(problem),
            'strategy': cfg['strategy'], 'config': cfg, 'hypotheses': list(HYPOTHESES),
            'cells': cells, 'comparisons': compare(cells), 'visited_conditions': visited,
            'unvisited_conditions': [c['id'] for c in request['conditions'] if c['id'] not in visited],
            'solver_calls': run.solver.log, 'max_solver_calls': cfg['max_solver_calls'],
            'closed_undeclared_transport_reactions': closed_ports,
            'objective_source': run.objective_source, 'wall_seconds': time.monotonic() - started,
            'gene_attribution': 'unverified', 'selected_explanation': None,
            'promotion': 'none', 'calibrated_probability': None, 'independent_proof_of_function': False,
            'assumptions': [
                'The focal reaction is available at the supplied model bounds for missing annotation and missing sequence, and disabled for pathway loss and external supply.',
                'Candidate-assignment knockout removes only the proposed focal activity, not all activities of a biological gene. Original GPRs are not consulted.',
                'Missing sequence assumes activity in the organism without a retained sequence; assembly and search metadata alone cannot establish or exclude it.',
                'External supply assumes pathway loss plus the declared bounded supply. These explanations overlap; they are not exclusive probability classes.',
                'Undeclared direct transport of the proposed supplied species is closed in every scenario; declared transport remains a model assumption.',
                'Condition interventions follow the supply assumption, including withdrawal and transport knockout.',
                'Growth classes are conditional on topology, biomass, bounds, maintenance and medium. Identical classes do not establish identical biology.',
                'No measured phenotype or function truth selects a hypothesis. Gene attribution requires independent evaluation.']}
