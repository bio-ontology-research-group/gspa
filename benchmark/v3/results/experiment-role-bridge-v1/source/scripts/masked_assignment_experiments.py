"""Competing gene assignments on a model with the target's GPR removed.

Only the masking boundary reads the original assignment. Planning and LP solves
use the masked model and explicit candidate GPR clauses, including a null case.
Results are conditional phenotype signatures, never self-confirming functions.
"""
from __future__ import annotations

import copy
import json
import math
import random
import re
import time
from dataclasses import replace

import numpy as np
import metabolic_experiments as ME

SCHEMA = 'gspa.masked_assignment_experiments/1'
TARGET = 'MASKED_TARGET'


def fingerprint(model):
    return ME.canonical_hash({
        'metabolites': model.met_ids, 'reactions': model.rxn_ids,
        'S': [model.S.shape, model.S.indptr.tolist(), model.S.indices.tolist(), model.S.data.tolist()],
        'lb': model.lb.tolist(), 'ub': model.ub.tolist(),
        'gpr': {r: sorted(sorted(set(c)) for c in model.gpr.get(r, [])) for r in model.rxn_ids},
        'gene_alias': model.gene_alias, 'compartments': model.met_compartment,
        'formula': model.met_formula, 'objective': model.objective_terms, 'sense': model.objective_sense})


def mask_model(model, target):
    """Remove entire AND clauses containing the target, preserving alternative OR clauses.

    An orphan reaction remains in the network: losing an annotation is not proof
    of pathway loss. No removed reaction list or original gene alias is returned.
    """
    if target not in model.gene_alias:
        raise ME.ManifestError('mask target is not a declared model gene')
    gene = model.gene_alias[target]
    if TARGET in model.gene_alias:
        raise ME.ManifestError('reserved masked target identifier already exists')
    masked = replace(model, path='masked', sha256='',
                     gpr={r: [list(c) for c in model.gpr.get(r, []) if gene not in c]
                          for r in model.rxn_ids},
                     gene_alias={k: v for k, v in model.gene_alias.items() if v != gene})
    masked.sha256 = fingerprint(masked)
    masked.path = 'masked:' + masked.sha256
    return masked


def public_problem(manifest, model):
    private = manifest['assignment_recovery']
    if set(private) != {'mask_gene', 'hypotheses', 'conditions'}:
        raise ME.ManifestError('assignment_recovery requires only mask_gene, hypotheses and conditions')
    target = private['mask_gene']
    if target not in model.gene_alias:
        raise ME.ManifestError('unknown target gene')
    aliases = {k for k, v in model.gene_alias.items() if v == model.gene_alias[target]}
    masked = mask_model(model, target)
    # Deliberately discard the original context, annotations, model path and metadata.
    public = {k: copy.deepcopy(manifest[k]) for k in ('baseline', 'objectives', 'condition_constraints') if k in manifest}
    public['model'] = masked.path
    public['hypotheses'] = copy.deepcopy(private['hypotheses'])
    public['conditions'] = copy.deepcopy(private['conditions'])
    def check(value):
        if isinstance(value, dict):
            if set(value) & {'truth', 'labels', 'expected', 'correct', 'is_true_target', 'mask_gene'}:
                raise ME.ManifestError('private answer field in public problem')
            for k, v in value.items():
                check(k)
                check(v)
        elif isinstance(value, list):
            for v in value:
                check(v)
        elif isinstance(value, str):
            if any(re.search(r'(?<![A-Za-z0-9_])' + re.escape(alias) + r'(?![A-Za-z0-9_])', value) for alias in aliases):
                raise ME.ManifestError('original target identity in public evidence')
    check(public)
    validate(public, masked)
    return public, masked


def validate(public, model):
    hypotheses, conditions = public['hypotheses'], public['conditions']
    if len(hypotheses) < 2 or not conditions:
        raise ME.ManifestError('need competing assignments and at least one condition')
    if len({h['id'] for h in hypotheses}) != len(hypotheses) or len({c['id'] for c in conditions}) != len(conditions):
        raise ME.ManifestError('duplicate hypothesis or condition ID')
    nulls = 0
    known_genes = set(model.gene_alias.values()) | {TARGET}
    for h in hypotheses:
        if set(h) - {'id', 'assignments', 'evidence'} or not re.fullmatch(r'H[0-9]+', h['id']):
            raise ME.ManifestError('hypothesis IDs must be opaque H<number>; extra metadata is not admitted')
        assignments = h['assignments']
        nulls += not bool(assignments)
        for rid, clauses in assignments.items():
            if rid not in model.rxn_index or model.is_single_metabolite(rid):
                raise ME.ManifestError('assignments must name internal model reactions')
            if not clauses or any(not c or TARGET not in c or not set(c) <= known_genes or len(c) != len(set(c)) for c in clauses):
                raise ME.ManifestError('each assignment clause must contain MASKED_TARGET and known partner genes')
        for key, value in h.get('evidence', {}).items():
            if key not in {'sequence', 'domain', 'neighborhood', 'recurrence', 'transport', 'bgc'}:
                raise ME.ManifestError('unsupported evidence channel')
            if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value) or not 0 <= value <= 1:
                raise ME.ManifestError('evidence must be a finite score in [0,1]')
    if nulls != 1:
        raise ME.ManifestError('exactly one explicit null assignment is required')
    external = ME.validate_boundary(model, public)
    for c in conditions:
        if set(c) - {'id', 'mods', 'evidence_channels', 'product_exchange'} or not re.fullmatch(r'C[0-9]+', c['id']):
            raise ME.ManifestError('conditions require opaque C<number> IDs and structured modifications')
        if c.get('product_exchange') and c['product_exchange'] not in external:
            raise ME.ManifestError('product must be a declared extracellular exchange')
        if any(channel not in {'sequence', 'domain', 'neighborhood', 'recurrence', 'transport', 'bgc'} for channel in c.get('evidence_channels', [])):
            raise ME.ManifestError('unsupported condition evidence channel')
        for mod in c['mods']:
            if TARGET in mod.get('genes', []):
                raise ME.ManifestError('target knockout is paired automatically, not part of condition selection')
        lb, ub = model.lb.copy(), model.ub.copy()
        ME.apply_mods(model, public, ME.DEFAULT_CONFIG, lb, ub, c['mods'], external)


def assigned_model(masked, hypothesis):
    gpr = copy.deepcopy(masked.gpr)
    for rid, clauses in hypothesis['assignments'].items():
        gpr[rid] = gpr.get(rid, []) + copy.deepcopy(clauses)
    return replace(masked, gpr=gpr, gene_alias={**masked.gene_alias, TARGET: TARGET})


def _interval(solver, model, rid, lb, ub, sign, label):
    c = np.zeros(len(model.rxn_ids))
    c[model.rxn_index[rid]] = sign
    low = solver.optimize(c, lb, ub, label + ':min')
    high = solver.optimize(-c, lb, ub, label + ':max')
    return {'minimum': low.objective if low.status == 'optimal' else None,
            'maximum': -high.objective if high.status == 'optimal' else None,
            'statuses': [low.status, high.status]}


def simulate_cell(masked, public, hypothesis, condition, cfg, solver):
    model = assigned_model(masked, hypothesis)
    run = ME.Runner(model, public, cfg)
    # All scenarios have identical stoichiometry; share one strict LP-call ledger.
    run.solver = solver
    first = solver.used
    base_lb, base_ub = run.base_lb.copy(), run.base_ub.copy()
    ME.apply_mods(model, public, cfg, base_lb, base_ub, condition['mods'], run.external)
    constraints = ME.apply_condition_constraints(model, run.condition_constraints, base_lb, base_ub)
    condition_knockouts = {gene for mod in condition['mods'] if mod['op'] == 'knockout'
                          for gene in mod.get('genes', [])}
    prefix = condition['id'] + ':' + hypothesis['id']
    output = {'hypothesis_id': hypothesis['id'], 'condition_id': condition['id'],
              'active_constraints': constraints, 'phenotypes': {}}
    for perturbation in ('wild_type', 'target_knockout'):
        lb, ub = base_lb.copy(), base_ub.copy()
        killed = []
        if perturbation == 'target_knockout':
            killed = ME.apply_mods(model, public, cfg, lb, ub,
                [{'op': 'knockout', 'genes': sorted(condition_knockouts | {TARGET}), 'reactions': []}], run.external)
        label = prefix + ':' + perturbation
        result, growth = run._growth(lb, ub, label)
        record = {'status': result.status, 'growth': growth if result.status == 'optimal' else None,
                  'knocked_out_reactions': killed, 'flux_intervals': {}}
        # Alternative feasible optima at the stated growth floor. Never infer expression.
        if result.status == 'optimal':
            growing = growth > cfg['growth_threshold']
            floor_lb = lb.copy()
            floor_lb[run.bio] = max(floor_lb[run.bio], cfg['growth_fraction'] * growth)
            if growing and perturbation == 'wild_type':
                for rid in sorted(hypothesis['assignments']):
                    record['flux_intervals'][rid] = _interval(solver, model, rid, floor_lb, ub, 1.0, label + ':' + rid)
            product = condition.get('product_exchange')
            if product and (growing or cfg['product_requirement'] == 'maintenance_only'):
                product_lb = lb if cfg['product_requirement'] == 'maintenance_only' else floor_lb
                record['product_interval'] = _interval(solver, model, product, product_lb, ub,
                                                       model.secretion_sign(product), label + ':product')
        output['phenotypes'][perturbation] = record
    output['solver_calls'] = list(range(first, solver.used))
    wt, ko = output['phenotypes']['wild_type'], output['phenotypes']['target_knockout']
    if wt['status'] == 'optimal' and wt['growth'] > cfg['growth_threshold']:
        if ko['status'] == 'infeasible':
            output['contrast'] = 'knockout_infeasible_with_maintenance'
        elif ko['status'] == 'optimal':
            output['contrast'] = 'knockout_no_growth' if ko['growth'] < cfg['growth_threshold'] else 'both_grow'
        else:
            output['contrast'] = 'inconclusive'
    else:
        output['contrast'] = 'no_growing_control' if wt['status'] in ('optimal', 'infeasible') else 'inconclusive'
    output['status'] = 'inconclusive' if output['contrast'] in ('inconclusive', 'no_growing_control') else 'completed'
    return output


def condition_cost(public, condition):
    # Reserve the full rectangular comparison; hypothesis ordering cannot decide
    # which candidates receive simulations when the budget is nearly exhausted.
    return sum(2 + 2 * len(h['assignments']) + (4 if condition.get('product_exchange') else 0)
               for h in public['hypotheses'])


def run_assignment_experiments(manifest, original, overrides=None, chat=None):
    start = time.monotonic()
    public, masked = public_problem(manifest, original)
    cfg = ME.resolve_config(manifest, overrides)
    if any(isinstance(cfg[k], bool) or not isinstance(cfg[k], int) or cfg[k] < 0
           for k in ('max_solver_calls', 'max_experiments')):
        raise ME.ManifestError('budgets must be nonnegative integers')
    ME.Runner(masked, public, cfg)  # validate medium/objective/growth floor even at zero budget
    solver = ME.BudgetedSolver(masked, cfg['max_solver_calls'], cfg['residual_tol'])
    remaining = sorted(copy.deepcopy(public['conditions']), key=lambda c: c['id'])
    strategy = cfg['strategy']
    if strategy not in ('fixed', 'random', 'evidence', 'qwen'):
        raise ME.ManifestError('unsupported assignment strategy')
    if strategy == 'random':
        random.Random(cfg['seed']).shuffle(remaining)
    elif strategy == 'fixed' and cfg.get('panel') is not None:
        by_id = {c['id']: c for c in remaining}
        panel = cfg['panel']
        if len(set(panel)) != len(panel) or not set(panel) <= set(by_id):
            raise ME.ManifestError('fixed panel must contain unique admissible conditions')
        remaining = [by_id[c] for c in panel]
    elif strategy == 'evidence':
        def support(c):
            return sum(max((h.get('evidence', {}).get(channel, 0) for h in public['hypotheses']), default=0)
                       for channel in set(c.get('evidence_channels', [])))
        remaining.sort(key=lambda c: (-support(c), c['id']))
    cells, selection, visited = [], [], []
    if strategy == 'qwen' and chat is None:
        chat = ME.make_qwen_chat(cfg)
    while remaining and len(visited) < cfg['max_experiments']:
        affordable = [c for c in remaining if condition_cost(public, c) <= solver.remaining]
        if not affordable:
            break
        chosen = affordable[0]
        if strategy == 'qwen':
            prompt = [{'role': 'system', 'content': 'Choose one admissible condition to distinguish candidate assignments. The target identity and original assignments are masked. Simulation is conditional, not independent confirmation. Reply only with JSON {"experiment_id":"C<number>"}.'},
                      {'role': 'user', 'content': json.dumps({'hypotheses': public['hypotheses'], 'conditions': affordable,
                          'previous_results': cells, 'solver_calls_remaining': solver.remaining}, sort_keys=True)}]
            try:
                reply = chat(prompt)
                parsed = json.loads(reply)
                cid = parsed.get('experiment_id')
                if cid not in {c['id'] for c in affordable}:
                    raise ValueError('inadmissible condition')
                chosen = next(c for c in affordable if c['id'] == cid)
                selection.append({'status': 'agent_selected', 'condition_id': cid, 'reply': reply})
            except Exception as error:
                # Preserve failures; do not credit a fallback choice to the agent.
                selection.append({'status': 'agent_failed', 'error': type(error).__name__})
                break
        else:
            selection.append({'status': 'selected', 'condition_id': chosen['id']})
        visited.append(chosen['id'])
        remaining.remove(chosen)
        for h in sorted(public['hypotheses'], key=lambda h: h['id']):
            cells.append(simulate_cell(masked, public, h, chosen, cfg, solver))
    return {'schema': SCHEMA, 'model': {'sha256': masked.sha256, 'path': masked.path},
            'input_model_sha256': original.sha256,
            'public_problem': public, 'public_problem_sha256': ME.canonical_hash(public),
            'strategy': strategy, 'config': cfg, 'visited_conditions': visited,
            'selection': selection, 'cells': cells, 'solver_calls': solver.log,
            'agent_usage': getattr(chat, 'usage', None),
            'unvisited_conditions': [c['id'] for c in public['conditions'] if c['id'] not in visited],
            'max_solver_calls': cfg['max_solver_calls'], 'wall_seconds': time.monotonic() - start,
            'promotion': 'none', 'calibrated_probability': None,
            'assumptions': ['Only the target annotation is masked; network reactions remain available.',
                'Candidate GPR clauses are added as alternatives to surviving non-target clauses.',
                'Exactly one null assignment is included.',
                'Growth and flux ranges are conditional on the model, medium and maintenance.',
                'Feasible product flux is not evidence of observed production or expression.',
                'No hidden function truth or measured outcome is supplied to the planner.',
                'Independent phenotypes and function truth must be joined in a separate evaluation.']}
