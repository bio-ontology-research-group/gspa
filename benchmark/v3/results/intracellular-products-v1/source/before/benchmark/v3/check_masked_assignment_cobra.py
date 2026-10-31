#!/usr/bin/env python3
"""Independent COBRApy/GLPK replay of masked assignment cells and flux intervals.

This checks execution semantics, not whether a candidate is biologically correct.
It parses the original SBML independently and uses COBRA's GPR removal/simplifier.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET

import cobra
from cobra.manipulation.delete import remove_genes


def simulate(model, public, hypothesis, condition, cfg, knockout):
    # SBML may declare unused genes whose _model is unset in COBRApy. Their
    # functional flag is not restored by a model context after knock_out().
    # Isolate every cell/perturbation, including solver and GPR state, so an
    # earlier null/partner-KO experiment cannot contaminate later hypotheses.
    model = model.copy()
    with model:
        for rid, clauses in hypothesis['assignments'].items():
            rxn = model.reactions.get_by_id(rid)
            added = ' or '.join('(' + ' and '.join(c) + ')' for c in clauses)
            rxn.gene_reaction_rule = f'({rxn.gene_reaction_rule}) or ({added})' if rxn.gene_reaction_rule else added
        medium = dict(public['baseline']['medium'])
        extra_kos = []
        for mod in condition['mods']:
            op = mod['op']
            if op in ('carbon_swap', 'nitrogen_swap'):
                for rid in public['baseline'].get('carbon_exchanges' if op == 'carbon_swap' else 'nitrogen_exchanges', []):
                    medium.pop(rid, None)
                medium[mod['exchange']] = cfg['uptake_rate']
            elif op == 'oxygen_removal':
                medium.pop(public['baseline']['oxygen_exchange'], None)
            elif op == 'close':
                medium.pop(mod['exchange'], None)
            elif op == 'open':
                medium[mod['exchange']] = mod.get('rate', cfg['uptake_rate'])
            elif op == 'knockout':
                extra_kos.append(mod)
            else:
                raise ValueError('unsupported modification')
        for rxn in model.reactions:
            if len(rxn.metabolites) != 1:
                continue
            coef = next(iter(rxn.metabolites.values()))
            lb, ub = rxn.bounds
            if coef < 0:
                lb = -medium[rxn.id] if rxn.id in medium else max(lb, 0)
            else:
                ub = medium[rxn.id] if rxn.id in medium else min(ub, 0)
            if lb > ub:
                return {'status': 'infeasible', 'growth': None, 'flux_intervals': {}}
            rxn.bounds = lb, ub
        maintenance = public['objectives']['maintenance']
        if maintenance != 'none':
            rxn = model.reactions.get_by_id(maintenance['reaction'])
            rxn.lower_bound = max(rxn.lower_bound, maintenance['lower_bound'])
        for con in public.get('condition_constraints', []):
            if medium.get(con['substrate_exchange'], 0) > 0:
                for rid in con['closed_reactions']:
                    rxn = model.reactions.get_by_id(rid)
                    if rxn.lower_bound > 0 or rxn.upper_bound < 0:
                        return {'status': 'infeasible', 'growth': None, 'flux_intervals': {}}
                    rxn.bounds = 0, 0
        for mod in extra_kos:
            for rid in mod.get('reactions', []):
                model.reactions.get_by_id(rid).knock_out()
            for gene in mod.get('genes', []):
                model.genes.get_by_id(gene).knock_out()
        if knockout and 'MASKED_TARGET' in model.genes:
            model.genes.get_by_id('MASKED_TARGET').knock_out()
        bio = model.reactions.get_by_id(public['objectives']['biomass'])
        model.objective = bio
        model.objective_direction = 'max'
        growth = model.slim_optimize()
        result = {'status': model.solver.status, 'growth': float(growth) if math.isfinite(growth) else None, 'flux_intervals': {}}
        if result['status'] != 'optimal':
            return result
        original_floor = bio.lower_bound
        growing = growth > cfg['growth_threshold']
        targets = [] if knockout or not growing else [(rid, False) for rid in sorted(hypothesis['assignments'])]
        if condition.get('product_exchange') and (growing or cfg['product_requirement'] == 'maintenance_only'):
            targets.append((condition['product_exchange'], True))
        for rid, product in targets:
            bio.lower_bound = original_floor if product and cfg['product_requirement'] == 'maintenance_only' else max(original_floor, cfg['growth_fraction'] * growth)
            rxn = model.reactions.get_by_id(rid)
            sign = (-1 if next(iter(rxn.metabolites.values())) > 0 else 1) if product else 1
            model.objective = model.problem.Objective(sign * rxn.flux_expression, direction='min')
            minimum = model.slim_optimize()
            statuses = [model.solver.status]
            model.objective_direction = 'max'
            maximum = model.slim_optimize()
            statuses.append(model.solver.status)
            interval = {'minimum': float(minimum) if statuses[0] == 'optimal' else None,
                        'maximum': float(maximum) if statuses[1] == 'optimal' else None, 'statuses': statuses}
            if product:
                result['product_interval'] = interval
            else:
                result['flux_intervals'][rid] = interval
        return result


def agrees(a, b):
    if isinstance(a, dict):
        return isinstance(b, dict) and all(k in b and agrees(v, b[k]) for k, v in a.items() if k != 'knocked_out_reactions')
    if isinstance(a, (int, float)):
        return isinstance(b, (int, float)) and math.isclose(a, b, rel_tol=1e-6, abs_tol=1e-7)
    return a == b


def compare(model_path, manifest, result):
    sha = hashlib.sha256(Path(model_path).read_bytes()).hexdigest()
    if result['input_model_sha256'] != sha:
        raise ValueError('model differs from production input')
    model = cobra.io.read_sbml_model(model_path, f_replace={})
    model.solver = 'glpk'
    target = manifest['assignment_recovery']['mask_gene']
    # COBRA preserves IDs but may omit FBC labels; resolve those independently
    # from XML, without using the production parser's alias table.
    target_ids = {target}
    for node in ET.parse(model_path).getroot().iter():
        if node.tag.rsplit('}', 1)[-1] == 'geneProduct':
            attrs = {k.rsplit('}', 1)[-1]: v for k, v in node.attrib.items()}
            if target in (attrs.get('id'), attrs.get('label')):
                target_ids.add(attrs['id'])
    genes = [g for g in model.genes if g.id in target_ids]
    if len(genes) != 1:
        raise ValueError('mask target is ambiguous or missing')
    remove_genes(model, [genes[0].id], remove_reactions=False)
    public = result['public_problem']
    hypotheses = {h['id']: h for h in public['hypotheses']}
    conditions = {c['id']: c for c in public['conditions']}
    rows = []
    for cell in result['cells']:
        h, c = hypotheses[cell['hypothesis_id']], conditions[cell['condition_id']]
        for kind in ('wild_type', 'target_knockout'):
            oracle = simulate(model, public, h, c, result['config'], kind == 'target_knockout')
            rows.append({'hypothesis_id': h['id'], 'condition_id': c['id'], 'perturbation': kind,
                         'agrees': agrees(cell['phenotypes'][kind], oracle), 'oracle': oracle})
    return {'schema': 'gspa.masked_assignment_oracle/1', 'model_sha256': sha,
            'cobra_version': cobra.__version__, 'solver': 'glpk', 'comparisons': rows,
            'n_comparisons': len(rows), 'n_agree': sum(r['agrees'] for r in rows),
            'interpretation': 'Independent implementation agreement, not function or phenotype accuracy.'}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('model', 'manifest', 'result', 'out'):
        p.add_argument('--' + name, required=True, type=Path)
    args = p.parse_args()
    report = compare(args.model, json.loads(args.manifest.read_text()), json.loads(args.result.read_text()))
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    raise SystemExit(0 if report['n_agree'] == report['n_comparisons'] else 1)
