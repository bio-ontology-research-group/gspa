#!/usr/bin/env python3
"""Independent SBML/GLPK check of absence scenarios; no production solver imports."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import cobra


def simulate(original, manifest, hypothesis, condition, knockout=False):
    # Each solve owns a model: no bound or GPR changes can leak between cells.
    model = original.copy()
    model.solver = 'glpk'
    request, base = manifest['absence_comparison'], manifest['baseline']
    bounds = {r.id: list(r.bounds) for r in model.reactions}

    def uptake(rid, rate, exact=False):
        reaction = model.reactions.get_by_id(rid)
        coeff = next(iter(reaction.metabolites.values()))
        if coeff < 0:
            bounds[rid][0] = -rate if exact else min(bounds[rid][0], -rate)
        else:
            bounds[rid][1] = rate if exact else max(bounds[rid][1], rate)

    def close(rid):
        bounds[rid][0] = max(0, bounds[rid][0])
        bounds[rid][1] = min(0, bounds[rid][1])

    for reaction in model.reactions:
        if len(reaction.metabolites) == 1:
            if reaction.id in base.get('medium', {}):
                uptake(reaction.id, base['medium'][reaction.id], exact=True)
            else:
                coeff = next(iter(reaction.metabolites.values()))
                if coeff < 0:
                    bounds[reaction.id][0] = max(0, bounds[reaction.id][0])
                else:
                    bounds[reaction.id][1] = min(0, bounds[reaction.id][1])
    maintenance = manifest['objectives']['maintenance']
    if maintenance != 'none':
        bounds[maintenance['reaction']][0] = max(bounds[maintenance['reaction']][0], maintenance['lower_bound'])
    supply = request['external_supply']
    if supply:
        exchange = model.reactions.get_by_id(supply['exchange'])
        species = next(iter(exchange.metabolites))
        compartments = base.get('extracellular_compartments', ['e'])
        for reaction in species.reactions:
            coeff = reaction.metabolites[species]
            if any(m.compartment not in compartments and v * coeff < 0 for m, v in reaction.metabolites.items()):
                if reaction.id not in supply['transport_reactions']:
                    close(reaction.id)
        if hypothesis == 'external_supply':
            uptake(supply['exchange'], supply['uptake_limit'], exact=True)
    for mod in condition['mods']:
        if mod['op'] == 'open':
            uptake(mod['exchange'], mod['rate'])
        elif mod['op'] == 'close':
            reaction = model.reactions.get_by_id(mod['exchange'])
            if next(iter(reaction.metabolites.values())) < 0:
                bounds[reaction.id][0] = max(0, bounds[reaction.id][0])
            else:
                bounds[reaction.id][1] = min(0, bounds[reaction.id][1])
        elif mod['op'] == 'knockout':
            for rid in mod['reactions']:
                close(rid)
        else:
            raise ValueError('unsupported oracle intervention')
    for constraint in manifest.get('condition_constraints', []):
        rid = constraint['substrate_exchange']
        coefficient = next(iter(model.reactions.get_by_id(rid).metabolites.values()))
        if bounds[rid][0] < 0 if coefficient < 0 else bounds[rid][1] > 0:
            for closed in constraint['closed_reactions']:
                close(closed)
    if hypothesis in ('pathway_loss', 'external_supply') or (hypothesis == 'missing_annotation' and knockout):
        close(request['reaction'])
    if any(lo > hi for lo, hi in bounds.values()):
        return {'status': 'infeasible', 'growth': None}
    for rid, values in bounds.items():
        model.reactions.get_by_id(rid).bounds = tuple(values)
    objective = manifest['objectives']['biomass']
    if objective != 'auto':
        model.objective = model.reactions.get_by_id(objective)
        model.objective_direction = 'max'
    growth = model.slim_optimize()
    return {'status': model.solver.status, 'growth': max(0, float(growth)) if math.isfinite(growth) else None}


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()


def compare(manifest_path, result_path):
    manifest, result = (json.loads(p.read_text()) for p in (manifest_path, result_path))
    model_path = manifest_path.parent / manifest['model']
    if sha(model_path) != result['input_model_sha256']:
        raise ValueError('model changed since production execution')
    problem = {k: manifest[k] for k in ('baseline', 'objectives', 'absence_comparison', 'condition_constraints') if k in manifest}
    if problem != result['problem']:
        raise ValueError('assumptions changed since production execution')
    model = cobra.io.read_sbml_model(model_path, f_replace={})
    conditions = {c['id']: c for c in manifest['absence_comparison']['conditions']}
    comparisons = []
    for cell in result['cells']:
        for perturbation, record in cell['phenotypes'].items():
            if record['status'] == 'unavailable':
                continue
            oracle = simulate(model, manifest, cell['hypothesis'], conditions[cell['condition_id']], perturbation != 'wild_type')
            agrees = oracle['status'] == record['status'] and (
                math.isclose(oracle['growth'], record['growth'], abs_tol=1e-7, rel_tol=1e-6)
                if oracle['status'] == 'optimal' else oracle['growth'] is None and record['growth'] is None)
            comparisons.append({'condition_id': cell['condition_id'], 'hypothesis': cell['hypothesis'],
                                'perturbation': perturbation, 'oracle': oracle, 'production': record, 'agrees': agrees})
    return {'schema': 'gspa.absence_cobra_check/1', 'interpretation': 'Independent implementation agreement, not biological accuracy.',
            'cobra_version': cobra.__version__, 'solver': 'glpk',
            'hashes': {'manifest': sha(manifest_path), 'result': sha(result_path), 'model': sha(model_path)},
            'comparisons': comparisons, 'agrees': bool(comparisons) and all(c['agrees'] for c in comparisons)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--result', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    report = compare(args.manifest, args.result)
    args.out.write_text(json.dumps(report, indent=2) + '\n')
    return 0 if report['agrees'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
