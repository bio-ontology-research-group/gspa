#!/usr/bin/env python3
"""Known-GPR perturbation and intracellular capacity on unmodified Sco-GEM.

Not an experiment-medium reconstruction, extracellular production prediction or
function recovery benchmark. Run on a compute node. No observation file is read.
"""
import argparse
import copy
import csv
import hashlib
import importlib.util
import importlib.metadata
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
import metabolic_experiments as ME
import numpy as np
import cobra
import libsbml

GENE = 'G_SCO7221'
BIO = 'R_BIOMASS_SCO_tRNA'
PRODUCTS = ['R_DM_germicidin' + c + '_c' for c in 'ABCD']
PROTOCOL = {'schema': 'gspa.germicidin_capacity_diagnostic/1',
    'conditions': ['published_model_minimal', 'phosphate_closed', 'carbon_closed'],
    'genotypes': ['native', 'gcs_deleted', 'native_repeat'],
    'requirements': ['growth_fraction', 'maintenance_only'], 'growth_fraction': .95,
    'growth_threshold': 1e-6, 'flux_tolerance': 1e-7, 'max_solver_calls': 500,
    'product_targets': PRODUCTS, 'biomass': BIO, 'maintenance': {'reaction': 'R_ATPM', 'lower_bound': 2.64},
    'native_gene': GENE,
    'limitations': ['Known source GPRs; no masked discovery or gene assignment accuracy.',
        'The bundled model minimal medium is not the experimental complex GYM medium.',
        'Phosphate/carbon closure are declared model interventions, not the published batch time course.',
        'Intracellular demand flux is synthesis capacity, not extracellular secretion or expression.',
        'Source-model A/B/C/D congeners are not assumed to match the publication A/B7 labels.',
        'native_repeat is a state-leakage/determinism control after deletion, not complementation or plasmid expression.',
        'Growth-floor capacities use 95% of each genotype-specific optimum, not one shared absolute growth rate.',
        'Flaviolin and the publication total-production quantity have no matched diagnostic endpoint.',
        'The source model is M145, not a complete reconstruction of M1152 or M1317.']}


def manifest(model, medium_path, condition):
    with medium_path.open() as handle:
        rows = list(csv.DictReader(handle, skipinitialspace=True))
    medium = {'R_' + row['exchange'].strip(): -float(row['uptake']) for row in rows}
    if condition == 'phosphate_closed':
        medium.pop('R_EX_pi_e')
    elif condition == 'carbon_closed':
        medium.pop('R_EX_glc__D_e')
    external = [r for r in model.rxn_ids if model.is_single_metabolite(r)
                and model.met_compartment[model.exchange_metabolite(r)] == 'e']
    return {'model': model.path, 'objectives': {'biomass': BIO, 'maintenance': PROTOCOL['maintenance']},
            'baseline': {'medium': medium, 'exchanges': external, 'carbon_exchanges': ['R_EX_glc__D_e']},
            'context': {}, 'config': {'max_solver_calls': PROTOCOL['max_solver_calls'],
                'max_focus_reactions': 0, 'growth_fraction': PROTOCOL['growth_fraction']}}


def oracle(original, m, genotype, requirement):
    """Independent parser/GPR deletion, bounds construction and GLPK objectives."""
    other = original.copy()
    medium = m['baseline']['medium']
    for r in other.reactions:
        if len(r.metabolites) == 1:
            coefficient = next(iter(r.metabolites.values()))
            lo, hi = r.bounds
            r.bounds = ((-medium[r.id] if r.id in medium else max(0, lo)), hi) if coefficient < 0 else (
                lo, medium[r.id] if r.id in medium else min(0, hi))
    atpm = other.reactions.get_by_id('R_ATPM')
    atpm.lower_bound = max(atpm.lower_bound, PROTOCOL['maintenance']['lower_bound'])
    if genotype == 'gcs_deleted':
        other.genes.get_by_id(GENE).knock_out()
    bio = other.reactions.get_by_id(BIO)
    other.objective = bio; other.objective_direction = 'max'
    growth = other.slim_optimize()
    out = {'growth_status': other.solver.status,
           'growth': float(growth) if other.solver.status == 'optimal' else None, 'capacities': {}}
    if out['growth_status'] != 'optimal' or (requirement == 'growth_fraction' and growth <= PROTOCOL['growth_threshold']):
        return out
    if requirement == 'growth_fraction':
        bio.lower_bound = max(bio.lower_bound, PROTOCOL['growth_fraction'] * growth)
    for rid in PRODUCTS:
        rxn = other.reactions.get_by_id(rid)
        sign = -float(np.sign(next(iter(rxn.metabolites.values()))))
        values, statuses = [], []
        for sense in ('min', 'max'):
            other.objective = other.problem.Objective(sign * rxn.flux_expression, direction=sense)
            value = other.slim_optimize(); statuses.append(other.solver.status)
            values.append(float(value) if other.solver.status == 'optimal' else None)
        out['capacities'][rid] = {'minimum': values[0], 'maximum': values[1], 'statuses': statuses}
    return out


def production(model, m, genotype, requirement, module=ME, gene=GENE):
    run = module.Runner(model, m, module.resolve_config(m))
    lb, ub = run.base_lb.copy(), run.base_ub.copy()
    mods = [{'op': 'knockout', 'genes': [gene], 'reactions': []}] if genotype == 'gcs_deleted' else []
    killed = module.apply_mods(model, m, run.cfg, lb, ub, mods, run.external)
    g, growth = run._growth(lb, ub, 'growth')
    out = {'growth_status': g.status, 'growth': growth if g.status == 'optimal' else None, 'capacities': {}}
    if g.status == 'optimal' and (requirement == 'maintenance_only' or growth > PROTOCOL['growth_threshold']):
        if requirement == 'growth_fraction':
            lb = run._fixed_growth_lb(lb, growth)
        for rid in PRODUCTS:
            col, _ = model.exchange(rid)
            assert model.met_compartment[model.exchange_metabolite(rid)] == 'c'
            c = np.zeros(len(model.rxn_ids)); c[col] = model.secretion_sign(rid)
            lo = run.solver.optimize(c, lb, ub, rid + ':min')
            hi = run.solver.optimize(-c, lb, ub, rid + ':max')
            out['capacities'][rid] = {'minimum': lo.objective if lo.status == 'optimal' else None,
                'maximum': -hi.objective if hi.status == 'optimal' else None, 'statuses': [lo.status, hi.status]}
    return out, killed, run.solver.log


def agree(a, b):
    if isinstance(a, dict):
        return isinstance(b, dict) and set(a) == set(b) and all(agree(a[k], b[k]) for k in a)
    if isinstance(a, (float, int)) and not isinstance(a, bool):
        return isinstance(b, (float, int)) and math.isclose(a, b, rel_tol=1e-6, abs_tol=1e-7)
    return a == b


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('model', 'medium', 'out'):
        p.add_argument('--' + name, type=Path, required=True)
    p.add_argument('--before', type=Path, required=True)
    args = p.parse_args(); args.out.mkdir(parents=True, exist_ok=False)
    def save(name, data):
        (args.out / name).write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')
    # Freeze the declared plan and actual input bytes before any LP is solved.
    save('protocol.json', {**PROTOCOL, 'model_sha256': hashlib.sha256(args.model.read_bytes()).hexdigest(),
                          'medium_sha256': hashlib.sha256(args.medium.read_bytes()).hexdigest(),
                          'before_source_sha256': hashlib.sha256(args.before.read_bytes()).hexdigest(),
                          'after_source_sha256': hashlib.sha256(Path(ME.__file__).read_bytes()).hexdigest(),
                          'diagnostic_source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    save('versions.json', {p: importlib.metadata.version(p) for p in ('numpy', 'scipy', 'python-libsbml', 'cobra', 'swiglpk')})
    spec = importlib.util.spec_from_file_location('before_knockout_identity', args.before)
    before = importlib.util.module_from_spec(spec); sys.modules[spec.name] = before; spec.loader.exec_module(before)
    model = ME.load_sbml_model(args.model)
    old_model = before.load_sbml_model(args.model)
    other = cobra.io.read_sbml_model(args.model, f_replace={}); other.solver = 'glpk'
    rows = []
    for condition in PROTOCOL['conditions']:
        m = manifest(model, args.medium, condition)
        save('manifest-' + condition + '.json', m)
        for genotype in PROTOCOL['genotypes']:
            for requirement in PROTOCOL['requirements']:
                actual, killed, calls = production(model, m, genotype, requirement)
                previous, old_killed, old_calls = production(old_model, m, genotype, requirement, before)
                expected = oracle(other, m, genotype, requirement)
                rows.append({'condition': condition, 'genotype': genotype, 'requirement': requirement,
                             'production': actual, 'oracle': expected, 'knocked_out_reactions': killed,
                             'solver_calls': calls, 'agrees': agree(actual, expected),
                             'before': previous, 'before_knocked_out_reactions': old_killed,
                             'before_solver_calls': old_calls, 'unchanged': actual == previous and killed == old_killed})
                save('results.json', rows)
    # Known-GPR expectations are mechanism tests, not independent function truth.
    for row in rows:
        if row['genotype'] == 'gcs_deleted':
            assert row['knocked_out_reactions'] == ['R_GERMA', 'R_GERMB', 'R_GERMC', 'R_GERMD']
            for capacity in row['production']['capacities'].values():
                if capacity['statuses'] == ['optimal', 'optimal']:
                    assert abs(capacity['maximum']) < PROTOCOL['flux_tolerance']
        elif row['genotype'] == 'native_repeat':
            native = next(r for r in rows if r['condition'] == row['condition'] and r['requirement'] == row['requirement'] and r['genotype'] == 'native')
            assert row['production'] == native['production']
    # Biological-mechanism regressions: loss of this secondary activity must not
    # be confused with loss of growth, and no-growth must not imply no synthesis.
    def result(condition, genotype, requirement):
        return next(r['production'] for r in rows if (r['condition'], r['genotype'], r['requirement']) ==
                    (condition, genotype, requirement))
    native = result('published_model_minimal', 'native', 'maintenance_only')
    deleted = result('published_model_minimal', 'gcs_deleted', 'maintenance_only')
    assert native['growth'] > PROTOCOL['growth_threshold']
    assert math.isclose(native['growth'], deleted['growth'], rel_tol=1e-6, abs_tol=1e-7)
    for capacity in native['capacities'].values():
        assert capacity['maximum'] > PROTOCOL['flux_tolerance'] and abs(capacity['minimum']) < PROTOCOL['flux_tolerance']
    stationary = result('phosphate_closed', 'native', 'maintenance_only')
    assert stationary['growth'] < PROTOCOL['growth_threshold']
    assert all(c['maximum'] > PROTOCOL['flux_tolerance'] for c in stationary['capacities'].values())
    assert result('phosphate_closed', 'native', 'growth_fraction')['capacities'] == {}
    assert all(r['production']['growth_status'] == 'infeasible' and not r['production']['capacities']
               for r in rows if r['condition'] == 'carbon_closed')

    # Representation-invariance fixture: only the SBML label changes. This is
    # not a new biological strain or an edited metabolic network.
    doc = libsbml.readSBMLFromFile(str(args.model))
    doc.getModel().getPlugin('fbc').getGeneProduct(GENE).setLabel('SCO7221')
    labelled_path = args.out / 'Sco-GEM-labelled.xml'
    libsbml.writeSBMLToFile(doc, str(labelled_path))
    labelled = ME.load_sbml_model(labelled_path)
    old_labelled = before.load_sbml_model(labelled_path)
    assert labelled.rxn_ids == model.rxn_ids and labelled.met_ids == model.met_ids
    assert (labelled.S != model.S).nnz == 0 and labelled.gpr == model.gpr
    assert np.array_equal(labelled.lb, model.lb) and np.array_equal(labelled.ub, model.ub)
    aliases = []
    for requirement in PROTOCOL['requirements']:
        m = manifest(labelled, args.medium, 'published_model_minimal')
        actual, killed, calls = production(labelled, m, 'gcs_deleted', requirement, gene='SCO7221')
        previous, old_killed, old_calls = production(old_labelled, m, 'gcs_deleted', requirement, before, 'SCO7221')
        expected = result('published_model_minimal', 'gcs_deleted', requirement)
        aliases.append({'requirement': requirement, 'manifest': m, 'before': previous, 'after': actual,
            'before_knocked_out_reactions': old_killed, 'after_knocked_out_reactions': killed,
            'before_solver_calls': old_calls, 'after_solver_calls': calls, 'matches_native_id': actual == expected,
            'derived_label_only_model_sha256': labelled.sha256})
        save('label-regressions.json', aliases)
        assert old_killed == [] and killed == ['R_GERMA', 'R_GERMB', 'R_GERMC', 'R_GERMD']
        assert actual == expected
        assert all(c['maximum'] > PROTOCOL['flux_tolerance'] for c in previous['capacities'].values())
    save('summary.json', {'scenarios': len(rows), 'oracle_agreements': sum(r['agrees'] for r in rows),
        'production_lp_calls': sum(len(r['solver_calls']) for r in rows),
        'native_before_after_pairs_unchanged': sum(r['unchanged'] for r in rows),
        'label_regressions': len(aliases), 'protocol': PROTOCOL,
        'scope': 'Known-GPR, intracellular-capacity mechanism diagnostic. No empirical production or function-recovery accuracy.'})
    assert all(r['agrees'] for r in rows), 'independent solver disagreement retained'
    assert all(r['unchanged'] for r in rows), 'valid native-ID behavior changed'
    print(json.dumps({'scenarios': len(rows), 'oracle_agreements': sum(r['agrees'] for r in rows)}))


if __name__ == '__main__':
    main()
