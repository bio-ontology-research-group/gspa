#!/usr/bin/env python3
"""Independent SBML/GPR/shared-pool GLPK replay; never imports the production LP builder."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import time

import cobra
import swiglpk
from check_community_cobra import Oracle


def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for chunk in iter(lambda: f.read(1048576), b''): h.update(chunk)
    return h.hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


class NutrientOracle(Oracle):
    def __init__(self, path, context, protocol):
        super().__init__(path, {'biomass_reaction': context['biomass'],
            'maintenance_reaction': context['maintenance']['reaction'],
            'ceiling': protocol['pool_removal_ceiling'], 'growth_threshold': protocol['growth_threshold']},
            {'cross_flux': 1e-7})
        self.protocol = protocol
        if set(context['ports']) != self.port_ids:
            raise ValueError('independent extracellular port mapping disagrees')
        for rid, port in context['ports'].items():
            met = next(iter(self.base.reactions.get_by_id(rid).metabolites))
            if port['pool'] != met.id or port['metabolite'] != met.id:
                raise ValueError('independent pool chemical mapping disagrees')

    def blocked(self, genes):
        if not set(genes) <= {g.id for g in self.base.genes}:
            raise ValueError('unknown canonical gene')
        return sorted(r.id for r in self.base.reactions if not r.gpr.eval(set(genes)))

    def simulate(self, plan, kind):
        selected = plan['members']
        if kind.startswith('monoculture_'):
            selected = [selected[int(kind[-1])]]
        # Frozen study strains each have one deletion; do not merge two genotypes.
        if any(len(m['genes']) != 1 for m in selected):
            raise ValueError('oracle requires the registered single-deletion strains')
        members = [(m['id'], m['genes'][0]) for m in selected]
        blocks = {m['id']: set(self.blocked(m['genes'])) for m in selected}
        comp = self.pool_model(members)
        supply = plan['supply_exchanges']
        forbid = self.port_ids - {r for r, rate in supply.items() if rate > 0} if kind == 'no_unsupplied_uptake' else set()
        native = kind == 'native_unsupplied'
        with comp:
            self.configure(comp, members, supply, None if native else self.protocol['port_ceiling'], forbid, blocks)
            if native:
                for rid, rate in supply.items():
                    if rate <= 0: continue
                    base = self.base.reactions.get_by_id(rid)
                    coefficient = next(iter(base.metabolites.values()))
                    for member, _ in members:
                        reaction = comp.reactions.get_by_id(f'{member}::{rid}')
                        low, high = reaction.bounds
                        reaction.bounds = (min(low, -self.protocol['port_ceiling']), high) if coefficient < 0 else (low, max(high, self.protocol['port_ceiling']))
                for member, rids in blocks.items():
                    for rid in rids: comp.reactions.get_by_id(f'{member}::{rid}').bounds = (0., 0.)
            comp.solver.configuration.tolerances.feasibility = 1e-9
            swiglpk.glp_std_basis(comp.solver.problem)
            value = comp.slim_optimize()
            status = comp.solver.status
            result = {'status': status, 'mu': None, 'blocked_reactions': {k: sorted(v) for k, v in blocks.items()}}
            if status == 'optimal':
                fluxes = {r.id: r.flux for r in comp.reactions}
                if not math.isfinite(value) or any(not math.isfinite(v) for v in fluxes.values()):
                    return {**result, 'status': 'nonfinite_primal'}
                residual = max(abs(c.primal) for c in comp.solver.constraints)
                bound = max(max(r.lower_bound - fluxes[r.id], fluxes[r.id] - r.upper_bound, 0.) for r in comp.reactions)
                result.update(mu=max(0., value), max_balance_residual=residual, max_bound_violation=bound)
                if residual > 1e-6 or bound > 1e-6: result['status'] = 'invalid_primal'
            return result


def agrees(production, oracle):
    if production['status'] != 'completed': return False
    if oracle['status'] == 'infeasible': return production['outcome'] == 'infeasible_with_maintenance'
    return (oracle['status'] == 'optimal' and production.get('mu') is not None
            and math.isfinite(production['mu']) and math.isclose(production['mu'], oracle['mu'], rel_tol=1e-6, abs_tol=1e-6))


def run(root, predictions, out, paths):
    frozen = json.loads((root / 'protocol.json').read_text())
    receipt = json.loads((predictions / 'run.json').read_text())
    if receipt['protocol_sha256'] != sha(root / 'protocol.json') or receipt['predictions_sha256'] != sha(predictions / 'predictions.jsonl'):
        raise ValueError('prediction provenance mismatch')
    if frozen['plans_sha256'] != sha(root / 'plans.jsonl') or frozen['contexts_sha256'] != sha(root / 'contexts.json'):
        raise ValueError('frozen input changed')
    contexts = json.loads((root / 'contexts.json').read_text())
    oracles = {}
    for sp, path in paths.items():
        if sha(path) != frozen['models'][sp]['sha256']: raise ValueError('independent model hash mismatch')
        oracles[sp] = NutrientOracle(path, contexts[sp], frozen['protocol'])
    plans = {(r['case_id'], r['arm']): r for r in map(json.loads, (root / 'plans.jsonl').read_text().splitlines())}
    count, matches, gpr_count, gpr_matches, started, seen = 0, 0, 0, 0, time.time(), set()
    out.mkdir(parents=True, exist_ok=False)
    with (predictions / 'predictions.jsonl').open() as source, (out / 'comparisons.jsonl').open('w') as dest:
        for i, line in enumerate(source):
            row = json.loads(line)
            key = row['case_id'], row['arm']
            if key in seen or key not in plans: raise ValueError('extra/duplicate prediction')
            seen.add(key)
            plan = plans[key]
            if row['plan_sha256'] != digest(plan): raise ValueError('plan identity mismatch')
            oracle = oracles[plan['species']]
            for m in plan['members']:
                equal = oracle.blocked(m['genes']) == m['blocked_reactions']
                gpr_count += 1
                gpr_matches += equal
                if not equal: raise ValueError('independent GPR disagreement: ' + str(key))
            by_kind = {r['id']: r for r in row.get('production', [])}
            for kind in frozen['protocol']['solve_order']:
                expected = oracle.simulate(plan, kind)
                actual = by_kind.get(kind, {'status': 'missing', 'mu': None})
                equal = agrees(actual, expected)
                dest.write(json.dumps({'case_id': key[0], 'arm': key[1], 'kind': kind, 'agrees': equal,
                    'production_status': actual['status'], 'production_mu': actual.get('mu'), 'oracle': expected}, sort_keys=True) + '\n')
                count += 1
                matches += equal
            if (i + 1) % 50 == 0:
                dest.flush()
                print(json.dumps({'cases': i + 1, 'comparisons': count, 'agree': matches, 'elapsed_s': time.time() - started}), flush=True)
    shard_index, shard_count = receipt.get('shard_index', 0), receipt.get('shard_count', 1)
    expected_keys = {k for i, k in enumerate(plans) if i % shard_count == shard_index}
    if seen != expected_keys: raise ValueError('incomplete shard replay')
    report = {'status': 'verified' if matches == count else 'disagreements', 'comparisons': count, 'agree': matches,
        'gpr_comparisons': gpr_count, 'gpr_agree': gpr_matches, 'elapsed_s': time.time() - started,
        'protocol_sha256': sha(root / 'protocol.json'), 'run_sha256': sha(predictions / 'run.json'),
        'comparisons_sha256': sha(out / 'comparisons.jsonl'), 'cobra_version': cobra.__version__,
        'tolerances': {'objective_absolute': 1e-6, 'objective_relative': 1e-6, 'primal_residual_and_bounds': 1e-6, 'GLPK_primal': 1e-9},
        'sources': {str(p): sha(p) for p in [Path(__file__), Path(__file__).with_name('check_community_cobra.py')]}}
    (out / 'verification.json').write_text(json.dumps(report, indent=2) + '\n')
    return report


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('plan', 'run', 'out', 'ec-model', 'bs-model'): p.add_argument('--' + key, type=Path, required=True)
    a = p.parse_args()
    report = run(a.plan, a.run, a.out, {'Ec': a.ec_model, 'Bs': a.bs_model})
    print(json.dumps(report))
    sys.exit(0 if report['status'] == 'verified' else 1)
