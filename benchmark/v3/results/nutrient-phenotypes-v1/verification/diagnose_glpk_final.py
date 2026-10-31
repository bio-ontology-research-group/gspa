"""Fresh GLPK replay of failed quantitative comparisons, without phenotype labels."""
from pathlib import Path
import hashlib
import json
import sys
import cobra

root = Path(sys.argv[1])
previous = Path(sys.argv[2])
plans = {(r['case_id'], r['arm']): r for r in
         (json.loads(line) for line in (root / 'plan/plans.jsonl').read_text().splitlines())}
oracle = json.loads((previous / 'oracle-v2.json').read_text())
rows = []
for case in oracle['results']:
    for comparison in case['comparisons']:
        if comparison['flux_agrees']:
            continue
        plan = plans[case['case_id'], case['arm']]
        manifest = plan['manifest']
        for tolerance in (None, 1e-9, 1e-10):
            model = cobra.io.read_sbml_model(root / 'models/iYO844.xml', f_replace={})
            model.solver = 'glpk'
            config = model.solver.configuration
            if tolerance is not None:
                config.tolerances.feasibility = tolerance
                # Optlang GLPK exposes no public optimality tolerance; this is
                # a diagnostic using the installed GLPK simplex control block.
                config._smcp.tol_dj = tolerance
            for reaction in model.reactions:
                if len(reaction.metabolites) == 1:
                    sign = next(iter(reaction.metabolites.values()))
                    medium = manifest['baseline']['medium']
                    lo, hi = reaction.bounds
                    if sign < 0:
                        lo = -medium[reaction.id] if reaction.id in medium else max(lo, 0.)
                    else:
                        hi = medium[reaction.id] if reaction.id in medium else min(hi, 0.)
                    reaction.bounds = lo, hi
            maintenance = manifest['objectives']['maintenance']
            atpm = model.reactions.get_by_id(maintenance['reaction'])
            atpm.lower_bound = max(atpm.lower_bound, maintenance['lower_bound'])
            if comparison['kind'] == 'knockout':
                for gene in plan['genes']:
                    model.genes.get_by_id(gene).knock_out()
            model.objective = model.reactions.get_by_id(manifest['objectives']['biomass'])
            model.objective_direction = 'max'
            solution = model.optimize()
            rows.append({'case_id': case['case_id'], 'arm': case['arm'], 'kind': comparison['kind'],
                         'tolerance_setting': tolerance, 'tol_bnd': config._smcp.tol_bnd,
                         'tol_dj': config._smcp.tol_dj, 'status': solution.status,
                         'growth': solution.objective_value,
                         'max_balance_residual': max(abs(c.primal) for c in model.constraints),
                         'max_bound_violation': max([0.] + [max(r.lower_bound - r.flux, r.flux - r.upper_bound) for r in model.reactions])})
print(json.dumps({'schema': 'gspa.glpk_numerical_diagnosis/1',
                  'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  'cobra_version': cobra.__version__, 'comparisons': rows,
                  'interpretation': 'Fresh model and tightened primal/dual LP tolerances only; no labels or replacement of registered predictions.'}, indent=2))
