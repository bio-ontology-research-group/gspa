"""Numerical diagnosis only; never replaces the frozen benchmark prediction."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
import scipy
from scipy.optimize import linprog

root = Path(sys.argv[1])
sys.path.insert(0, str(root / 'scripts'))
import metabolic_experiments as ME

case_id, arm = 'NUTRIENT2026:Bs:GNAc + NoC:M+Met', 'additional_inorganic:oxygen=20'
plans = [json.loads(line) for line in (root / 'plan/plans.jsonl').read_text().splitlines()]
plan = next(r for r in plans if r['case_id'] == case_id and r['arm'] == arm)
model = ME.load_sbml_model(root / 'models/iYO844.xml')
manifest = plan['manifest']
run = ME.Runner(model, manifest, ME.resolve_config(manifest))
menu, rejected = ME.build_menu(model, manifest, 0, run.external)
assert len(menu) == 1 and not rejected
lb, ub = run.base_lb.copy(), run.base_ub.copy()
killed = []  # Diagnosing the matched baseline; no gene deletion.
c = np.zeros(len(model.rxn_ids))
c[run.bio] = -1.
results = []
for method, presolve in [('highs', True), ('highs', False), ('highs-ds', True), ('highs-ipm', True)]:
    result = linprog(c, A_eq=model.S, b_eq=np.zeros(model.S.shape[0]), bounds=list(zip(lb, ub)),
                     method=method, options={'presolve': presolve})
    row = {'method': method, 'presolve': presolve, 'status': result.status,
           'message': result.message, 'success': result.success,
           'objective': float(result.fun) if result.fun is not None else None}
    if result.x is not None:
        row['max_balance_residual'] = float(np.max(np.abs(model.S @ result.x)))
        row['max_bound_violation'] = float(max(np.max(lb - result.x), np.max(result.x - ub), 0.))
    results.append(row)
print(json.dumps({'case_id': case_id, 'arm': arm, 'model_sha256': model.sha256,
                  'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  'scipy_version': scipy.__version__, 'knocked_out_reactions': killed,
                  'plan_sha256': hashlib.sha256(json.dumps(plan, sort_keys=True, allow_nan=False).encode()).hexdigest(),
                  'interpretation': 'Post-run solver diagnosis; no observation read and no benchmark prediction replaced.',
                  'results': results}, indent=2))
