"""Report fixed constant-class baselines without changing any prediction."""
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path.cwd() / 'benchmark/v3'))
import run_nutrient_phenotypes as N

root = Path(sys.argv[1])
plans = {(r['case_id'], r['arm']): r for r in N.read_rows(root / 'plan/plans.jsonl')}
evaluation = json.loads((root / 'evaluation.json').read_text())
groups, pairs = defaultdict(list), defaultdict(dict)
for row in evaluation['cases']:
    plan = plans[row['case_id'], row['arm']]
    scope = 'primary' if plan['primary_environment'] else 'no_nitrogen_control'
    key = row['arm'], scope, plan['species']
    groups[key].append(row['observed_growth'])
    if plan['kind'] != 'wild_type':
        pairs[key + (plan['strain_id'], plan['environment'])][plan['kind']] = row['observed_growth']
rescue = defaultdict(list)
for key, assays in pairs.items():
    a, b = assays.get('monoculture'), assays.get('supplemented_monoculture')
    rescue[key[:3]].append(None if a is None or b is None else b and not a)

def baseline(table):
    return [{'arm': k[0], 'scope': k[1], 'species': k[2], 'constant_prediction': value,
             **N.confusion([(observed, value) for observed in outcomes])}
            for k, outcomes in sorted(table.items()) for value in (False, True)]

report = {'schema': 'gspa.nutrient_constant_baselines/1', 'evaluation_sha256': N.sha(root / 'evaluation.json'),
          'script_sha256': N.sha(Path(__file__)), 'assays': baseline(groups), 'paired_rescue': baseline(rescue),
          'interpretation': 'Descriptive constant-class comparators on the same full catalogue; no fitted parameter, prediction replacement or arm selection.'}
out = root / 'constant-baselines.json'
if out.exists():
    raise ValueError('baseline report already exists')
N.write_json(out, report)
print(json.dumps([r for r in report['paired_rescue'] if r['scope'] == 'primary' and r['constant_prediction'] is True], indent=2))
