"""Independent check on real data that the N1 gate behaves as the plan prescribes."""
import gzip, json, sys
sys.path.insert(0, '/mnt/data2/leechuck/gspa/benchmark/v3')
from offgrowth2.rules import reaction_calls

path = '/home/leechuck/tmp/gspa-v3-offgrowth2/score-tierA-six/staging/GCF_000146165.2/genome_results2.json.gz'
d = json.load(gzip.open(path))
draft = {r for r in d['draft_reactions']}
run = [r for r in d['runs'] if r['run_kind'] == 'unmasked'][0]
calls = reaction_calls(run, fill_max=3, require_stable=True, test_filter=lambda t: True, draft_reactions=draft)
base = set(run['baseline_fill']['fill'])
in_draft = [r for r in calls if r in draft]
in_base = [r for r in calls if r in base]
bad = []
for r, c in calls.items():
    ok = False
    for t in run['tests']:
        if r in set(t['stable_fill']) and len(t['stable_fill']) <= 3:
            ok = True
            break
    if not ok:
        bad.append(r)
print('called reactions', len(calls))
print('called although in the draft model:', len(in_draft))
print('called although in the baseline fill:', len(in_base))
print('called without a stable fill of at most 3:', len(bad))
ungated = reaction_calls(run, fill_max=float('inf'), require_stable=False, test_filter=lambda t: True, draft_reactions=draft)
print('reactions any test adds (no gate):', len(ungated), '-> gate keeps', round(100 * len(calls) / max(1, len(ungated))), '%')
