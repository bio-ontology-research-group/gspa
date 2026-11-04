#!/usr/bin/env python3
"""Room for improvement on the maskable loci: protein-centric and reaction-centric ranks at stage 0 / raw,
score windows, and split membership. Design-time only."""
import csv, json, os, re
from collections import defaultdict, Counter
ROOT = os.path.expanduser('~/tmp/row5-design')
REPO = '/mnt/data2/leechuck/gspa'
per = json.load(open(os.path.join(ROOT, 'per_locus.json')))
mask = [tuple(k.split('|')) for k, v in per.items() if v['maskable']]
price22 = [l.strip() for l in open(os.path.join(ROOT, 'price22.txt')) if l.strip()]
cohort18 = [l.strip() for l in open(os.path.join(ROOT, 'cohort18.txt')) if l.strip()]

# splits
csplit = {}
with open(os.path.join(REPO, 'benchmark/v3/results/training-cohort/splits.tsv')) as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        g = r.get('genome_id') or r.get('genome') or r.get('assembly')
        sp = r.get('split') or r.get('role') or r.get('arm')
        csplit[g] = sp
ps = json.load(open(os.path.join(REPO, 'benchmark/v3/results/price-calibration-context-v1/row1/split.json')))
psplit = {g: ('price_' + a) for g, a in ps['assignment'].items()}
print('split membership of maskable loci:')
c = Counter()
for s, g, p, ec in mask:
    c[(s, csplit.get(g) or psplit.get(g) or '?')] += 1
for k, v in sorted(c.items()): print('  ', k, v)

def load_claims_ec(path):
    d = defaultdict(dict)
    with open(path) as fh:
        for line in fh:
            j = json.loads(line)
            if j.get('function_type') == 'EC':
                d[j['protein_id']][j['function_id']] = max(j.get('raw_score', 0.0), d[j['protein_id']].get(j['function_id'], 0.0))
    return d
sets = {}
for regime in ('homology', 'family_excluded'):
    sets['cohort_diamond_' + regime] = {g: load_claims_ec(os.path.expanduser(f'~/.cache/gspa-v3-data/training-cohort/baseline/claims/{regime}/{g}.claims.jsonl')) for g in cohort18}
sets['price_diamond_contract'] = {g: load_claims_ec(f'/mnt/data2/leechuck/gspa-v3-data/price-contract-learned-v1/claims-run/claims/{g}.claims.jsonl') for g in price22 if os.path.exists(f'/mnt/data2/leechuck/gspa-v3-data/price-contract-learned-v1/claims-run/claims/{g}.claims.jsonl')}

print('\nreaction-centric view on maskable loci (rank of P among the genome proteins claiming E, raw score):')
for name, cs in sets.items():
    n = top = sub = none = 0; window = 0; comp = []
    for s, g, p, ec in mask:
        if g not in cs: continue
        n += 1
        claimers = [(pp, d[ec]) for pp, d in cs[g].items() if ec in d]
        if not any(pp == p for pp, _ in claimers): none += 1; continue
        claimers.sort(key=lambda x: -x[1])
        rank = [pp for pp, _ in claimers].index(p) + 1
        comp.append(len(claimers))
        if rank == 1: top += 1
        else: sub += 1
        sc = cs[g][p][ec]
        if 0.18 <= sc <= 0.5: window += 1
    print(f'  {name:30s} loci={n} P_top_for_E={top} P_not_top={sub} P_no_claim={none} '
          f'P_score_in_prior_window[0.18,0.5]={window} median_claimers_per_E={sorted(comp)[len(comp)//2] if comp else None}')

# six Price predictors: scores of the reference proteins only (stage0); window membership of the true EC
reg = json.load(open('/mnt/data2/leechuck/gspa-v3-data/price-calibration-context-v1/inputs/rc1-evaluation/registration.json'))
pred = defaultdict(dict)
for cell in reg['cells']:
    f = f"/mnt/data2/leechuck/gspa-v3-data/price-calibration-context-v1/inputs/rc1-evaluation/retained/{cell['index']}.json"
    d = json.load(open(f))
    pred[cell['predictor']][cell['assembly']] = d['scores']['stage0']
print('\nsix Price predictors, maskable Price loci: true-EC stage-0 score when it is a candidate but not top-1:')
for name, cs in pred.items():
    rows = []
    for s, g, p, ec in mask:
        if g not in cs: continue
        c = cs[g].get(p) or {}
        if ec in c and c[ec] < max(c.values()):
            rows.append((round(c[ec], 2), round(max(c.values()), 2)))
    win = sum(1 for a, b in rows if 0.18 <= a <= 0.5)
    print(f'  {name:10s} n={len(rows):2d} in_window={win:2d} (true, top) pairs: {sorted(rows)}')
