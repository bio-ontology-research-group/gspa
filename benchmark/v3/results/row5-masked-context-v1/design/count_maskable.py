#!/usr/bin/env python3
"""Feasibility count for the row-5 masking design (design-time only, no scores read for any arm).

For every experimentally labelled (genome, protein, EC) on the 40 genomes that have a gapsmith model
(22 Price + 18 cohort), classify what happens to the reaction(s) carrying that EC in the gapsmith
pathway-search tables if the protein's hits are removed ("masked"):

  already_gap         pathway predicted, reaction group has no good_blast row (gap without masking)
  maskable_gap        pathway predicted, every good_blast row of the group names the labelled protein
                      -> removing the protein turns the group into a gap under PathwayGapBuilder's rules
  redundant           pathway predicted, another protein also has a good_blast row -> stays strong_hit_present
  pathway_not_predicted
  ec_not_in_tables    the EC appears in no reaction row of that genome's table
  excluded_rule       spontaneous / complex / reference_unknown exclusions of PathwayGapBuilder

A locus is "maskable" when at least one (pathway, reaction) group is maskable_gap or already_gap AND the
pathway has at least one other reaction group with a good_blast row (context present).
Candidate availability is then reported per predictor set.
"""
import csv, glob, json, os, re, sys
from collections import defaultdict, Counter

ROOT = os.path.expanduser('~/tmp/row5-design')
TABLES = os.path.join(ROOT, 'tables')
REPO = '/mnt/data2/leechuck/gspa'
EC_RE = re.compile(r'^\d+\.\d+\.\d+\.\d+$')

price22 = [l.strip() for l in open(os.path.join(ROOT, 'price22.txt')) if l.strip()]
cohort18 = [l.strip() for l in open(os.path.join(ROOT, 'cohort18.txt')) if l.strip()]

# ---------------------------------------------------------------- labels
labels = []  # (set, genome, protein, ec)
with open(os.path.join(REPO, 'benchmark/manifests/suite6_price149_manifest.tsv')) as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        for ec in r['ground_truth_ecs'].split(';'):
            ec = ec.strip()
            if EC_RE.match(ec):
                labels.append(('price148', r['host_assembly'], r['proteome_protein_id'], ec))
with open(os.path.join(REPO, 'benchmark/v3/offgrowth/data/swissprot_experimental_ecs_price22.tsv')) as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        if EC_RE.match(r['function_id']):
            labels.append(('price_swissprot', r['genome_id'], r['protein_id'], r['function_id']))
with open(os.path.expanduser('~/.cache/gspa-v3-data/training-cohort/labels/labels.tsv')) as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        fid = r['function_id'].replace('EC:', '')
        if r['tier'] == 'experimental' and EC_RE.match(fid):
            labels.append(('cohort_swissprot', r['genome_id'], r['protein_id'], fid))
labels = sorted(set(labels))
print('labelled (set, genome, protein, EC) rows:', len(labels), Counter(l[0] for l in labels))

# ---------------------------------------------------------------- tables
def load_tables(g):
    d = os.path.join(TABLES, g)
    rx = os.path.join(d, f'{g}-all-Reactions.tbl')
    pw = os.path.join(d, f'{g}-all-Pathways.tbl')
    if not (os.path.exists(rx) and os.path.exists(pw)):
        return None
    pred = {}
    conflict = set()
    with open(pw) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            pid = r['ID'].strip('|')
            p = r['Prediction'].strip().upper() == 'TRUE'
            if pid in pred and pred[pid] != p:
                conflict.add(pid)
            pred[pid] = p
    for c in conflict:  # row-5 converter rule: conflicting pathway ids are dropped (not predicted)
        pred[c] = False
    groups = defaultdict(list)  # (pathway, rxn) -> rows
    with open(rx) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            pw_id = r['pathway'].strip('|'); rxn = r['rxn'].strip('|')
            prot = (r.get('stitle') or '').split(' ')[0].strip()
            ecs = [e.replace('EC:', '') for e in re.split(r'[/,;\s]+', (r.get('ec') or '').strip()) if e and e not in ('NA', '-')]
            groups[(pw_id, rxn)].append(dict(status=r['status'], prot=prot, ecs=ecs,
                                             spont=(r.get('spont') or '').upper() == 'TRUE',
                                             complex=(r.get('is_complex') or '').upper() == 'TRUE'))
    return pred, groups

def group_state(rows, masked_prot=None):
    rows_eff = [dict(r, status=('no_blast' if (r['prot'] == masked_prot and r['status'] == 'good_blast') else r['status'])) for r in rows]
    if any(r['status'] == 'good_blast' for r in rows_eff): return 'strong_hit_present'
    if any(r['status'] == 'no_seq_data' for r in rows_eff): return 'reference_unknown'
    if any(r['spont'] or r['status'] == 'spontaneous' for r in rows_eff): return 'spontaneous'
    if any(r['complex'] for r in rows_eff): return 'complex_requires_subunit_inference'
    return 'gap'

tables = {}
for g in price22 + cohort18:
    t = load_tables(g)
    if t is None: print('MISSING tables', g); continue
    tables[g] = t

# ---------------------------------------------------------------- classify
per_locus = {}
cls_counter = Counter()
detail_rows = []
for s, g, p, ec in labels:
    if g not in tables:
        cls_counter['no_tables'] += 1; continue
    pred, groups = tables[g]
    hit_groups = [(k, rows) for k, rows in groups.items() if any(ec in r['ecs'] for r in rows)]
    hit_by_prot = [(k, rows) for k, rows in groups.items() if any(r['prot'] == p and r['status'] == 'good_blast' for r in rows)]
    if not hit_groups:
        cls = 'ec_not_in_tables'
        cls_counter[cls] += 1
        per_locus[(s, g, p, ec)] = dict(cls=cls, maskable=False, n_pred_groups=0, protein_hit_groups=len(hit_by_prot))
        continue
    states = Counter()
    maskable_groups = []; context_ok = 0
    for (pw_id, rxn), rows in hit_groups:
        if not pred.get(pw_id, False):
            states['pathway_not_predicted'] += 1; continue
        before = group_state(rows)
        after = group_state(rows, masked_prot=p)
        p_present = any(r['prot'] == p and r['status'] == 'good_blast' for r in rows)
        if before == 'gap': st = 'already_gap'
        elif before == 'strong_hit_present' and after == 'gap': st = 'maskable_gap'
        elif before == 'strong_hit_present' and after == 'strong_hit_present':
            st = 'redundant' if p_present else 'other_protein_only'
        else: st = 'excluded_' + after
        states[st] += 1
        if st in ('already_gap', 'maskable_gap'):
            # context: other reaction groups of this pathway with a strong hit from another protein
            others = [k2 for k2 in groups if k2[0] == pw_id and k2[1] != rxn and
                      any(r['status'] == 'good_blast' and r['prot'] != p for r in groups[k2])]
            maskable_groups.append((pw_id, rxn, st, len(others)))
            if others: context_ok += 1
    if any(st in ('already_gap', 'maskable_gap') for st in states):
        cls = 'maskable_with_context' if context_ok else 'maskable_no_context'
    elif 'redundant' in states or 'other_protein_only' in states:
        cls = 'redundant'
    elif 'pathway_not_predicted' in states and len(states) == states['pathway_not_predicted']:
        cls = 'pathway_not_predicted'
    else:
        cls = 'excluded_rule'
    cls_counter[cls] += 1
    per_locus[(s, g, p, ec)] = dict(cls=cls, maskable=cls == 'maskable_with_context', states=dict(states),
                                    groups=maskable_groups, protein_hit_groups=len(hit_by_prot))

print('\nclassification of labelled loci:')
for k, v in cls_counter.most_common(): print(f'  {k:28s} {v}')
by_set = defaultdict(Counter)
for (s, g, p, ec), v in per_locus.items(): by_set[s][v['cls']] += 1
print('\nby label set:')
for s, c in by_set.items(): print(' ', s, dict(c))
print('\nmaskable_with_context per genome:')
pg = Counter(g for (s, g, p, ec), v in per_locus.items() if v['maskable'])
for g, n in sorted(pg.items(), key=lambda x: -x[1]): print(f'  {g} {n}')
print('  genomes with >=1:', len(pg))

# ---------------------------------------------------------------- candidates
# DIAMOND Swiss-Prot claims: cohort (homology, family_excluded) and Price (other session, read-only)
def load_claims_ec(path):
    d = defaultdict(dict)
    with open(path) as fh:
        for line in fh:
            try: j = json.loads(line)
            except Exception: continue
            if j.get('function_type') == 'EC':
                d[j['protein_id']][j['function_id']] = max(j.get('raw_score', 0.0), d[j['protein_id']].get(j['function_id'], 0.0))
    return d
cand_sets = {}
for regime in ('homology', 'family_excluded'):
    cand_sets['cohort_diamond_' + regime] = {}
    for g in cohort18:
        f = os.path.expanduser(f'~/.cache/gspa-v3-data/training-cohort/baseline/claims/{regime}/{g}.claims.jsonl')
        if os.path.exists(f): cand_sets['cohort_diamond_' + regime][g] = load_claims_ec(f)
cand_sets['price_diamond_contract'] = {}
for g in price22:
    f = f'/mnt/data2/leechuck/gspa-v3-data/price-contract-learned-v1/claims-run/claims/{g}.claims.jsonl'
    if os.path.exists(f): cand_sets['price_diamond_contract'][g] = load_claims_ec(f)
# six Price predictors: rc1 retained raw_export scores of the reference proteins
reg = json.load(open('/mnt/data2/leechuck/gspa-v3-data/price-calibration-context-v1/inputs/rc1-evaluation/registration.json'))
for cell in reg['cells']:
    f = f"/mnt/data2/leechuck/gspa-v3-data/price-calibration-context-v1/inputs/rc1-evaluation/retained/{cell['index']}.json"
    if not os.path.exists(f): continue
    d = json.load(open(f))
    key = 'price_' + cell['predictor']
    cand_sets.setdefault(key, {})[cell['assembly']] = d['scores'].get('stage0') or d['scores'].get('raw_export')

print('\ncandidate availability on maskable_with_context loci (E among candidates of P; E top-1 at stage 0 / raw):')
for name, cs in cand_sets.items():
    n = has = top1 = 0; sub = 0
    for (s, g, p, ec), v in per_locus.items():
        if not v['maskable'] or g not in cs: continue
        n += 1
        c = cs[g].get(p) or {}
        if ec in c:
            has += 1
            best = max(c.values())
            if c[ec] >= best: top1 += 1
            else: sub += 1
    print(f'  {name:32s} loci={n:4d} has_candidate={has:4d} already_top1={top1:4d} candidate_but_not_top1={sub:4d}')

json.dump({'%s|%s|%s|%s' % k: v for k, v in per_locus.items()}, open(os.path.join(ROOT, 'per_locus.json'), 'w'), indent=1)
