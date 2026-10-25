#!/usr/bin/env python3
"""Live-pipeline check on E. coli K-12 MG1655, run inside a container of the deployed service image.

annotate twice on the example's GFF3 + genome: (a) as the service does now (--min-score 0.5: every GO prediction
below 0.5 is dropped before integration) and (b) without the global floor (DeepGO++ Light keeps its own 0.5 floor
through the config; DIAMOND GO transfers from 30 % identity up enter as claims with score = identity). Each claim set
(+ the example job's PSORTb claims) is integrated with all priors and the module definitions in module mode, and
with the coherence prior off. GO annotations at posterior >= 0.50 are judged against the curated E. coli annotation
(truth_all / truth_exp, propagated; RefSeq -> UniProt map). Writes live_summary.tsv and the integrated outputs.
"""
import csv, json, os, subprocess, sys, time
from collections import defaultdict
W = sys.argv[1]; JOB = sys.argv[2]; R = '/opt/gspa-resources'; A = '/opt/dgpp_assets'; BIN = '/opt/gspa/bin/gspa-cli'
os.makedirs(W, exist_ok=True)

def sh(cmd, log):
    t0 = time.time()
    with open(log, 'w') as fh: rc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT).returncode
    print(f'{os.path.basename(log)} rc={rc} {time.time()-t0:.0f}s', flush=True)
    if rc: raise SystemExit(rc)

config = {'predictors': {'domains': {'enabled': False}, 'similarity': {'enabled': True, 'database': f'{R}/diamond/swissprot.dmnd', 'assignProductNames': True},
          'operons': {'enabled': True}, 'neural': {'deepGoPlusPlusLight': {'minScore': 0.5}}}, 'integration': {'enabled': False}, 'quality': {'enabled': False}}
json.dump(config, open(f'{W}/config.json', 'w'))
aspects = {}
for line in open(f'{R}/go/go_aspect_map.tsv') if os.path.exists(f'{R}/go/go_aspect_map.tsv') else []:
    f = line.rstrip('\n').split('\t')
    if len(f) >= 2: aspects[f[0]] = f[1]
psort = [l for l in open(f'{JOB}/claims.jsonl') if '"psortb"' in l]
for arm, extra in (('floor', ['--min-score', '0.5']), ('keep', [])):
    out = f'{W}/{arm}/annot'
    if not os.path.exists(f'{W}/{arm}/claims.jsonl'):
        sh([BIN, 'annotate', '--genome', f'{JOB}/sample.fna', '--gff3', f'{JOB}/sample.gff3', '-o', out, '-t', '4', '-c', f'{W}/config.json',
            '--base-predictor', 'light', '--neural-sidecar', '/app/benchmark/neural/run_neural_predictors.py', '--deepgo-plusplus-light-assets', A,
            '--preset', 'bacteria'] + extra, f'{W}/{arm}_annotate.log')
        ann = next(p for p in os.listdir(out) if p.endswith('_annotations.tsv'))
        seen = set()
        with open(f'{W}/{arm}/claims.jsonl', 'w') as fh:
            for row in csv.DictReader(open(f'{out}/{ann}'), delimiter='\t'):
                key = row['protein_id'], row['type'], row['value'], row['source']
                if key in seen: continue
                seen.add(key); s = float(row['score'])
                fh.write(json.dumps(dict(protein_id=key[0], function_type=key[1], function_id=key[2], source=key[3], raw_score=s,
                         calibrated_prob=min(1-1e-6, max(1e-6, s)), go_aspect=aspects.get(key[2], row.get('aspect', '')))) + '\n')
            fh.writelines(psort)
    for v, pri in (('off', 'essentiality,consistency,anatomical,gap_filling,genomic_context'),
                   ('module', 'essentiality,coherence,consistency,anatomical,gap_filling,genomic_context')):
        o = f'{W}/{arm}/{v}_integrated.tsv'
        if os.path.exists(o): continue
        sh([BIN, 'integrate', '--claims', f'{W}/{arm}/claims.jsonl', '--out', o, '--go-owl', f'{R}/go/go.owl', '--ec2go', f'{R}/ec2go/ec2go',
            '--preset', 'bacteria', '--enable-priors', pri, '--pathways', f'{R}/kegg_pathways.tsv', '--modules', f'{R}/kegg_modules.tsv',
            '--module-defs', f'{R}/kegg_module_defs.json', '--module-exclude-organisms', 'eco,ece', '--operons', f'{out}/operons_for_integrate.tsv'],
           f'{W}/{arm}_{v}_integrate.log')

# evaluation
parents, alt, cur = defaultdict(set), {}, None
for line in open(f'{A}/go.obo'):
    line = line.rstrip('\n')
    if line == '[Term]': cur = None
    elif line.startswith('id: GO:'): cur = line[4:]
    elif cur and line.startswith('is_a: '): parents[cur].add(line[6:].split(' ')[0])
    elif cur and line.startswith('relationship: part_of '): parents[cur].add(line.split(' ')[2])
    elif cur and line.startswith('alt_id: '): alt[line[8:]] = cur
cache = {}
def anc(t):
    t = alt.get(t, t)
    if t in cache: return cache[t]
    o, st = {t}, [t]
    while st:
        u = st.pop()
        for p in parents.get(u, ()):
            if p not in o: o.add(p); st.append(p)
    cache[t] = o; return o
m = {}
for line in open(f'{W}/ecoli.refseq_to_uniprot.tsv'):
    a, b = line.rstrip('\n').split('\t')[:2]; m[a] = b
truth = {}
for k in ('all', 'exp'):
    d = defaultdict(set)
    for r in csv.DictReader(open(f'{W}/ecoli_truth_{k}.tsv'), delimiter='\t'): d[r['accession']].add(r['go_term'])
    truth[k] = {p: set().union(*(anc(t) for t in ts)) for p, ts in d.items()}
def pairs(path):
    s = set()
    for x in csv.DictReader(open(path), delimiter='\t'):
        if x['type'] == 'GO' and float(x['posterior_prob']) >= 0.5 and x['protein_id'] in m: s.add((m[x['protein_id']], x['function_id']))
    return s
rows = []
for arm in ('floor', 'keep'):
    P = {v: pairs(f'{W}/{arm}/{v}_integrated.tsv') for v in ('off', 'module')}
    for v in ('off', 'module'):
        for k in ('all', 'exp'):
            T = truth[k]; sel = [p for p in P[v] if p[0] in T]
            tp = sum(1 for p, g in sel if g in T[p]); fp = len(sel) - tp
            tot = sum(len(ts) for p, ts in T.items())
            covered = set()
            for p, g in sel:
                if g in T[p]: covered.add((p, g))
            rows.append(dict(arm=arm, variant=v, truth=k, annotations=len(sel), tp=tp, fp=fp, precision=tp / max(1, len(sel))))
    add = P['module'] - P['off']; rem = P['off'] - P['module']
    for k in ('all', 'exp'):
        T = truth[k]; a = [p for p in add if p[0] in T]
        rows.append(dict(arm=arm, variant='module_minus_off', truth=k, annotations=len(a), tp=sum(1 for p, g in a if g in T[p]), fp=sum(1 for p, g in a if g not in T[p]),
                         precision=sum(1 for p, g in a if g in T[p]) / max(1, len(a)), removed=len([p for p in rem if p[0] in T])))
    resc = f'{W}/{arm}/module_integrated.tsv.rescues.tsv'
    n = sum(1 for _ in open(resc)) - 1 if os.path.exists(resc) else 0
    mj = json.load(open(f'{W}/{arm}/module_integrated.tsv.modules.json'))['summary']
    rows.append(dict(arm=arm, variant='module_report', truth='', annotations=n, tp='', fp='', precision='', removed=json.dumps(mj)))
with open(f'{W}/live_summary.tsv', 'w') as fh:
    cols = ['arm', 'variant', 'truth', 'annotations', 'tp', 'fp', 'precision', 'removed']
    fh.write('\t'.join(cols) + '\n')
    for r in rows: fh.write('\t'.join(str(r.get(c, '')) for c in cols) + '\n')
print(open(f'{W}/live_summary.tsv').read())
