"""Sequence-only comparator Q on the seen labels L0 only (plan section 8): per genome best EC per protein,
call when score >= t; find the largest t with precision >= 0.75 on >= 20 labelled decisions (pooled over the 22
genomes). Reads the six-predictor candidate tables. No held-out label is read here."""
import csv, glob, collections, os, re
R='/mnt/data2/leechuck/gspa'
L0=collections.defaultdict(set)
for r in csv.DictReader(open(R+'/benchmark/manifests/suite6_price149_manifest.tsv'),delimiter='\t'):
    for ec in r['ground_truth_ecs'].split(';'):
        if re.fullmatch(r'\d+\.\d+\.\d+\.\d+',ec): L0[(r['host_assembly'],r['proteome_protein_id'])].add(ec)
for r in csv.DictReader(open(R+'/benchmark/v3/offgrowth/data/swissprot_experimental_ecs_price22.tsv'),delimiter='\t'):
    L0[(r['genome_id'],r['protein_id'])].add(r['function_id'])
print('L0 labelled proteins',len(L0))
dec=[]  # (score, supported) for labelled proteins' best EC
for f in sorted(glob.glob(os.path.expanduser('~/tmp/gspa-v3-offgrowth2/candidates-six/*.candidates.tsv'))):
    g=os.path.basename(f).split('.candidates')[0]
    best={}
    for row in csv.DictReader(open(f),delimiter='\t'):
        p=row['protein_id']; s=float(row['score'])
        if (g,p) in L0 and s>best.get(p,(None,-1))[1]: best[p]=(row['ec'],s)
    for p,(ec,s) in best.items():
        dec.append((s, ec in L0[(g,p)]))
dec.sort(reverse=True)
print('labelled proteins with a candidate', len(dec), 'top-1 exact precision overall', round(sum(d[1] for d in dec)/len(dec),3))
best_t=None
for i in range(20,len(dec)+1):
    t=dec[i-1][0]; sup=sum(d[1] for d in dec[:i]); prec=sup/i
    if prec>=0.75: best_t=(t,i,round(prec,3))
print('fixed-precision (>=0.75 on >=20 decisions): largest threshold with the property =', best_t)
for t in (0.6,0.55,0.5,0.45,0.4,0.3):
    sel=[d for d in dec if d[0]>=t]; print(f'  t>={t}: decisions {len(sel)} precision {sum(d[1] for d in sel)/max(1,len(sel)):.3f}')
