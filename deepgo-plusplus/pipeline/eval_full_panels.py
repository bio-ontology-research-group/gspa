import os, sys
from collections import defaultdict
import numpy as np
_HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,_HERE)
from train_integrator import load_dag, load_aspect_from_dag_ia, load_component, load_gt
from ablation import oof_predictions, write_preds, score_dir_with_cafaeval, log

COMPS="full_comps"; GT="gt/gt_no_cleanA.tsv"; DAG="go-dag.tsv"
IA=os.path.expanduser("~/cafa6_cleaneval/IA.tsv"); OBO=os.path.expanduser("~/cafa6_cleaneval/go.obo")
TT="train_terms.tsv"; OUT="full_panels_out"

HOM=['diam','foldseek','interpro','net','clean','lit']
LEAN=['diam','interpro','cnn','net_union','esm2_knn','proteinfer']
PANELS=[
 ('cpu_lean (6, CPU-only)', LEAN),
 ('cpu_lean + esm2_35m head (7)', LEAN+['esm2_head']),
 ('full DG++ BCE heads', HOM+['esm2_650m_bce','prostt5_bce','esm2_3b_bce']),
 ('full DG++ MCM heads', HOM+['esm2_650m_mcm','prostt5_mcm','esm2_3b_mcm']),
 ('full MCM + cpu aux (kitchen sink)', HOM+['esm2_650m_mcm','prostt5_mcm','esm2_3b_mcm','cnn','esm2_knn','proteinfer','esm2_head','interpro_lr']),
]
os.makedirs(OUT,exist_ok=True)
anc=load_dag(DAG); aspect_of=load_aspect_from_dag_ia(IA,anc)
ia={}
for line in open(IA):
    g,_,v=line.rstrip('\n').partition('\t')
    try: ia[g]=float(v)
    except ValueError: pass
freq=defaultdict(int)
fh=open(TT); next(fh,None)
for line in fh:
    p=line.rstrip('\n').split('\t')
    if len(p)>=2: freq[p[1]]+=1
logfreq={t:np.log10(c+1) for t,c in freq.items()}
gt=load_gt(GT,anc,aspect_of); proteins=sorted(gt.keys()); keep=set(proteins)
log(f'GT proteins: {len(proteins)}')
need=sorted(set(c for _,cl in PANELS for c in cl))
comp={}
for c in need:
    p=os.path.join(COMPS,f'{c}.tsv.gz')
    if not os.path.exists(p): p=os.path.join(COMPS,f'{c}.tsv')
    if not os.path.exists(p): log(f'WARN missing {c}'); continue
    comp[c]=load_component(p,keep,anc,aspect_of)
pd=os.path.join(OUT,'P'); os.makedirs(pd,exist_ok=True)
f2c={}
for i,(label,cl) in enumerate(PANELS):
    cl=[c for c in cl if c in comp]
    log(f'OOF [{i+1}/{len(PANELS)}] {label}')
    rows=oof_predictions(comp,proteins,gt,aspect_of,cl,'logreg','scores',ia,logfreq,folds=5)
    fn=f'p{i:02d}.tsv'; write_preds(rows,os.path.join(pd,fn)); f2c[fn[:-4]]=label
log('cafaeval ...')
sc=score_dir_with_cafaeval(pd,os.path.join(OUT,'res'),OBO,GT,IA)
print("\n===== Ensemble: full DG++ (BCE vs MCM heads) vs CPU-only cpu_lean (no-knowledge f_w) =====")
print("%-38s %6s  %5s %5s %5s"%("panel","mean","MF","BP","CC"))
print("-"*64)
res=[(f2c[k],)+sc[k] for k in f2c if k in sc]
for label,mf,bp,cc,mean in [(l,)+tuple(v) for l,v in [(f2c[k],sc[k]) for k in sorted(f2c)]]:
    print("%-38s %6.3f  %5.3f %5.3f %5.3f"%(label,mean,mf,bp,cc))
