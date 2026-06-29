#!/bin/bash
# Reproduce the 2026-06-23 leak-free no-knowledge evaluation.
#
# Builds the leak-free GT (fixing the GAF last-modified-date contamination, see
# deepgo-plusplus/pipeline/build_clean_gt.py + TRAINING.md §1.0a) and scores the
# 2025 submission, the net/mlp components, and the 6+net integrator with the
# official cafaeval on the DIRTY vs CLEAN ground truths.
#
# cafaeval is CPU-heavy -> run on ws (56 cores), per project policy. This script
# builds the small inputs locally, ships them to ws, and runs cafaeval there.
set -eo pipefail
cd "$(dirname "$0")"
CAFA=~/Public/software/cafa6
DAG="$CAFA/go-dag.tsv"; TT="$CAFA/kaggle-official/train_terms.tsv"
REMOTE=ws; RDIR='~/cafa6_cleaneval'

echo "[1] build leak-free GTs (per-aspect train_terms filter + novel-protein subset)"
python3 ../../../deepgo-plusplus/pipeline/build_clean_gt.py \
  --gt-no gt/gt_no.tsv --train-terms "$TT" --obo "$CAFA/go.obo" \
  --out gt/gt_no_cleanA.tsv --novel-proteins-out gt/gt_no_cleanB.tsv

echo "[2] filter predictions to the no-knowledge proteins (small, for transfer)"
python3 - <<'PY'
import gzip
nk={l.split("\t")[0] for l in open("gt/gt_no.tsv") if l.strip()}
import os; os.makedirs("cleaneval_preds",exist_ok=True)
def filt(src,dst,op=open,cols3=False):
    with op(src,'rt') as fh, open(dst,'w') as out:
        for l in fh:
            c=l.rstrip("\n").split("\t")
            if len(c)<3 or c[0] not in nk or not c[1].startswith("GO:"): continue
            out.write(f"{c[0]}\t{c[1]}\t{c[2]}\n")
filt("preds/hoehndorflab.tsv","cleaneval_preds/submission.tsv")
filt("run_7comp_net/preds/ltr.tsv","cleaneval_preds/integrator_6net.tsv")
filt("components/net.tsv.gz","cleaneval_preds/net.tsv",gzip.open)
filt("components/mlp.tsv.gz","cleaneval_preds/mlp.tsv",gzip.open)
PY

echo "[3] ship to $REMOTE and run cafaeval there"
ssh $REMOTE "mkdir -p $RDIR/{gt,P_submission,P_integrator,P_net,P_mlp}"
scp -q "$CAFA/go.obo" "$CAFA/kaggle-official/IA.tsv" $REMOTE:$RDIR/
scp -q gt/gt_no.tsv gt/gt_no_cleanA.tsv gt/gt_no_cleanB.tsv $REMOTE:$RDIR/gt/
scp -q cleaneval_preds/submission.tsv      $REMOTE:$RDIR/P_submission/
scp -q cleaneval_preds/integrator_6net.tsv $REMOTE:$RDIR/P_integrator/
scp -q cleaneval_preds/net.tsv             $REMOTE:$RDIR/P_net/
scp -q cleaneval_preds/mlp.tsv             $REMOTE:$RDIR/P_mlp/

ssh $REMOTE 'cd ~/cafa6_cleaneval; CE=~/cnn_work/.venv/bin/cafaeval; PY=~/cnn_work/.venv/bin/python
declare -A PD=( [submission]=P_submission [integrator]=P_integrator [net]=P_net [mlp]=P_mlp )
printf "%-12s %-12s %s\n" PRED GT "MF    BP    CC    mean (IA-weighted f_w, no-knowledge)"
for name in submission net mlp integrator; do for gt in gt_no gt_no_cleanA gt_no_cleanB; do
  out=res_${name}_${gt}
  $CE go.obo ${PD[$name]} gt/${gt}.tsv -ia IA.tsv -norm cafa -prop max -th_step 0.01 -out_dir $out >/dev/null 2>&1
  $PY - $out/evaluation_best_f_w.tsv $name $gt <<PYEOF
import sys
rows=open(sys.argv[1]).read().splitlines();h=rows[0].split("\t");fi=h.index("f_w");ni=h.index("ns")
m={r.split("\t")[ni]:float(r.split("\t")[fi]) for r in rows[1:]}
mf=m.get("molecular_function",0);bp=m.get("biological_process",0);cc=m.get("cellular_component",0)
print("%-12s %-12s %.3f %.3f %.3f %.3f"%(sys.argv[2],sys.argv[3],mf,bp,cc,(mf+bp+cc)/3))
PYEOF
done; done'

echo "[4] FULL standalone-component ablation (every component incl. cnn, dirty vs clean)"
python3 - <<'PY'
import gzip, os, glob
nk={l.split("\t")[0] for l in open("gt/gt_no.tsv") if l.strip()}
os.makedirs("cleaneval_comps", exist_ok=True)
for src in sorted(glob.glob("components/*.tsv.gz")):
    name=os.path.basename(src)[:-7]
    with gzip.open(src,'rt') as fh, open(f"cleaneval_comps/{name}.tsv","w") as out:
        for l in fh:
            c=l.rstrip("\n").split("\t")
            if len(c)>=3 and c[0] in nk and c[1].startswith("GO:"): out.write(f"{c[0]}\t{c[1]}\t{c[2]}\n")
PY
ssh $REMOTE "rm -rf $RDIR/P_comps; mkdir -p $RDIR/P_comps"
scp -q cleaneval_comps/*.tsv $REMOTE:$RDIR/P_comps/
ssh $REMOTE 'cd ~/cafa6_cleaneval; CE=~/cnn_work/.venv/bin/cafaeval; PY=~/cnn_work/.venv/bin/python
for gt in gt_no gt_no_cleanA; do $CE go.obo P_comps gt/${gt}.tsv -ia IA.tsv -norm cafa -prop max -th_step 0.01 -out_dir abl_comps_${gt} >/dev/null 2>&1; done
$PY - <<PYEOF
def load(gt):
    rows=open(f"abl_comps_{gt}/evaluation_best_f_w.tsv").read().splitlines();h=rows[0].split("\t")
    fn=h.index("filename");ns=h.index("ns");fw=h.index("f_w");NS={"molecular_function":"MF","biological_process":"BP","cellular_component":"CC"}
    o={}
    for r in rows[1:]:
        c=r.split("\t");o.setdefault(c[fn].replace(".tsv",""),{})[NS[c[ns]]]=float(c[fw])
    return o
d=load("gt_no");c=load("gt_no_cleanA")
mean=lambda x:sum(x.get(a,0) for a in("MF","BP","CC"))/3
print("%-13s %6s %6s"%("component","dirty","clean"))
for n in sorted(d,key=lambda n:-mean(c.get(n,{}))): print("%-13s %6.3f %6.3f"%(n,mean(d.get(n,{})),mean(c.get(n,{}))))
PYEOF
'
