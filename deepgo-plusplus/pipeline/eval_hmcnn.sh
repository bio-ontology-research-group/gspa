#!/bin/bash
# Pull the ORIX hierarchy-aware head preds (bce/mcm/softreg × prostt5/esm2_3b),
# score each standalone with cafaeval on the leak-free clean GTs, print a per-aspect
# A/B table. Run locally after job train_head_hmcnn.sbatch completes.
set -uo pipefail
SCRATCH=/tmp/claude-1000/-home-leechuck-Public-software-gspa/45a645bf-2c55-4082-9718-d452c4bced25/scratchpad/hmcnn_preds
mkdir -p "$SCRATCH"

echo "[1] filter to clean-GT proteins on ORIX, then pull (avoids multi-GB transfer)"
ssh orix 'cd /mnt/data/u/hohndor/cafa6-plm/preds/hmcnn
mkdir -p flt
for f in head_prostt5_bce head_prostt5_mcm head_prostt5_softreg head_esm2_3b_bce head_esm2_3b_mcm head_esm2_3b_softreg; do
  [ -s $f.tsv ] && awk -F"\t" "NR==FNR{k[\$1];next} \$1 in k" eval_ids.txt $f.tsv > flt/$f.tsv && echo "  $f -> $(wc -l < flt/$f.tsv) rows"
done'
scp -q 'orix:/mnt/data/u/hohndor/cafa6-plm/preds/hmcnn/flt/head_*.tsv' "$SCRATCH/" || { echo "pull failed"; exit 1; }
ls -l "$SCRATCH"

echo "[2] push to ws P_hmcnn/"
ssh ws 'rm -rf ~/cafa6_cleaneval/P_hmcnn; mkdir -p ~/cafa6_cleaneval/P_hmcnn'
for f in "$SCRATCH"/head_*.tsv; do
  scp -q "$f" "ws:~/cafa6_cleaneval/P_hmcnn/$(basename "$f")"
done

echo "[3] cafaeval each file on the clean GTs"
ssh ws 'cd ~/cafa6_cleaneval; CE=~/cnn_work/.venv/bin/cafaeval; PY=~/cnn_work/.venv/bin/python
for gt in gt_no_cleanA gt_no_cleanB; do
  $CE go.obo P_hmcnn gt/${gt}.tsv -ia IA.tsv -norm cafa -prop max -th_step 0.01 -out_dir res_hmcnn_${gt} >/dev/null 2>&1 || true
done
$PY - <<PYEOF
import os
NS={"molecular_function":"MF","biological_process":"BP","cellular_component":"CC"}
def load(d):
    p=f"{d}/evaluation_best_f_w.tsv"
    if not os.path.exists(p): return {}
    rows=open(p).read().splitlines(); h=rows[0].split("\t")
    fn=h.index("filename"); ns=h.index("ns"); fw=h.index("f_w"); o={}
    for r in rows[1:]:
        c=r.split("\t"); o.setdefault(c[fn].replace(".tsv",""),{})[NS[c[ns]]]=float(c[fw])
    return o
mean=lambda x:sum(x.get(a,0) for a in("MF","BP","CC"))/3
for gt in ("gt_no_cleanA","gt_no_cleanB"):
    d=load(f"res_hmcnn_{gt}")
    if not d: continue
    print(f"\n===== {gt} (IA-weighted f_w, standalone) =====")
    print("%-26s %6s  %5s %5s %5s"%("head","mean","MF","BP","CC"))
    for c in sorted(d,key=lambda c:-mean(d[c])):
        v=d[c]; print("%-26s %6.3f  %5.3f %5.3f %5.3f"%(c,mean(v),v.get("MF",0),v.get("BP",0),v.get("CC",0)))
PYEOF'
