#!/bin/bash
# Score a 3-column prediction TSV (protein\tterm\tscore) with the official
# cafaeval against all three CAFA6 knowledge classes and report IA-weighted
# f_w per class + the official 3-class mean.
# Usage: eval_classes.sh <preds.tsv> <label>
set -eo pipefail
cd "$(dirname "$0")"
PREDS="$1"; LABEL="${2:-preds}"
CAFA=~/Public/software/cafa6
OBO="$CAFA/go.obo"; IA="$CAFA/kaggle-official/IA.tsv"
PY=.venv/bin/python
WORK=$(mktemp -d)
mkdir -p "$WORK/preds"; cp "$PREDS" "$WORK/preds/p.tsv"

declare -A M
for cls in no limited partial; do
  cafaeval "$OBO" "$WORK/preds" "gt/gt_${cls}.tsv" -ia "$IA" -norm cafa -prop max \
    -th_step 0.01 -out_dir "$WORK/res_${cls}" >/dev/null 2>&1
  M[$cls]=$($PY - "$WORK/res_${cls}/evaluation_best_f_w.tsv" <<'PY'
import sys
rows=open(sys.argv[1]).read().splitlines(); h=rows[0].split('\t'); fi=h.index('f_w'); ni=h.index('ns')
m={r.split('\t')[ni]:float(r.split('\t')[fi]) for r in rows[1:]}
mf,bp,cc=m.get('molecular_function',0),m.get('biological_process',0),m.get('cellular_component',0)
print(f"{mf:.3f} {bp:.3f} {cc:.3f} {(mf+bp+cc)/3:.4f}")
PY
)
done
rm -rf "$WORK"
echo "=== $LABEL  (IA-weighted f_w: MF BP CC mean) ==="
echo "  no-knowledge : ${M[no]}"
echo "  limited      : ${M[limited]}"
echo "  partial      : ${M[partial]}"
$PY - "${M[no]}" "${M[limited]}" "${M[partial]}" <<'PY'
import sys
means=[float(a.split()[3]) for a in sys.argv[1:4]]
print(f"  OFFICIAL 3-class mean = {sum(means)/3:.4f}  (no {means[0]:.3f} | lim {means[1]:.3f} | part {means[2]:.3f})")
PY
