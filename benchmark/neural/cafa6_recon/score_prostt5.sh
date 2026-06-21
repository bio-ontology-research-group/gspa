#!/bin/bash
# Run the leak-free logreg LTR integrator for a given component list and score
# the OOF predictions with the official cafaeval on the no-knowledge GT.
# Usage: score_prostt5.sh <out_dir> <component_list>
set -eo pipefail
cd "$(dirname "$0")"

OUT="$1"; COMPS="$2"
CAFA=~/Public/software/cafa6
DAG="$CAFA/go-dag.tsv"
OBO="$CAFA/go.obo"
IA="$CAFA/kaggle-official/IA.tsv"
TRAIN_TERMS="$CAFA/kaggle-official/train_terms.tsv"
TAXON="$CAFA/kaggle-official/testsuperset-taxon-list.tsv"
GT=gt/gt_no.tsv
PY=.venv/bin/python

echo "=== integrator: $COMPS -> $OUT ==="
$PY ../train_ltr_integrator.py \
  --components components --gt "$GT" \
  --dag "$DAG" --ia "$IA" --train-terms "$TRAIN_TERMS" --taxon "$TAXON" \
  --model logreg --features scores \
  --component-list "$COMPS" --out "$OUT"

echo "=== cafaeval ==="
cafaeval "$OBO" "$OUT/preds" "$GT" -ia "$IA" -norm cafa -prop max -th_step 0.01 -out_dir "$OUT/results"

echo "=== f_w by aspect (mean over MF/BP/CC) ==="
$PY - "$OUT/results/evaluation_best_f_w.tsv" <<'PY'
import sys
rows = open(sys.argv[1]).read().splitlines()
hdr = rows[0].split('\t'); fi = hdr.index('f_w'); ni = hdr.index('ns')
m = {}
for r in rows[1:]:
    c = r.split('\t'); m[c[ni]] = float(c[fi])
mf = m.get('molecular_function', 0.0); bp = m.get('biological_process', 0.0); cc = m.get('cellular_component', 0.0)
print(f"MF={mf:.3f} BP={bp:.3f} CC={cc:.3f} mean={(mf+bp+cc)/3:.3f}")
PY
