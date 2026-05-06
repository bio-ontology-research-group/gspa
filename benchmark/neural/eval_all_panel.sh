#!/usr/bin/env bash
#SBATCH --job-name=gspa-panel-eval
#SBATCH --partition=debug
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH -o /data/hohndor/gspa-neural/logs/panel-eval-%j.out
#SBATCH -e /data/hohndor/gspa-neural/logs/panel-eval-%j.err

# For each (genome, predictor) pair, compute F-max / CAFA F-max / Smin
# vs SwissProt truth. Emit one JSON row per pair; concatenate to a master
# TSV for the results table.

set -euo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate gapfix

PANEL=/data/hohndor/gspa-neural/panel
EVAL=/data/hohndor/gspa-neural/benchmark
ASPECT_MAP=/data/hohndor/gspa/reference/go_aspect_map.tsv

OUT=$PANEL/results
mkdir -p $OUT

PREDICTORS=(esm2-deepgoplus proteinfer clean esm2-centroid foldseek)
TRUTHS=(truth_all_refseq_prop truth_all_refseq truth_exp_refseq)

: > $OUT/eval_all.jsonl

while IFS=$'\t' read tag rest; do
    [[ "$tag" == "tag" ]] && continue
    for pred in "${PREDICTORS[@]}"; do
        pred_tsv=$PANEL/preds/$tag/${tag}.${pred}.tsv
        [[ -s $pred_tsv ]] || continue
        for truth_name in "${TRUTHS[@]}"; do
            truth_tsv=$PANEL/truth/${tag}_${truth_name}.tsv
            [[ -s $truth_tsv ]] || continue
            ann=GO
            [[ "$pred" == "clean" ]] && ann=EC
            echo "  $tag $pred $truth_name ($ann)"
            python $EVAL/evaluate_panel.py \
                --predictor-tsv $pred_tsv \
                --truth $truth_tsv \
                --annotation-type $ann \
                --go-aspect-map $ASPECT_MAP \
                --tag $tag \
                --predictor $pred \
                --truth-name $truth_name \
                --n-bootstrap 0 \
                >> $OUT/eval_all.jsonl 2>>$OUT/eval_all.err \
                || echo "  FAILED $tag/$pred/$truth_name"
        done
    done
done < $PANEL/panel_manifest.tsv

# Flatten JSONL to TSV
python - <<'EOF' > $OUT/eval_all.tsv
import json
from pathlib import Path
rows = []
for line in Path('/data/hohndor/gspa-neural/panel/results/eval_all.jsonl').open():
    line = line.strip()
    if not line:
        continue
    try:
        rows.append(json.loads(line))
    except Exception:
        # multi-line JSON objects (indent=2): accumulate until valid
        pass

# Since evaluate_panel.py prints indent=2, we need to parse the full file
# as a sequence of concatenated JSON objects. Use a decoder loop.
txt = Path('/data/hohndor/gspa-neural/panel/results/eval_all.jsonl').read_text()
rows = []
dec = json.JSONDecoder()
i = 0
while i < len(txt):
    while i < len(txt) and txt[i] in ' \n\t\r':
        i += 1
    if i >= len(txt):
        break
    obj, end = dec.raw_decode(txt[i:])
    rows.append(obj)
    i += end

# Collect columns
cols = []
seen = set()
for r in rows:
    for k in r.keys():
        if k not in seen:
            cols.append(k)
            seen.add(k)

# Write
import sys
out = sys.stdout
out.write('\t'.join(cols) + '\n')
for r in rows:
    out.write('\t'.join(str(r.get(c, '')) for c in cols) + '\n')
EOF

echo DONE; date
wc -l $OUT/eval_all.tsv $OUT/eval_all.jsonl
