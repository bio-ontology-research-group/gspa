#!/usr/bin/env bash
#SBATCH --job-name=gspa-panel-ensemble
#SBATCH --partition=debug
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH -o /data/hohndor/gspa-neural/logs/panel-ensemble-%j.out
#SBATCH -e /data/hohndor/gspa-neural/logs/panel-ensemble-%j.err

set -euo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate gapfix
export PYTHONUNBUFFERED=1

PANEL=/data/hohndor/gspa-neural/panel
NEURAL=/data/hohndor/gspa-neural/benchmark/neural

# For each genome and each mode, fuse all available per-predictor TSVs.
MODES=(max mean rank)

while IFS=$'\t' read tag rest; do
    [[ "$tag" == "tag" ]] && continue
    preds_dir=$PANEL/preds/$tag
    [[ -d $preds_dir ]] || continue
    PRED_ARGS=()
    for p in esm2-deepgoplus proteinfer clean esm2-centroid foldseek; do
        f=$preds_dir/${tag}.${p}.tsv
        [[ -s $f ]] && PRED_ARGS+=(--pred "$f")
    done
    [[ ${#PRED_ARGS[@]} -lt 2 ]] && { echo "skip $tag (<2 predictors)"; continue; }
    for mode in "${MODES[@]}"; do
        out=$preds_dir/${tag}.ensemble-${mode}.tsv
        python $NEURAL/build_ensemble_preds.py \
            "${PRED_ARGS[@]}" \
            --out $out \
            --mode $mode \
            --min-score 0.0
        echo "  $tag/$mode: $(wc -l < $out) rows"
    done
done < $PANEL/panel_manifest.tsv

echo DONE; date
