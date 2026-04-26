#!/usr/bin/env bash
#SBATCH --job-name=gspa-eval-one
#SBATCH --partition=debug
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 00:30:00
#SBATCH -o /data/hohndor/gspa-neural/logs/eval-%A_%a.out
#SBATCH -e /data/hohndor/gspa-neural/logs/eval-%A_%a.err

set -euo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate gapfix
export PYTHONUNBUFFERED=1

PANEL=/data/hohndor/gspa-neural/panel
EVAL=/data/hohndor/gspa-neural/benchmark
ASPECT_MAP=/data/hohndor/gspa/reference/go_aspect_map.tsv
MANIFEST=$PANEL/panel_manifest.tsv
OUT=$PANEL/results
mkdir -p $OUT

ROW=$(( SLURM_ARRAY_TASK_ID + 1 ))
TAG=$(awk -F'\t' -v r=$ROW 'NR==r{print $1}' $MANIFEST)
[[ -z "$TAG" ]] && { echo "no tag for row $ROW"; exit 0; }

GO_PREDICTORS=(esm2-deepgoplus proteinfer esm2-centroid foldseek ensemble-max ensemble-mean ensemble-rank)
GO_TRUTHS=(truth_all_refseq_prop truth_all_refseq truth_exp_refseq truth_sprot_refseq_prop truth_sprot_refseq)
EC_PREDICTORS=(clean proteinfer ensemble-max ensemble-mean ensemble-rank)
EC_TRUTHS=(ec_refseq ec_sprot_refseq)

OUT_JSONL=$OUT/${TAG}.jsonl
: > $OUT_JSONL

eval_pair() {
    local pred=$1 truth_name=$2 ann=$3
    local pred_tsv=$PANEL/preds/$TAG/${TAG}.${pred}.tsv
    local truth_tsv=$PANEL/truth/${TAG}_${truth_name}.tsv
    [[ -s $pred_tsv && -s $truth_tsv ]] || return 0
    echo "  $TAG $pred $truth_name ($ann)"
    python $EVAL/evaluate_panel.py \
        --predictor-tsv $pred_tsv \
        --truth $truth_tsv \
        --annotation-type $ann \
        --go-aspect-map $ASPECT_MAP \
        --tag $TAG \
        --predictor $pred \
        --truth-name $truth_name \
        --n-bootstrap 0 \
        >> $OUT_JSONL 2>> $OUT/${TAG}.err \
        || echo "  FAILED $TAG/$pred/$truth_name"
}

for pred in "${GO_PREDICTORS[@]}"; do
    for truth_name in "${GO_TRUTHS[@]}"; do
        eval_pair "$pred" "$truth_name" GO
    done
done
for pred in "${EC_PREDICTORS[@]}"; do
    for truth_name in "${EC_TRUTHS[@]}"; do
        eval_pair "$pred" "$truth_name" EC
    done
done

echo DONE $TAG; date
wc -l $OUT_JSONL
