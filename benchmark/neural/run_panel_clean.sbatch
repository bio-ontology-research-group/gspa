#!/usr/bin/env bash
#SBATCH --job-name=gspa-panel-clean
#SBATCH --partition=debug
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 08:00:00
#SBATCH -o /data/hohndor/gspa-neural/logs/panel-clean-%A_%a.out
#SBATCH -e /data/hohndor/gspa-neural/logs/panel-clean-%A_%a.err

set -euo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate gapfix
export PYTHONUNBUFFERED=1

PANEL=/data/hohndor/gspa-neural/panel
NEURAL=/data/hohndor/gspa-neural/benchmark/neural
MANIFEST=$PANEL/panel_manifest.tsv
MODEL=/data/hohndor/gspa-neural/models/CLEAN/app

ROW=$(( SLURM_ARRAY_TASK_ID + 1 ))
TAG=$(awk -F'\t' -v r=$ROW 'NR==r{print $1}' $MANIFEST)
FASTA=$PANEL/fastas/${TAG}_protein.faa
OUT=$PANEL/preds/$TAG
mkdir -p $OUT

[[ -s $FASTA ]] || { echo "FASTA missing: $FASTA"; exit 2; }
echo "=== $TAG ==="

cat > $OUT/manifest.tsv <<EOF
tag	fasta_path	output_dir
$TAG	$FASTA	$OUT
EOF

python $NEURAL/run_neural_predictors.py \
    --predictor clean \
    --manifest $OUT/manifest.tsv \
    --model-dir $MODEL \
    --min-score 0.1

wc -l $OUT/${TAG}.clean.tsv
head -3 $OUT/${TAG}.clean.tsv
echo DONE $TAG; date
