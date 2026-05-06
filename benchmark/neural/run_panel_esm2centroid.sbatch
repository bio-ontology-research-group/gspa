#!/usr/bin/env bash
#SBATCH --job-name=gspa-panel-esm2c
#SBATCH --partition=debug
#SBATCH --gres=gpu:1
#SBATCH --exclude=node001,node002,node003,node007
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -o /data/hohndor/gspa-neural/logs/panel-esm2c-%A_%a.out
#SBATCH -e /data/hohndor/gspa-neural/logs/panel-esm2c-%A_%a.err

set -euo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate gapfix
export PYTHONUNBUFFERED=1
export TORCH_HOME=/data/hohndor/gapfix/data/esm2

PANEL=/data/hohndor/gspa-neural/panel
NEURAL=/data/hohndor/gspa-neural/benchmark/neural
MANIFEST=$PANEL/panel_manifest.tsv
CENTROIDS=/data/hohndor/gspa-neural/models/esm2_centroids_t12.npz

ROW=$(( SLURM_ARRAY_TASK_ID + 1 ))
TAG=$(awk -F'\t' -v r=$ROW 'NR==r{print $1}' $MANIFEST)
FASTA=$PANEL/fastas/${TAG}_protein.faa
OUT=$PANEL/preds/$TAG
mkdir -p $OUT

[[ -s $FASTA ]]     || { echo "FASTA missing: $FASTA"; exit 2; }
[[ -s $CENTROIDS ]] || { echo "centroid NPZ missing: $CENTROIDS"; exit 2; }

cat > $OUT/manifest.tsv <<EOF
tag	fasta_path	output_dir
$TAG	$FASTA	$OUT
EOF

python $NEURAL/run_neural_predictors.py \
    --predictor esm2-centroid \
    --manifest $OUT/manifest.tsv \
    --centroid-db $CENTROIDS \
    --esm2-model esm2_t12_35M_UR50D \
    --batch-size 32 \
    --min-score 0.1 \
    --seed 42

wc -l $OUT/${TAG}.esm2-centroid.tsv
echo DONE $TAG; date
