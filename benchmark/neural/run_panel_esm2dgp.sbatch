#!/usr/bin/env bash
#SBATCH --job-name=gspa-panel-esm2dgp
#SBATCH --partition=debug
#SBATCH --nodelist=node005,node006
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 02:00:00
#SBATCH -o /data/hohndor/gspa-neural/logs/panel-esm2dgp-%A_%a.out
#SBATCH -e /data/hohndor/gspa-neural/logs/panel-esm2dgp-%A_%a.err

# Run ESM2 DeepGO-Plus sidecar on a single panel genome, selected by
# SLURM_ARRAY_TASK_ID (1-based row index into panel_manifest.tsv).
# Output: /data/hohndor/gspa-neural/panel/preds/<tag>/esm2-deepgoplus.tsv

set -euo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate gapfix
export PYTHONUNBUFFERED=1

PANEL=/data/hohndor/gspa-neural/panel
NEURAL=/data/hohndor/gspa-neural/benchmark/neural
WORK=/data/hohndor/gspa-neural/work
GAPFIX=/data/hohndor/gapfix
MANIFEST=$PANEL/panel_manifest.tsv

# Header is line 1; row N is on line N+1.
ROW=$(( SLURM_ARRAY_TASK_ID + 1 ))
TAG=$(awk -F'\t' -v r=$ROW 'NR==r{print $1}' $MANIFEST)
FASTA=$PANEL/fastas/${TAG}_protein.faa
OUT=$PANEL/preds/$TAG
mkdir -p $OUT

if [[ ! -s $FASTA ]]; then
    echo "FASTA missing: $FASTA"; exit 2
fi

N_PROT=$(grep -c '^>' $FASTA)
echo "=== $TAG ($N_PROT proteins) ==="

# Write a manifest with $OUT as output_dir so the sidecar writes there.
cat > $OUT/manifest.tsv <<EOF
tag	fasta_path	output_dir
$TAG	$FASTA	$OUT
EOF

export TORCH_HOME=$GAPFIX/data/esm2
python $NEURAL/run_neural_predictors.py \
    --predictor esm2-deepgoplus \
    --manifest $OUT/manifest.tsv \
    --checkpoint $GAPFIX/results/esm2-train-1950/best.pt \
    --terms $WORK/go_terms_5707.txt \
    --esm2-model esm2_t33_650M_UR50D \
    --batch-size 8 \
    --min-score 0.05 \
    --seed 42

echo "=== $TAG output TSV ==="
wc -l $OUT/${TAG}.esm2-deepgoplus.tsv
head -3 $OUT/${TAG}.esm2-deepgoplus.tsv
echo DONE $TAG; date
