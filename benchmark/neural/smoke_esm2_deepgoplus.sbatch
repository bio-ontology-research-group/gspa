#!/usr/bin/env bash
#SBATCH --job-name=gspa-esm2dgp-smoke
#SBATCH --partition=debug
#SBATCH --nodelist=node005
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 02:00:00
#SBATCH -o /data/hohndor/gspa-neural/logs/esm2dgp-smoke-%j.out
#SBATCH -e /data/hohndor/gspa-neural/logs/esm2dgp-smoke-%j.err

# Smoke test: run the ESM2 DeepGO-Plus sidecar on the first N rows of the
# DeepGO-Plus temporal holdout and emit predictions into
# /data/hohndor/gspa-neural/work/ .

set -euo pipefail

source /storage/miniforge3/etc/profile.d/conda.sh
conda activate gapfix

NEURAL=/data/hohndor/gspa-neural/benchmark/neural
WORK=/data/hohndor/gspa-neural/work
GAPFIX=/data/hohndor/gapfix

cd "$WORK"

# 1) Dump test FASTA + truth TSVs (fast)
python "$NEURAL/dump_test_fasta.py" \
    --in "$GAPFIX/data/deepgoplus-real/data/test_data.parquet" \
    --out-prefix "$WORK/dgp_test"

# 2) Take the first 200 proteins for the smoke test
head -400 "$WORK/dgp_test.faa" > "$WORK/dgp_test_200.faa"   # ~200 seqs
head -1001 "$WORK/dgp_test_truth_exp.tsv" > "$WORK/dgp_test_200_truth_exp.tsv"

# 3) Build a one-row manifest pointing at the 200-protein FASTA
cat > "$WORK/manifest_smoke.tsv" <<EOF
tag	fasta_path	output_dir
smoke	$WORK/dgp_test_200.faa	$WORK
EOF

# 4) Run the sidecar
export TORCH_HOME=/data/hohndor/gapfix/data/esm2
python "$NEURAL/run_neural_predictors.py" \
    --predictor esm2-deepgoplus \
    --manifest "$WORK/manifest_smoke.tsv" \
    --checkpoint "$GAPFIX/results/esm2-train-1950/best.pt" \
    --terms "$WORK/go_terms_5707.txt" \
    --esm2-model esm2_t33_650M_UR50D \
    --batch-size 8 \
    --min-score 0.05 \
    --seed 42

# 5) Report
echo "=== predictions TSV head ==="
head -5 "$WORK/smoke.esm2-deepgoplus.tsv" || true
echo "=== counts ==="
wc -l "$WORK/smoke.esm2-deepgoplus.tsv" "$WORK/dgp_test_200.faa" "$WORK/dgp_test_200_truth_exp.tsv"
echo DONE; date
