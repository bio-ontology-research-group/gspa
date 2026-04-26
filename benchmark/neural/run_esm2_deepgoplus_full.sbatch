#!/usr/bin/env bash
#SBATCH --job-name=gspa-esm2dgp-full
#SBATCH --partition=debug
#SBATCH --nodelist=node005
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH -o /data/hohndor/gspa-neural/logs/esm2dgp-full-%j.out
#SBATCH -e /data/hohndor/gspa-neural/logs/esm2dgp-full-%j.err

# Run the ESM2 DeepGO-Plus sidecar on the full DeepGO-Plus temporal holdout
# (4180 proteins), transform to integrated-TSV shape, and report F-max
# against experimental-evidence ground truth.

set -euo pipefail

source /storage/miniforge3/etc/profile.d/conda.sh
conda activate gapfix

NEURAL=/data/hohndor/gspa-neural/benchmark/neural
WORK=/data/hohndor/gspa-neural/work
GAPFIX=/data/hohndor/gapfix

cd "$WORK"

# Use the already-cached ESM2 weights in the user's home cache (skip the
# 2.5 GB redownload the first smoke triggered by setting TORCH_HOME).
export TORCH_HOME=/home/hohndor/.cache/gapfix/torch-hub

# 1) Manifest pointing at the full FASTA
cat > "$WORK/manifest_full.tsv" <<EOF
tag	fasta_path	output_dir
full	$WORK/dgp_test.faa	$WORK
EOF

# 2) Sidecar inference
python "$NEURAL/run_neural_predictors.py" \
    --predictor esm2-deepgoplus \
    --manifest "$WORK/manifest_full.tsv" \
    --checkpoint "$GAPFIX/results/esm2-train-1950/best.pt" \
    --terms "$WORK/go_terms_5707.txt" \
    --esm2-model esm2_t33_650M_UR50D \
    --batch-size 8 \
    --min-score 0.01 \
    --seed 42

# 3) Transform into integrated-TSV shape (with GO aspect)
python "$NEURAL/sidecar_to_integrated.py" \
    --in "$WORK/full.esm2-deepgoplus.tsv" \
    --out "$WORK/full.integrated.tsv" \
    --go-obo "$GAPFIX/data/deepgoplus-real/data/go.obo"

# 4) F-max via benchmark_pgap_v2 (exp-only truth)
python /data/hohndor/gspa/benchmark/benchmark_pgap_v2.py \
    --gspa "$WORK/full.integrated.tsv" \
    --truth "exp:$WORK/dgp_test_truth_exp.tsv" \
    --tag esm2-deepgoplus \
    --annotation-type GO \
    --n-bootstrap 100 \
    > "$WORK/full_fmax.json"

echo "=== F-max summary ==="
python -c "
import json
with open('$WORK/full_fmax.json') as f: r = json.load(f)
for name, block in r['by_truth'].items():
    print(f'truth={name}: {block[\"truth_proteins\"]} proteins, {block[\"truth_annotations\"]} annotations')
    for res in block['results']:
        print(f'  method={res[\"method\"]}')
        print(f'    overall micro F-max = {res[\"fmax_overall\"]:.4f} CI {res[\"fmax_ci\"]}')
        print(f'    overall CAFA F-max  = {res[\"fmax_cafa_overall\"]:.4f} CI {res[\"fmax_cafa_ci\"]}')
        for asp in ('MF','BP','CC'):
            k = f'fmax_{asp}'
            if k in res:
                print(f'    {asp} micro F-max    = {res[k]:.4f} CI {res[k+\"_ci\"]}')
"

echo DONE; date
