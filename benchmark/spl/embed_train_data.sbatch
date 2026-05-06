#!/usr/bin/env bash
#SBATCH --job-name=gspa-spl-embed
#SBATCH --partition=debug
#SBATCH --nodelist=node005
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH -o /data/hohndor/gspa-neural/logs/spl-embed-%j.out
#SBATCH -e /data/hohndor/gspa-neural/logs/spl-embed-%j.err

set -euo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate gapfix

cd /data/hohndor/gspa-neural/work
export TORCH_HOME=/home/hohndor/.cache/gapfix/torch-hub
export PYTHONUNBUFFERED=1

# Train pool (labels = prop_annotations)
python -u /data/hohndor/gspa-neural/benchmark/spl/embed_train_data.py \
    --parquet /data/hohndor/gapfix/data/deepgoplus-real/data/train_data.parquet \
    --terms /data/hohndor/gspa-neural/work/go_terms_5707.txt \
    --label-col prop_annotations \
    --model esm2_t33_650M_UR50D \
    --batch-size 8 \
    --out /data/hohndor/gspa-neural/work/pooled_train.npz

# Test pool (labels = exp_annotations for eval; prop_annotations for F-max scoring)
python -u /data/hohndor/gspa-neural/benchmark/spl/embed_train_data.py \
    --parquet /data/hohndor/gapfix/data/deepgoplus-real/data/test_data.parquet \
    --terms /data/hohndor/gspa-neural/work/go_terms_5707.txt \
    --label-col prop_annotations \
    --model esm2_t33_650M_UR50D \
    --batch-size 8 \
    --out /data/hohndor/gspa-neural/work/pooled_test.npz

echo DONE; date
ls -la /data/hohndor/gspa-neural/work/pooled_*.npz
