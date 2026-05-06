#!/usr/bin/env bash
#SBATCH --job-name=gspa-mcm-train
#SBATCH --partition=debug
#SBATCH --nodelist=node005
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 06:00:00
#SBATCH -o /data/hohndor/gspa-neural/logs/mcm-train-%j.out
#SBATCH -e /data/hohndor/gspa-neural/logs/mcm-train-%j.err

set -euo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate gapfix
export PYTHONUNBUFFERED=1

cd /data/hohndor/gspa-neural/work

python -u /data/hohndor/gspa-neural/benchmark/spl/train_mcm.py \
    --train-pooled /data/hohndor/gspa-neural/work/pooled_train.npz \
    --hierarchy /data/hohndor/gspa-neural/work/go_hierarchy_5707.npz \
    --out-dir /data/hohndor/gspa-neural/work/mcm-run \
    --hidden 2048 \
    --dropout 0.2 \
    --batch-size 32 \
    --lr 1e-3 \
    --wd 1e-5 \
    --epochs 100 \
    --patience 10 \
    --seed 0

echo DONE; date
