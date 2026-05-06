#!/usr/bin/env bash
#SBATCH --job-name=gspa-spl-cc
#SBATCH --partition=debug
#SBATCH --nodelist=node004
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 06:00:00
#SBATCH -o /data/hohndor/gspa-neural/logs/spl-cc-train-%j.out
#SBATCH -e /data/hohndor/gspa-neural/logs/spl-cc-train-%j.err

set -euo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate gapfix
export PYTHONUNBUFFERED=1

cd /data/hohndor/gspa-neural/work

python -u /data/hohndor/gspa-neural/benchmark/spl/train_spl.py \
    --train-pooled /data/hohndor/gspa-neural/work/pooled_train.npz \
    --hierarchy /data/hohndor/gspa-neural/work/go_hierarchy_CC.npz \
    --sdd /data/hohndor/gspa-neural/work/constraints/go_true_path_CC.sdd \
    --vtree /data/hohndor/gspa-neural/work/constraints/go_true_path_CC.vtree \
    --hmc-utils /data/hohndor/gspa-neural/benchmark/spl/hmc-utils \
    --out-dir /data/hohndor/gspa-neural/work/spl-cc-run \
    --gate-layers 1280,512,512 \
    --num-reps 1 \
    --lr 1e-3 \
    --wd 1e-5 \
    --batch-size 128 \
    --epochs 50 \
    --patience 8 \
    --seed 0

echo DONE; date
