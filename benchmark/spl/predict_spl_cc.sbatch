#!/usr/bin/env bash
#SBATCH --job-name=gspa-spl-cc-pred
#SBATCH --partition=debug
#SBATCH --nodelist=node005
#SBATCH --gres=gpu:1
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH -o /data/hohndor/gspa-neural/logs/spl-cc-pred-%j.out
#SBATCH -e /data/hohndor/gspa-neural/logs/spl-cc-pred-%j.err

set -euo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate gapfix
export PYTHONUNBUFFERED=1

W=/data/hohndor/gspa-neural/work

python -u /data/hohndor/gspa-neural/benchmark/spl/predict_spl.py \
    --checkpoint $W/spl-cc-run/best.pt \
    --sdd $W/constraints/go_true_path_CC.sdd \
    --vtree $W/constraints/go_true_path_CC.vtree \
    --hierarchy $W/go_hierarchy_CC.npz \
    --pooled $W/pooled_test.npz \
    --hmc-utils /data/hohndor/gspa-neural/benchmark/spl/hmc-utils \
    --out-tsv $W/spl_cc.predictions.tsv \
    --batch-size 64

echo DONE; date
