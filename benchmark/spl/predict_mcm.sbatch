#!/usr/bin/env bash
#SBATCH --job-name=gspa-mcm-pred
#SBATCH --partition=debug
#SBATCH --nodelist=node005
#SBATCH --gres=gpu:1
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH -o /data/hohndor/gspa-neural/logs/mcm-pred-%j.out
#SBATCH -e /data/hohndor/gspa-neural/logs/mcm-pred-%j.err

set -euo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate gapfix
export PYTHONUNBUFFERED=1

W=/data/hohndor/gspa-neural/work

python -u /data/hohndor/gspa-neural/benchmark/spl/predict_mcm.py \
    --checkpoint $W/mcm-run/best.pt \
    --pooled $W/pooled_test.npz \
    --terms $W/go_terms_5707.txt \
    --out-tsv $W/mcm.predictions.tsv \
    --batch-size 32 \
    --min-score 0.01

echo DONE; date
