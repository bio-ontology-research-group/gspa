#!/bin/bash
#SBATCH --job-name=gspa_sota_wg
#SBATCH --partition=freecycle-h200,freecycle-h100
#SBATCH --qos=freecycle
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --array=0-4
#SBATCH --output=/mnt/data/u/hohndor/gspa_wg_cafa/logs/sota_wg_cafa_%A_%a.out
#SBATCH --error=/mnt/data/u/hohndor/gspa_wg_cafa/logs/sota_wg_cafa_%A_%a.err

set -euo pipefail

export PATH="/home/hohndor/pixi-envs/gspa-glm/.pixi/envs/default/bin:${PATH}"
export OMP_NUM_THREADS=8

cd /mnt/data/u/hohndor/gspa_wg_cafa
echo "=== Starting External SOTA WG CAFA split ${SLURM_ARRAY_TASK_ID} on $(hostname) at $(date) ==="
nvidia-smi

python3 scripts/train_eval_external_sota_wg_cafa.py --split-idx "${SLURM_ARRAY_TASK_ID}"

echo "=== Completed External SOTA WG CAFA split ${SLURM_ARRAY_TASK_ID} at $(date) ==="
