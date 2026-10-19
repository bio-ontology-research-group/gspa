#!/bin/bash
#SBATCH --partition=freecycle-h200
#SBATCH --qos=freecycle
#SBATCH --requeue
#SBATCH --job-name=gspa-wg-gpu
#SBATCH --output=/mnt/data/u/hohndor/gspa-wg-cafa/logs/%x-%A_%a.out
#SBATCH --error=/mnt/data/u/hohndor/gspa-wg-cafa/logs/%x-%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --gres=gpu:nvidia_h200:1
#SBATCH --array=0-4

set -euo pipefail
mkdir -p /mnt/data/u/hohndor/gspa-wg-cafa/logs

TAGS=(wg_time_rollover wg_pid100 wg_pid70 wg_pid50 wg_pid30)
export GSPA_CACHE_TAG="${TAGS[$SLURM_ARRAY_TASK_ID]}"
export GSPA_ROOT=/mnt/data/u/hohndor/gspa-wg-cafa
export GSPA_CAFA_DIR=/mnt/data/u/hohndor/gspa-wg-cafa/cafa6
export PYTHONUNBUFFERED=1

echo "[$(date -Is)] node=$(hostname) jobid=${SLURM_JOB_ID} task=${SLURM_ARRAY_TASK_ID} tag=${GSPA_CACHE_TAG}"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv

cd /home/hohndor/pixi-envs/gspa-glm
CONDA_OVERRIDE_CUDA=12.6 /home/hohndor/.pixi/bin/pixi run --frozen python \
    /mnt/data/u/hohndor/gspa-wg-cafa/scripts/run_cafa_whole_genome_splits.py --gpu-only
