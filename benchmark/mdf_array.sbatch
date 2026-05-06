#!/bin/bash
#SBATCH --job-name=mdf_panel
#SBATCH --partition=debug
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=1-13%14
#SBATCH --output=/data/hohndor/mdf-runs/array/%A_%a.out

set -euo pipefail
mkdir -p /data/hohndor/mdf-runs/array
echo "host=$(hostname) idx=${SLURM_ARRAY_TASK_ID} jobid=${SLURM_JOB_ID:-?}"
bash /data/hohndor/gspa-v1.5.0/benchmark/mdf_array.sh
