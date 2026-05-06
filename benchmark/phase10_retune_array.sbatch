#!/bin/bash
#SBATCH --job-name=p10arr
#SBATCH --partition=debug
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --array=1-40%14
#SBATCH --output=/data/hohndor/gspa/proteomes/bench10/phase10_retune/array/%A_%a.out

set -euo pipefail
mkdir -p /data/hohndor/gspa/proteomes/bench10/phase10_retune/array
echo "host=$(hostname) idx=${SLURM_ARRAY_TASK_ID} jobid=${SLURM_JOB_ID:-?}"
bash /data/hohndor/gspa-v1.5.0/benchmark/phase10_retune_array.sh
