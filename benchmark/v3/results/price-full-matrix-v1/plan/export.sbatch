#!/bin/bash
#SBATCH --job-name=gspa-price-full-export
#SBATCH --partition=debug
#SBATCH --nodelist=node003
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --array=0-131%2
#SBATCH --output=/data/hohndor/gspa-v3/price-full-matrix-v2/logs/export-%A_%a.log
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export TMPDIR="${SLURM_TMPDIR:-/scratch/$SLURM_JOB_ID}"
mkdir -p "$TMPDIR"
/usr/bin/time -v /data/hohndor/gspa-v3/tools/metabolic-venv/bin/python /data/hohndor/gspa-v3/price-full-matrix-v2/code/export_price_matrix_job.py --plan /data/hohndor/gspa-v3/price-full-matrix-v2/plan.json --index "$SLURM_ARRAY_TASK_ID" --scratch "$TMPDIR"
