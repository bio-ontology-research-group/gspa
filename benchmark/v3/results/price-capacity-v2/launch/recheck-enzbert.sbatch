#!/bin/bash
#SBATCH --job-name=gspa-price-recheck
#SBATCH --partition=debug
#SBATCH --nodelist=node003
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=/data/hohndor/gspa-v3/price-capacity-v2/recheck-%j.log
set -euo pipefail
export TMPDIR="${SLURM_TMPDIR:-/scratch/$SLURM_JOB_ID}"
mkdir -p "$TMPDIR"
/data/hohndor/gspa-v3/tools/metabolic-venv/bin/python /data/hohndor/gspa-v3/price-capacity-v2/recheck-enzbert.py
