#!/bin/bash
#SBATCH --job-name=gspa-index-dpz-check
#SBATCH --partition=debug
#SBATCH --nodelist=node003
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=01:00:00
#SBATCH --output=/data/hohndor/gspa-v3/price-index-comparison-dpz-v1/recheck-%j.log
set -euo pipefail
export TMPDIR="${SLURM_TMPDIR:-/scratch/$SLURM_JOB_ID}"
mkdir -p "$TMPDIR"
cd /data/hohndor/gspa-v3/price-index-comparison-dpz-v1
/data/hohndor/gspa-v3/tools/metabolic-venv/bin/python recheck-dpz.py
