#!/bin/bash
#SBATCH --job-name=gspa-dpz-recheck
#SBATCH --partition=debug
#SBATCH --nodelist=node003
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:30:00
#SBATCH --output=/data/hohndor/gspa-v3/price-capacity-v2/hash-recheck-%j.log
set -euo pipefail
export TMPDIR="${SLURM_TMPDIR:-/scratch/$SLURM_JOB_ID}"
python3 /data/hohndor/gspa-v3/price-capacity-v2/hash-recheck.py
