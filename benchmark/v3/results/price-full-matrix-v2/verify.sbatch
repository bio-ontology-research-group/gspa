#!/bin/bash
#SBATCH --job-name=gspa-full-export-audit
#SBATCH --partition=debug
#SBATCH --nodelist=node003
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=06:00:00
#SBATCH --output=/data/hohndor/gspa-v3/price-full-collection-audit-v1/audit-%j.log
set -euo pipefail
export TMPDIR="${SLURM_TMPDIR:-/scratch/$SLURM_JOB_ID}"
mkdir -p "$TMPDIR"
python3 /data/hohndor/gspa-v3/price-full-collection-audit-v1/verify_price_export_collection.py --plan /data/hohndor/gspa-v3/price-full-matrix-v2/plan.json --out /data/hohndor/gspa-v3/price-full-collection-audit-v1/verification.json
