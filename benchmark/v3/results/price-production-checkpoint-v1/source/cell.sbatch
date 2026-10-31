#!/bin/bash
#SBATCH --job-name=gspa-price-cell
#SBATCH --partition=debug
#SBATCH --cpus-per-task=2
#SBATCH --mem=36G
#SBATCH --time=06:00:00
#SBATCH --output=/data/hohndor/gspa-v3/price-production-v1/logs/cell-%A_%a.log
set -euo pipefail
ROOT=/data/hohndor/gspa-v3/price-production-v1
export TMPDIR=/scratch/$SLURM_JOB_ID
mkdir -p "$TMPDIR"
unset JAVA_TOOL_OPTIONS JDK_JAVA_OPTIONS _JAVA_OPTIONS
export PYTHONPATH="$ROOT/code"
cd "$ROOT"
INDEX=${SLURM_ARRAY_TASK_ID:-0}
"$ROOT/venv/bin/python" "$ROOT/code/run_price_production_collection.py" run --out "$ROOT/collection" --index "$INDEX" --attempt "job-$SLURM_JOB_ID" --heap 24g --timeout 7200 --scratch "$TMPDIR"
