#!/bin/bash
#SBATCH --job-name=gspa-price-local-medium
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --exclude=node001,node005,node006,node007,node008,unimatrix
#SBATCH --cpus-per-task=4
#SBATCH --mem=112G
#SBATCH --time=10:00:00
#SBATCH --output=/data/hohndor/gspa-v3/price-production-v3/logs/medium-%A_%a.log
set -euo pipefail
ROOT=/data/hohndor/gspa-v3/price-production-v3
export TMPDIR=/scratch/$SLURM_JOB_ID
mkdir -p "$TMPDIR"
unset JAVA_TOOL_OPTIONS JDK_JAVA_OPTIONS _JAVA_OPTIONS
export PYTHONPATH="$ROOT/code"
INDEX=${SLURM_ARRAY_TASK_ID:-${1:?cell index required}}
cd "$ROOT"
/data/hohndor/gspa-v3/price-production-v1/venv/bin/python "$ROOT/code/run_price_local.py" --out "$ROOT/collection" --index "$INDEX" --attempt "job-$SLURM_JOB_ID" --heap 96g --timeout 14400 --scratch "$TMPDIR" --container /data/hohndor/gspa-v3/observed-gap-context-modes-v1-relocation/container-root
