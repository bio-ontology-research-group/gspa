#!/bin/bash
#SBATCH --job-name=gspa-clean-capacity
#SBATCH --partition=debug
#SBATCH --nodelist=node005
#SBATCH --cpus-per-task=4
#SBATCH --mem=192G
#SBATCH --time=12:00:00
#SBATCH --output=/data/hohndor/gspa-v3/price-production-v1/logs/dense-%j.log
set -euo pipefail
ROOT=/data/hohndor/gspa-v3/price-production-v1
export TMPDIR=/scratch/$SLURM_JOB_ID
mkdir -p "$TMPDIR"
unset JAVA_TOOL_OPTIONS JDK_JAVA_OPTIONS _JAVA_OPTIONS
export PYTHONPATH="$ROOT/code"
cd "$ROOT"
/usr/bin/time -v -o "$ROOT/logs/dense-$SLURM_JOB_ID.resources.txt" "$ROOT/venv/bin/python" "$ROOT/dense_profile.py"
