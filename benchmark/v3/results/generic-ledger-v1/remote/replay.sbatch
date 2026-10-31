#!/bin/bash
#SBATCH --job-name=gspa-generic-ledger
#SBATCH --partition=debug
#SBATCH --nodelist=node004
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/data/hohndor/gspa-v3/generic-ledger-v2/logs/%j.log
set -euo pipefail
ROOT=/data/hohndor/gspa-v3/generic-ledger-v2
export TMPDIR=/scratch/$SLURM_JOB_ID
mkdir -p "$TMPDIR"
unset GSPA_V3_ENGINE_SMOKE
mkdir "$ROOT/source"
tar xzf "$ROOT/source.tar.gz" -C "$ROOT/source"
export PYTHONPATH="$ROOT/source/benchmark/v3:$ROOT/source/benchmark/v3/tests"
PYTHON=/data/hohndor/gspa-v3/price-production-v1/venv/bin/python
cd "$ROOT/source"
"$PYTHON" -m unittest -v test_streaming_ledger test_runner test_hardening test_verified_batch test_price_production_collection test_price_capacity > "$ROOT/logs/$SLURM_JOB_ID.tests.log" 2>&1
/usr/bin/time -v -o "$ROOT/logs/$SLURM_JOB_ID.resources.txt" "$PYTHON" "$ROOT/replay.py" "$SLURM_JOB_ID"
