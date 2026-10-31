#!/bin/bash
#SBATCH --job-name=gspa-ledger-replay
#SBATCH --partition=debug
#SBATCH --nodelist=node002
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/data/hohndor/gspa-v3/price-production-v1/logs/replay-%j.log
set -euo pipefail
ROOT=/data/hohndor/gspa-v3/price-production-v1
export TMPDIR=/scratch/$SLURM_JOB_ID
mkdir -p "$TMPDIR"
export PYTHONPATH="$ROOT/code"
cd "$ROOT/code"
"$ROOT/venv/bin/python" -m unittest -v tests.test_streaming_ledger tests.test_price_production_collection tests.test_price_capacity
/usr/bin/time -v -o "$ROOT/logs/replay-$SLURM_JOB_ID.resources.txt" "$ROOT/venv/bin/python" "$ROOT/replay_ledger.py" "$ROOT/replay-$SLURM_JOB_ID"
