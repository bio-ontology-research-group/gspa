#!/usr/bin/env bash
#SBATCH --job-name=eq_gapsmith_batch
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --time=24:00:00
#SBATCH --exclusive
#SBATCH --exclude=nid00022
#SBATCH --no-requeue

# Stage 2: run gapsmith 0.1.0 across a staged batch of ~500 genomes on an exclusive
# 192-core Shaheen III workq node (48 concurrent genomes x 4 threads/genome).
# Resumable: any genome already DONE in done_status.tsv is skipped.
# Outputs per genome (<mag>.xml, <mag>-draft.xml, tables, logs, /usr/bin/time -v
# resources, input_sha256.txt) are packed per batch into <BATCH_DIR>/outputs.tar.

set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
if [ -f "$SCRIPT_DIR/common.sh" ]; then
    # shellcheck source=/dev/null
    source "$SCRIPT_DIR/common.sh"
else
    # shellcheck source=/dev/null
    source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"
fi

BATCH_IDX="${SLURM_ARRAY_TASK_ID:-${1:-0}}"
BATCH_ID="$(printf 'batch_%03d' "$BATCH_IDX")"
BATCH_DIR="$BATCHES_DIR/$BATCH_ID"

mkdir -p "$LOG_DIR"
NODE_TMP="${NODE_TMP_ROOT:-/dev/shm}/eq_gapsmith_${SLURM_JOB_ID:-local}_${BATCH_ID}"
mkdir -p "$NODE_TMP"
trap 'rm -rf "$NODE_TMP"' EXIT

DATA_DIR="$GAPSMITH_DATA_MERGED"
if [ -n "${GAPSMITH_DATA_TAR:-}" ]; then
    # October reference (dereferenced data_merged) unpacked once per node into memory; checked against its sha256.
    [ "$(sha256sum "$GAPSMITH_DATA_TAR" | cut -d' ' -f1)" = "$(cat "$GAPSMITH_DATA_TAR.sha256")" ] || { echo "data_merged.tar sha256 mismatch" >&2; exit 3; }
    mkdir -p "$NODE_TMP/ref" && tar -xf "$GAPSMITH_DATA_TAR" -C "$NODE_TMP/ref"
    DATA_DIR="$NODE_TMP/ref/data_merged"
fi

echo "[$(date -Is)] Starting run_gapsmith_batch for $BATCH_ID on $(hostname) (workers=$JOBS_PER_NODE, threads=$THREADS_PER_GENOME)"

"$PYTHON_BIN" "$SHAHEEN_SCRIPTS_DIR/run_gapsmith_batch.py" \
    --batch-dir "$BATCH_DIR" \
    --gapsmith-bin "$GAPSMITH_BIN" \
    --expected-gapsmith-version "$EXPECTED_GAPSMITH_VERSION" \
    --data-dir "$DATA_DIR" \
    --bacteria-ref "$GAPSMITH_BACTERIA_REF" \
    --archaea-ref "$GAPSMITH_ARCHAEA_REF" \
    --diamond-bin "$DIAMOND_BIN" \
    --workers "$JOBS_PER_NODE" \
    --threads-per-genome "$THREADS_PER_GENOME" \
    --scratch-tmp "$NODE_TMP"

echo "[$(date -Is)] Completed run_gapsmith_batch for $BATCH_ID"
