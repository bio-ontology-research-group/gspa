#!/usr/bin/env bash
#SBATCH --job-name=eq_stage_batch
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --time=24:00:00
#SBATCH --exclusive
#SBATCH --exclude=nid00022
#SBATCH --no-requeue

# Stage 1: selectively extract needed bin FASTA + Bakta .faa/.gff3/.fna from site tars
# (or RAK-PR / benchmark inputs), run pyrodigal 3.7.1 -p single where protein_source == pyrodigal,
# strictly verify n_proteins against eq_mags_per_bin.tsv, and pack verified inputs into
# <WORK_ROOT>/batches/batch_NNN/inputs.tar.

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
NODE_TMP="${NODE_TMP_ROOT:-/dev/shm}/eq_stage_${SLURM_JOB_ID:-local}_${BATCH_ID}"
mkdir -p "$NODE_TMP"
trap 'rm -rf "$NODE_TMP"' EXIT

echo "[$(date -Is)] Starting stage_batch for $BATCH_ID on $(hostname)"

"$PYTHON_BIN" "$SHAHEEN_SCRIPTS_DIR/stage_batch.py" \
    --batch-dir "$BATCH_DIR" \
    --pyrodigal-bin "$PYRODIGAL_BIN" \
    --expected-pyrodigal-version "$EXPECTED_PYRODIGAL_VERSION" \
    --workers "$JOBS_PER_NODE" \
    --scratch-tmp "$NODE_TMP" \
    --ensure-merged-ref \
    --bacteria-ref "$GAPSMITH_BACTERIA_REF" \
    --archaea-ref "$GAPSMITH_ARCHAEA_REF" \
    --merged-ref "$GAPSMITH_DATA_MERGED"

echo "[$(date -Is)] Completed stage_batch for $BATCH_ID"
