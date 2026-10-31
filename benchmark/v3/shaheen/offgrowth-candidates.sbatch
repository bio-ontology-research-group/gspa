#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=gspa-v3-offgrowth-candidates
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --exclusive
#SBATCH --time=00:45:00
#SBATCH --no-requeue
#SBATCH --exclude=nid00022,nid00497,nid01033
#SBATCH --output=/scratch/project/k10018/gspa-v3/out/offgrowth-v1/slurm/%x-%j.log
#
# Pre-pass job: stream the 6 predictors' stage1 `<assembly>.integrated.tsv` files
# from `/scratch/project/k10018/gspa-v3/out/price-production-rc1/collection`
# (mapped via `plan.json` -> `attempts/<idx>/job-*/production/stage1/`) in bounded
# memory across all 22 Price genomes in parallel on one Shaheen III node (192 cores),
# writing one compact candidate TSV (`<genome_id>.candidates.tsv`) and SHA-256
# receipt (`<genome_id>.candidates.receipt.json`) per genome.
#
# Usage outside SLURM (launcher mode):
#   bash benchmark/v3/shaheen/offgrowth-candidates.sbatch           # prints sbatch command (default dry-run)
#   bash benchmark/v3/shaheen/offgrowth-candidates.sbatch --submit  # submits the 1-node pre-pass job

set -euo pipefail

REPO_ROOT=${REPO_ROOT:-$(cd "$(dirname "$0")/../../.." && pwd)}
MANIFEST=${MANIFEST:-$REPO_ROOT/benchmark/manifests/suite6_price149_manifest.tsv}
OUT_ROOT=${OUT_ROOT:-/scratch/project/k10018/gspa-v3/out/offgrowth-v1}
RC1_COLLECTION=${RC1_COLLECTION:-/scratch/project/k10018/gspa-v3/out/price-production-rc1/collection}
CANDIDATES_DIR=${CANDIDATES_DIR:-$OUT_ROOT/candidates}
PYTHON_BIN=${PYTHON_BIN:-/project/k10018/sw/envs/gspa-v3/bin/python3}
MAX_PER_EC=${MAX_PER_EC:-64}
WALLTIME=${WALLTIME:-00:45:00}

# ---------------------------------------------------------------------------
# Launcher mode: when not running inside a SLURM job (or when --dry-run/--submit passed)
# ---------------------------------------------------------------------------
if [ -z "${SLURM_JOB_ID:-}" ] || [ "${1:-}" = "--dry-run" ] || [ "${1:-}" = "--submit" ]; then
    submit=0
    case "${1:-}" in
        --submit)      submit=1 ;;
        --dry-run|'')  submit=0 ;;
        *)
            echo "usage: $0 [--submit|--dry-run]" >&2
            exit 2
            ;;
    esac

    script_path=$(cd "$(dirname "$0")" && pwd)/$(basename "$0")
    log_file="$OUT_ROOT/slurm/offgrowth-candidates-%j.log"

    cmds=(
        "mkdir -p $(printf '%q ' "$OUT_ROOT/slurm" "$CANDIDATES_DIR")"
        "sbatch --parsable --account=k10018 --partition=workq --nodes=1 --ntasks=1 --cpus-per-task=192 --exclusive --no-requeue --exclude=nid00022,nid00497,nid01033 --time=$(printf '%q' "$WALLTIME") --job-name=gspa-v3-offgrowth-candidates --output=$(printf '%q' "$log_file") --export=ALL,REPO_ROOT=$(printf '%q' "$REPO_ROOT"),MANIFEST=$(printf '%q' "$MANIFEST"),OUT_ROOT=$(printf '%q' "$OUT_ROOT"),RC1_COLLECTION=$(printf '%q' "$RC1_COLLECTION"),CANDIDATES_DIR=$(printf '%q' "$CANDIDATES_DIR"),PYTHON_BIN=$(printf '%q' "$PYTHON_BIN"),MAX_PER_EC=$(printf '%q' "$MAX_PER_EC") $(printf '%q' "$script_path")"
    )

    if [ "$submit" -eq 1 ]; then
        command -v sbatch >/dev/null 2>&1 || {
            echo "sbatch not found: run with --submit only on a Shaheen login node" >&2
            exit 2
        }
        printf '%s\n' "${cmds[@]}" | bash
    else
        echo "# Dry-run (default): 1-node parallel candidate pre-pass over 22 Price genomes."
        echo "# Pass --submit to submit to Slurm on Shaheen III."
        printf '%s\n' "${cmds[@]}"
    fi
    exit 0
fi

# ---------------------------------------------------------------------------
# Worker mode: executing inside a SLURM allocation
# ---------------------------------------------------------------------------
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="$REPO_ROOT/benchmark/v3${PYTHONPATH:+:$PYTHONPATH}"

WORK_DIR="/dev/shm/offgrowth-candidates-${SLURM_JOB_ID}"
mkdir -p "$WORK_DIR/candidates" "$CANDIDATES_DIR" "$OUT_ROOT/slurm"

cleanup() {
    local code=$?
    trap - EXIT
    set +e
    rm -rf "$WORK_DIR"
    exit "$code"
}
trap 'exit 143' TERM INT
trap cleanup EXIT

date --iso-8601=seconds > "$WORK_DIR/candidates/started_at.txt"
hostname > "$WORK_DIR/candidates/hostname.txt"
cp "$0" "$WORK_DIR/candidates/executed.sbatch"
sha256sum "$MANIFEST" "$RC1_COLLECTION/plan.json" > "$WORK_DIR/candidates/inputs.sha256"

"$PYTHON_BIN" -c '
import sys, cobra
from offgrowth.arms import detect_solvers
print("python:", sys.version)
print("cobra:", cobra.__version__)
print("solvers:", detect_solvers())
' > "$WORK_DIR/candidates/versions.txt" 2>&1

WORKERS=${SLURM_CPUS_ON_NODE:-192}

/usr/bin/time -v -o "$WORK_DIR/candidates/time_v.txt" \
    "$PYTHON_BIN" -m offgrowth.runner prepass \
    --collection-dir "$RC1_COLLECTION" \
    --manifest "$MANIFEST" \
    --out-dir "$WORK_DIR/candidates" \
    --workers "$WORKERS" \
    --max-per-ec "$MAX_PER_EC" \
    > "$WORK_DIR/candidates/prepass.log" 2>&1

date --iso-8601=seconds > "$WORK_DIR/candidates/finished_at.txt"

# Copy compact candidate TSVs and receipts to $CANDIDATES_DIR and pack an archive receipt
cp -f "$WORK_DIR/candidates/"*.candidates.tsv "$CANDIDATES_DIR/"
cp -f "$WORK_DIR/candidates/"*.candidates.receipt.json "$CANDIDATES_DIR/"
cp -f "$WORK_DIR/candidates/prepass_summary.json" "$CANDIDATES_DIR/"
cp -f "$WORK_DIR/candidates/time_v.txt" "$CANDIDATES_DIR/prepass_time_v.txt"
cp -f "$WORK_DIR/candidates/versions.txt" "$CANDIDATES_DIR/prepass_versions.txt"
cp -f "$WORK_DIR/candidates/inputs.sha256" "$CANDIDATES_DIR/prepass_inputs.sha256"
sha256sum "$CANDIDATES_DIR/"*.candidates.tsv > "$CANDIDATES_DIR/candidates.sha256"
