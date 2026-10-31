#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=gspa-v3-biology-extract
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --no-requeue
#SBATCH --exclude=nid00022,nid00497,nid01033
#SBATCH --output=/scratch/project/k10018/gspa-v3/out/biology-rc1/slurm/%x-%j.log
#
# Parallel biology summary extraction, rule checking, and Markdown reporting
# across all 132 cells of the release-candidate Price collection on Shaheen III.
#
# Usage outside SLURM (launcher mode):
#   bash benchmark/v3/shaheen/biology-extract.sbatch           # prints commands (default dry-run)
#   bash benchmark/v3/shaheen/biology-extract.sbatch --submit  # submits job to Slurm
#
# Usage inside SLURM (worker mode):
#   Streams stage0 and stage1 <assembly>.integrated.tsv files for all 132 cells
#   in parallel, writes compact per-cell JSON summaries and input SHA-256 receipts
#   to /scratch/project/k10018/gspa-v3/out/biology-rc1, runs checks.py and
#   report.py, and packs biology-rc1-summaries.tar.gz for transfer back to
#   benchmark/v3/results/biology-rc1/.

set -euo pipefail

REPO_ROOT=${REPO_ROOT:-$(cd "$(dirname "$0")/../../.." && pwd)}
COLLECTION=${COLLECTION:-/scratch/project/k10018/gspa-v3/out/price-production-rc1/collection}
OUT_ROOT=${OUT_ROOT:-/scratch/project/k10018/gspa-v3/out/biology-rc1}
PYTHON_BIN=${PYTHON_BIN:-/project/k10018/sw/envs/gspa-tools/bin/python3.12}
WORKERS=${WORKERS:-64}
WALLTIME=${WALLTIME:-02:00:00}

# ---------------------------------------------------------------------------
# Launcher mode: default prints commands; submits only when --submit is passed
# ---------------------------------------------------------------------------
if [ -z "${SLURM_JOB_ID:-}" ] || [ "${1:-}" = "--dry-run" ] || [ "${1:-}" = "--submit" ]; then
    submit=0
    case "${1:-}" in
        --submit)
            submit=1
            ;;
        --dry-run|'')
            submit=0
            ;;
        *)
            echo "usage: $0 [--submit|--dry-run]" >&2
            exit 2
            ;;
    esac

    script_path=$(cd "$(dirname "$0")" && pwd)/$(basename "$0")
    log_file="$OUT_ROOT/slurm/gspa-v3-biology-extract-%j.log"

    sbatch_flags=(
        --account=k10018
        --partition=workq
        --job-name=gspa-v3-biology-extract
        --nodes=1
        --ntasks=1
        --cpus-per-task=192
        --exclusive
        --no-requeue
        --exclude=nid00022,nid00497,nid01033
        --time="$WALLTIME"
        --output="$log_file"
        --export="ALL,REPO_ROOT=$REPO_ROOT,COLLECTION=$COLLECTION,OUT_ROOT=$OUT_ROOT,PYTHON_BIN=$PYTHON_BIN,WORKERS=$WORKERS"
    )

    mkdir_cmd="mkdir -p $(printf '%q ' "$OUT_ROOT/slurm" "$OUT_ROOT/cells" "$OUT_ROOT/receipts")"
    submit_cmd="sbatch --parsable $(printf '%q ' "${sbatch_flags[@]}") $(printf '%q' "$script_path")"

    if [ "$submit" -eq 1 ]; then
        command -v sbatch >/dev/null 2>&1 || {
            echo "sbatch not found: run with --submit only on a Shaheen login node" >&2
            exit 2
        }
        eval "$mkdir_cmd"
        job_id=$(eval "$submit_cmd")
        echo "Submitted job $job_id -> $OUT_ROOT"
    else
        echo "# Dry-run (default): prints commands only. Pass --submit to submit on Shaheen III."
        printf '%s\n' "$mkdir_cmd"
        printf '%s\n' "$submit_cmd"
    fi
    exit 0
fi

# ---------------------------------------------------------------------------
# Worker mode: executing inside a SLURM allocation on a Shaheen compute node
# ---------------------------------------------------------------------------
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="$REPO_ROOT/benchmark/v3${PYTHONPATH:+:$PYTHONPATH}"

JOB_DIR="$OUT_ROOT/jobs/${SLURM_JOB_ID}"
mkdir -p "$OUT_ROOT/cells" "$OUT_ROOT/receipts" "$OUT_ROOT/slurm" "$JOB_DIR"

finish() {
    local code=$?
    trap - EXIT
    printf '%s\n' "$code" > "$JOB_DIR/exit-code.txt"
    date --iso-8601=seconds > "$JOB_DIR/finished-at.txt"
    exit "$code"
}
trap finish EXIT

date --iso-8601=seconds > "$JOB_DIR/started-at.txt"
hostname > "$JOB_DIR/hostname.txt"
cp "$0" "$JOB_DIR/executed.sbatch"

if [ ! -x "$PYTHON_BIN" ]; then
    PYTHON_BIN=$(command -v python3)
fi
"$PYTHON_BIN" --version > "$JOB_DIR/python-version.txt" 2>&1

(
    cd "$REPO_ROOT"
    sha256sum \
        benchmark/v3/ec_identity.py \
        benchmark/v3/biology/__init__.py \
        benchmark/v3/biology/extract.py \
        benchmark/v3/biology/checks.py \
        benchmark/v3/biology/report.py \
        benchmark/v3/shaheen/biology-extract.sbatch \
        benchmark/manifests/suite6_price149_manifest.tsv \
        gspa-core/src/main/resources/taxon-constraints/go-taxon-constraints.tsv \
        gspa-core/src/main/resources/taxon-constraints/ncbi-taxon-hierarchy.tsv
) > "$OUT_ROOT/source.sha256"
cp "$OUT_ROOT/source.sha256" "$JOB_DIR/source.sha256"

# 1. Run unit tests first
"$PYTHON_BIN" -m unittest discover -s "$REPO_ROOT/benchmark/v3/biology/tests" -p 'test_*.py' -v \
    > "$JOB_DIR/unit-tests.log" 2>&1

# 2. Extract compact summaries and input SHA-256 receipts across all 132 cells in parallel
/usr/bin/time -v "$PYTHON_BIN" "$REPO_ROOT/benchmark/v3/biology/extract.py" \
    --collection "$COLLECTION" \
    --out "$OUT_ROOT" \
    --workers "$WORKERS" \
    > "$JOB_DIR/extract.log" 2>&1

# 3. Run fixed biological checks
"$PYTHON_BIN" "$REPO_ROOT/benchmark/v3/biology/checks.py" \
    --summaries "$OUT_ROOT" \
    --out "$OUT_ROOT/evaluation.json" \
    > "$JOB_DIR/checks.log" 2>&1

# 4. Render Markdown report tables and plain-word summary
"$PYTHON_BIN" "$REPO_ROOT/benchmark/v3/biology/report.py" \
    --input "$OUT_ROOT/evaluation.json" \
    --out "$OUT_ROOT/report.md" \
    > "$JOB_DIR/report.log" 2>&1

# 5. Package compact summaries and receipts for transfer to benchmark/v3/results/biology-rc1/
(
    cd "$OUT_ROOT"
    tar -czf biology-rc1-summaries.tar.gz \
        source.sha256 \
        receipts.json \
        evaluation.json \
        report.md \
        cells \
        receipts
    sha256sum biology-rc1-summaries.tar.gz > biology-rc1-summaries.tar.gz.sha256
)
