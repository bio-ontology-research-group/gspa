#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=gspa-v3-offgrowth
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --exclusive
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH --exclude=nid00022,nid00497,nid01033
#SBATCH --output=/scratch/project/k10018/gspa-v3/out/offgrowth-v1/slurm/%x-%j.log
#
# Pre-registered off-growth enzyme discovery experiment (Price-148, 22 genomes).
#
# Modes (set via MODE=pilot|timing|aggregate):
#   MODE=pilot     (default) One Shaheen node (192 cores) per genome across all 22 Price genomes.
#   MODE=timing    Runs a fixed small subset (default MAX_CASES=3 x all arms on 1 genome,
#                  default GENOME_ID=GCF_000236665.1) and writes per-test LP timings
#                  (`timing_per_test.tsv` and `timing_summary.json`) to size the pilot.
#   MODE=aggregate Merges all per-genome `.outputs.tar` archives, runs the cohort-wide
#                  statistical tests and controls, and writes `report.md` + `summary.json`.
#
# Usage outside SLURM (launcher mode):
#   bash benchmark/v3/shaheen/offgrowth-genome.sbatch                      # prints pilot sbatch commands
#   bash benchmark/v3/shaheen/offgrowth-genome.sbatch --submit             # submits 22 pilot jobs
#   MODE=timing bash benchmark/v3/shaheen/offgrowth-genome.sbatch          # prints timing sbatch command
#   MODE=timing bash benchmark/v3/shaheen/offgrowth-genome.sbatch --submit # submits 1-node timing job
#   MODE=aggregate bash benchmark/v3/shaheen/offgrowth-genome.sbatch       # prints aggregate sbatch command

set -euo pipefail

MODE=${MODE:-pilot}
REPO_ROOT=${REPO_ROOT:-$(cd "$(dirname "$0")/../../.." && pwd)}
MANIFEST=${MANIFEST:-$REPO_ROOT/benchmark/manifests/suite6_price149_manifest.tsv}
OUT_ROOT=${OUT_ROOT:-/scratch/project/k10018/gspa-v3/out/offgrowth-v1}
MODELS_TAR=${MODELS_TAR:-/scratch/project/k10018/gspa-v3/benchmark-models-v1/batches/batch_000/outputs.tar}
RC1_COLLECTION=${RC1_COLLECTION:-/scratch/project/k10018/gspa-v3/out/price-production-rc1/collection}
CANDIDATES_DIR=${CANDIDATES_DIR:-$OUT_ROOT/candidates}
REF_DIR=${REF_DIR:-/project/k10018/gspa-v3/gapsmith-reference-bacteria-v1/reference}
KEGG_MODULES=${KEGG_MODULES:-/project/k10018/gspa-v3/price-production-v1/reference/kegg_modules.tsv}
PANEL_TSV=${PANEL_TSV:-$REPO_ROOT/benchmark/v3/offgrowth/data/biolog_pm1_pm4_panel.tsv}
UNMAPPED_TSV=${UNMAPPED_TSV:-$REPO_ROOT/benchmark/v3/offgrowth/data/biolog_pm1_pm4_unmapped.tsv}
EXTRA_CASES_TSV=${EXTRA_CASES_TSV:-}
CASES_ONLY=${CASES_ONLY:-}
if [ -z "$EXTRA_CASES_TSV" ] && [ "${CASES_ONLY:-}" = "extra" ]; then
    EXTRA_CASES_TSV="$REPO_ROOT/benchmark/v3/offgrowth/data/swissprot_experimental_ecs_price22.tsv"
fi
PYTHON_BIN=${PYTHON_BIN:-/project/k10018/sw/envs/gspa-v3/bin/python3}
MAX_CASES=${MAX_CASES:-3}
MAX_TESTS_PER_ARM=${MAX_TESTS_PER_ARM:-0}

if [ "$MODE" = "timing" ]; then
    WALLTIME=${WALLTIME:-01:00:00}
elif [ "$MODE" = "aggregate" ]; then
    WALLTIME=${WALLTIME:-00:30:00}
else
    WALLTIME=${WALLTIME:-04:00:00}
fi

# ---------------------------------------------------------------------------
# Launcher mode: when not running inside a SLURM job (or when --dry-run/--submit passed)
# ---------------------------------------------------------------------------
if [ -z "${SLURM_JOB_ID:-}" ] || [ "${1:-}" = "--dry-run" ] || [ "${1:-}" = "--submit" ]; then
    submit=0
    case "${1:-}" in
        --submit)     submit=1 ;;
        --dry-run|'') submit=0 ;;
        *)
            echo "usage: [MODE=pilot|timing|aggregate] $0 [--submit|--dry-run]" >&2
            exit 2
            ;;
    esac

    if [ ! -f "$MANIFEST" ]; then
        echo "Manifest not found: $MANIFEST" >&2
        exit 2
    fi

    script_path=$(cd "$(dirname "$0")" && pwd)/$(basename "$0")

    common_flags=(
        --account=k10018
        --partition=workq
        --nodes=1
        --ntasks=1
        --cpus-per-task=192
        --exclusive
        --no-requeue
        --exclude=nid00022,nid00497,nid01033
        --time="$WALLTIME"
    )

    cmds=()
    if [ "$MODE" = "timing" ]; then
        TIMING_GENOME=${GENOME_ID:-GCF_000236665.1}
        cmds+=("mkdir -p $(printf '%q ' "$OUT_ROOT/slurm" "$OUT_ROOT/timing")")
        job_name="offgrowth-timing-${TIMING_GENOME}"
        log_file="$OUT_ROOT/slurm/timing-${TIMING_GENOME}-%j.log"
        cmd="sbatch --parsable ${common_flags[*]} --job-name=$(printf '%q' "$job_name") --output=$(printf '%q' "$log_file") --export=ALL,MODE=timing,GENOME_ID=$(printf '%q' "$TIMING_GENOME"),MAX_CASES=$(printf '%q' "$MAX_CASES"),MAX_TESTS_PER_ARM=$(printf '%q' "$MAX_TESTS_PER_ARM"),REPO_ROOT=$(printf '%q' "$REPO_ROOT"),MANIFEST=$(printf '%q' "$MANIFEST"),OUT_ROOT=$(printf '%q' "$OUT_ROOT"),MODELS_TAR=$(printf '%q' "$MODELS_TAR"),RC1_COLLECTION=$(printf '%q' "$RC1_COLLECTION"),CANDIDATES_DIR=$(printf '%q' "$CANDIDATES_DIR"),REF_DIR=$(printf '%q' "$REF_DIR"),KEGG_MODULES=$(printf '%q' "$KEGG_MODULES"),PANEL_TSV=$(printf '%q' "$PANEL_TSV"),UNMAPPED_TSV=$(printf '%q' "$UNMAPPED_TSV"),PYTHON_BIN=$(printf '%q' "$PYTHON_BIN") $(printf '%q' "$script_path")"
        cmds+=("$cmd")
    elif [ "$MODE" = "aggregate" ]; then
        REPORT_DIR=${REPORT_DIR:-$OUT_ROOT/report}
        cmds+=("mkdir -p $(printf '%q ' "$OUT_ROOT/slurm" "$REPORT_DIR")")
        job_name="offgrowth-aggregate"
        log_file="$OUT_ROOT/slurm/aggregate-%j.log"
        export_vars="ALL,MODE=aggregate,REPO_ROOT=$(printf '%q' "$REPO_ROOT"),MANIFEST=$(printf '%q' "$MANIFEST"),OUT_ROOT=$(printf '%q' "$OUT_ROOT"),REPORT_DIR=$(printf '%q' "$REPORT_DIR"),MODELS_TAR=$(printf '%q' "$MODELS_TAR"),RC1_COLLECTION=$(printf '%q' "$RC1_COLLECTION"),CANDIDATES_DIR=$(printf '%q' "$CANDIDATES_DIR"),REF_DIR=$(printf '%q' "$REF_DIR"),KEGG_MODULES=$(printf '%q' "$KEGG_MODULES"),PANEL_TSV=$(printf '%q' "$PANEL_TSV"),UNMAPPED_TSV=$(printf '%q' "$UNMAPPED_TSV"),PYTHON_BIN=$(printf '%q' "$PYTHON_BIN")"
        if [ -n "${EXTRA_CASES_TSV:-}" ]; then
            export_vars="${export_vars},EXTRA_CASES_TSV=$(printf '%q' "$EXTRA_CASES_TSV")"
        fi
        if [ -n "${CASES_ONLY:-}" ]; then
            export_vars="${export_vars},CASES_ONLY=$(printf '%q' "$CASES_ONLY")"
        fi
        cmd="sbatch --parsable ${common_flags[*]} --job-name=$(printf '%q' "$job_name") --output=$(printf '%q' "$log_file") --export=${export_vars} $(printf '%q' "$script_path")"
        cmds+=("$cmd")
    elif [ "$MODE" = "pilot" ]; then
        if [ -n "${GENOME_ID:-}" ]; then
            GENOMES=("$GENOME_ID")
        elif [ -n "${GENOMES:-}" ]; then
            read -r -a GENOMES <<< "$GENOMES"
        elif [ "${CASES_ONLY:-}" = "extra" ] && [ -n "${EXTRA_CASES_TSV:-}" ] && [ -f "$EXTRA_CASES_TSV" ]; then
            mapfile -t GENOMES < <(awk -F'\t' 'NR > 1 && $1 != "" { print $1 }' "$EXTRA_CASES_TSV" | sort -u)
        else
            mapfile -t GENOMES < <(awk -F'\t' 'NR > 1 && $2 != "" && $3 !~ /^(NC|NZ|AC|NT|NW)_/ { print $2 }' "$MANIFEST" | sort -u)
        fi
        cmds+=("mkdir -p $(printf '%q ' "$OUT_ROOT/slurm" "$OUT_ROOT/genomes")")
        for gid in "${GENOMES[@]}"; do
            if [ "${CASES_ONLY:-}" = "extra" ]; then
                job_name="offgrowth-sp-${gid}"
                log_file="$OUT_ROOT/slurm/sp-${gid}-%j.log"
            else
                job_name="offgrowth-${gid}"
                log_file="$OUT_ROOT/slurm/${gid}-%j.log"
            fi
            export_vars="ALL,MODE=pilot,GENOME_ID=$(printf '%q' "$gid"),REPO_ROOT=$(printf '%q' "$REPO_ROOT"),MANIFEST=$(printf '%q' "$MANIFEST"),OUT_ROOT=$(printf '%q' "$OUT_ROOT"),MODELS_TAR=$(printf '%q' "$MODELS_TAR"),RC1_COLLECTION=$(printf '%q' "$RC1_COLLECTION"),CANDIDATES_DIR=$(printf '%q' "$CANDIDATES_DIR"),REF_DIR=$(printf '%q' "$REF_DIR"),KEGG_MODULES=$(printf '%q' "$KEGG_MODULES"),PANEL_TSV=$(printf '%q' "$PANEL_TSV"),UNMAPPED_TSV=$(printf '%q' "$UNMAPPED_TSV"),PYTHON_BIN=$(printf '%q' "$PYTHON_BIN")"
            if [ -n "${EXTRA_CASES_TSV:-}" ]; then
                export_vars="${export_vars},EXTRA_CASES_TSV=$(printf '%q' "$EXTRA_CASES_TSV")"
            fi
            if [ -n "${CASES_ONLY:-}" ]; then
                export_vars="${export_vars},CASES_ONLY=$(printf '%q' "$CASES_ONLY")"
            fi
            cmd="sbatch --parsable ${common_flags[*]} --job-name=$(printf '%q' "$job_name") --output=$(printf '%q' "$log_file") --export=${export_vars} $(printf '%q' "$script_path")"
            cmds+=("$cmd")
        done
    else
        echo "Unsupported MODE='$MODE' (expected pilot, timing, or aggregate)" >&2
        exit 2
    fi

    if [ "$submit" -eq 1 ]; then
        command -v sbatch >/dev/null 2>&1 || {
            echo "sbatch not found: run with --submit only on a Shaheen login node" >&2
            exit 2
        }
        printf '%s\n' "${cmds[@]}" | bash
    else
        echo "# Dry-run (default): MODE=$MODE (${#cmds[@]} commands)."
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

if [ "$MODE" = "aggregate" ]; then
    REPORT_DIR=${REPORT_DIR:-$OUT_ROOT/report}
    mkdir -p "$REPORT_DIR" "$OUT_ROOT/slurm"
    AGG_EXTRA=()
    if [ -n "${CASES_ONLY:-}" ]; then
        AGG_EXTRA+=(--cases-only "$CASES_ONLY")
    fi
    SP_LABELS=${SWISSPROT_LABELS:-${EXTRA_CASES_TSV:-$REPO_ROOT/benchmark/v3/offgrowth/data/swissprot_experimental_ecs_price22.tsv}}
    if [ -f "$SP_LABELS" ]; then
        AGG_EXTRA+=(--swissprot-labels "$SP_LABELS")
    fi
    /usr/bin/time -v -o "$REPORT_DIR/aggregate_time_v.txt" \
        "$PYTHON_BIN" -m offgrowth.runner aggregate \
        --results-dir "$OUT_ROOT/genomes" \
        --candidates-dir "$CANDIDATES_DIR" \
        --seed-reactions "$REF_DIR/seed_reactions_corrected.tsv" \
        --ec-aliases "$REF_DIR/seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv" \
        --unmapped-tsv "$UNMAPPED_TSV" \
        --out-dir "$REPORT_DIR" \
        "${AGG_EXTRA[@]}" \
        > "$REPORT_DIR/aggregate.log" 2>&1
    exit 0
fi

: "${GENOME_ID:?GENOME_ID must be set in worker mode for MODE=$MODE}"

WORK_DIR="/dev/shm/offgrowth-${SLURM_JOB_ID}-${GENOME_ID}"
STAGE_OUT="$WORK_DIR/out/${GENOME_ID}"
mkdir -p "$WORK_DIR/models" "$WORK_DIR/candidates" "$WORK_DIR/ref" "$STAGE_OUT" "$OUT_ROOT/genomes" "$OUT_ROOT/timing" "$OUT_ROOT/slurm"

cleanup() {
    local code=$?
    trap - EXIT
    set +e
    if [ "$code" -ne 0 ] && [ -d "$STAGE_OUT" ]; then
        # Keep the evidence of a failed run: logs and receipts (no bulk outputs) go to scratch.
        mkdir -p "$OUT_ROOT/failed/${SLURM_JOB_ID}-${GENOME_ID}"
        find "$STAGE_OUT" -maxdepth 1 -type f \( -name '*.txt' -o -name '*.log' -o -name '*.sha256' -o -name '*.json' \) \
            -exec cp {} "$OUT_ROOT/failed/${SLURM_JOB_ID}-${GENOME_ID}/" \;
        echo "exit $code; logs kept in $OUT_ROOT/failed/${SLURM_JOB_ID}-${GENOME_ID}" >&2
    fi
    rm -rf "$WORK_DIR"
    exit "$code"
}
trap 'exit 143' TERM INT
trap cleanup EXIT

date --iso-8601=seconds > "$STAGE_OUT/started_at.txt"
hostname > "$STAGE_OUT/hostname.txt"
cp "$0" "$STAGE_OUT/executed.sbatch"

# 1. Stage reference tables and this genome's gapsmith model members into /dev/shm
cp -f "$MANIFEST" "$WORK_DIR/ref/manifest.tsv"
cp -f "$PANEL_TSV" "$WORK_DIR/ref/biolog_pm1_pm4_panel.tsv"
cp -f "$UNMAPPED_TSV" "$WORK_DIR/ref/biolog_pm1_pm4_unmapped.tsv"
cp -f "$REF_DIR/seed_metabolites_edited.tsv" "$WORK_DIR/ref/seed_metabolites_edited.tsv"
cp -f "$REF_DIR/seed_reactions_corrected.tsv" "$WORK_DIR/ref/seed_reactions_corrected.tsv"
cp -f "$REF_DIR/seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv" "$WORK_DIR/ref/seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv"
cp -f "$KEGG_MODULES" "$WORK_DIR/ref/kegg_modules.tsv"
EXTRA_CASES_ARG=()
EXTRA_HASH_FILES=()
if [ -n "${EXTRA_CASES_TSV:-}" ] && [ -f "$EXTRA_CASES_TSV" ]; then
    cp -f "$EXTRA_CASES_TSV" "$WORK_DIR/ref/extra_cases.tsv"
    EXTRA_CASES_ARG=(--extra-cases-tsv "$WORK_DIR/ref/extra_cases.tsv")
    EXTRA_HASH_FILES=("$WORK_DIR/ref/extra_cases.tsv")
fi

tar -xf "$MODELS_TAR" -C "$WORK_DIR/models" "${GENOME_ID}/"

# 2. Stage compact candidate TSV into /dev/shm if present; otherwise resolve raw stage1 TSVs via plan.json
CAND_INPUTS=()
if [ -f "$CANDIDATES_DIR/${GENOME_ID}.candidates.tsv" ]; then
    cp -f "$CANDIDATES_DIR/${GENOME_ID}.candidates.tsv" "$WORK_DIR/candidates/${GENOME_ID}.candidates.tsv"
    CAND_INPUTS=("$WORK_DIR/candidates/${GENOME_ID}.candidates.tsv")
else
    mapfile -t CAND_INPUTS < <(find "$RC1_COLLECTION/attempts" -type f -path "*/production/stage1/${GENOME_ID}.integrated.tsv" | sort)
    if [ "${#CAND_INPUTS[@]}" -eq 0 ]; then
        echo "No compact candidate TSV or stage1 integrated TSVs found for $GENOME_ID" >&2
        exit 3
    fi
fi

# 3. Record input SHA-256 receipt and software/solver versions
sha256sum \
    "$WORK_DIR/ref/manifest.tsv" \
    "$WORK_DIR/ref/biolog_pm1_pm4_panel.tsv" \
    "$WORK_DIR/ref/biolog_pm1_pm4_unmapped.tsv" \
    "$WORK_DIR/ref/seed_metabolites_edited.tsv" \
    "$WORK_DIR/ref/seed_reactions_corrected.tsv" \
    "$WORK_DIR/ref/seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv" \
    "$WORK_DIR/ref/kegg_modules.tsv" \
    "${EXTRA_HASH_FILES[@]}" \
    "$WORK_DIR/models/${GENOME_ID}/${GENOME_ID}-draft.xml" \
    "$WORK_DIR/models/${GENOME_ID}/${GENOME_ID}-filled.xml" \
    "$WORK_DIR/models/${GENOME_ID}/${GENOME_ID}-medium.csv" \
    "$WORK_DIR/models/${GENOME_ID}/${GENOME_ID}-Transporter.tbl" \
    "${CAND_INPUTS[@]}" \
    > "$STAGE_OUT/inputs.sha256"

"$PYTHON_BIN" -c '
import sys, cobra, optlang
from offgrowth.arms import detect_solvers
print("python:", sys.version)
print("cobra:", cobra.__version__)
print("optlang:", optlang.__version__)
print("solvers:", detect_solvers())
' > "$STAGE_OUT/versions.txt" 2>&1

# SLURM_CPUS_ON_NODE counts hardware threads (384 on workq); one LP worker per physical core.
WORKERS=${WORKERS:-${SLURM_CPUS_PER_TASK:-192}}
if [ -f "$WORK_DIR/candidates/${GENOME_ID}.candidates.tsv" ]; then
    CAND_ARGS=(--candidates-tsv "$WORK_DIR/candidates/${GENOME_ID}.candidates.tsv")
else
    CAND_ARGS=(--stage1-tsvs "${CAND_INPUTS[@]}")
fi
EXTRA_ARGS=(--workers "$WORKERS")
if [ -n "${CASES_ONLY:-}" ]; then
    EXTRA_ARGS+=(--cases-only "$CASES_ONLY")
fi
if [ "$MODE" = "timing" ]; then
    EXTRA_ARGS+=(--run-mode timing --max-cases "$MAX_CASES")
    if [ "${MAX_TESTS_PER_ARM:-0}" -gt 0 ]; then
        EXTRA_ARGS+=(--max-tests-per-arm "$MAX_TESTS_PER_ARM")
    fi
else
    EXTRA_ARGS+=(--run-mode pilot)
fi

# 4. Run off-growth discovery for GENOME_ID under /usr/bin/time -v
/usr/bin/time -v -o "$STAGE_OUT/time_v.txt" \
    "$PYTHON_BIN" -m "${RUNNER_MODULE:-offgrowth.runner}" genome \
    --genome-id "$GENOME_ID" \
    --manifest "$WORK_DIR/ref/manifest.tsv" \
    --genome-model-dir "$WORK_DIR/models/${GENOME_ID}" \
    "${CAND_ARGS[@]}" \
    --seed-metabolites "$WORK_DIR/ref/seed_metabolites_edited.tsv" \
    --seed-reactions "$WORK_DIR/ref/seed_reactions_corrected.tsv" \
    --ec-aliases "$WORK_DIR/ref/seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv" \
    --panel-tsv "$WORK_DIR/ref/biolog_pm1_pm4_panel.tsv" \
    --kegg-modules "$WORK_DIR/ref/kegg_modules.tsv" \
    "${EXTRA_CASES_ARG[@]}" \
    --out-dir "$WORK_DIR/out" \
    "${EXTRA_ARGS[@]}" \
    2>&1 | tee "$STAGE_OUT/runner.log" > "$OUT_ROOT/slurm/runner-${SLURM_JOB_ID}-${GENOME_ID}.log"

date --iso-8601=seconds > "$STAGE_OUT/finished_at.txt"

# 5. Pack outputs into a single tar archive on /scratch (and copy timing TSV/JSON in MODE=timing)
if [ "$MODE" = "timing" ]; then
    mkdir -p "$OUT_ROOT/timing/${GENOME_ID}"
    cp -f "$STAGE_OUT/timing_per_test.tsv" "$OUT_ROOT/timing/${GENOME_ID}/"
    cp -f "$STAGE_OUT/timing_summary.json" "$OUT_ROOT/timing/${GENOME_ID}/"
    cp -f "$STAGE_OUT/time_v.txt" "$OUT_ROOT/timing/${GENOME_ID}/"
    cp -f "$STAGE_OUT/versions.txt" "$OUT_ROOT/timing/${GENOME_ID}/"
    cp -f "$STAGE_OUT/inputs.sha256" "$OUT_ROOT/timing/${GENOME_ID}/"
    TAR_DEST="$OUT_ROOT/timing/${GENOME_ID}.timing.outputs.tar"
elif [ "${CASES_ONLY:-}" = "extra" ]; then
    TAR_DEST="$OUT_ROOT/genomes/${GENOME_ID}.swissprot.outputs.tar"
else
    TAR_DEST="$OUT_ROOT/genomes/${GENOME_ID}.outputs.tar"
fi

tar -cf "${TAR_DEST}.partial" -C "$WORK_DIR/out" "${GENOME_ID}"
mv "${TAR_DEST}.partial" "$TAR_DEST"
sha256sum "$TAR_DEST" > "${TAR_DEST}.sha256"
