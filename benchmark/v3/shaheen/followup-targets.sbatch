#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=gspa-followup-targets
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --exclusive
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH --exclude=nid00022,nid00497,nid01033
#SBATCH --output=/scratch/project/k10018/gspa-v3/out/stage2b-followups-v1/slurm/targets-%j.log
#
# Experiment 2: unblocking / dead-end tests as a gated evidence source; writes sparse update records per cell and arm
# from the stored unmasked A2/A3 runs and each predictor's stage-0 output. No linear programs.
set -euo pipefail
: "${REPO_ROOT:?REPO_ROOT must point at the code checkout}"
OUT_ROOT=${OUT_ROOT:-/scratch/project/k10018/gspa-v3/out/stage2b-followups-v1}
COLLECTION=${COLLECTION:-/scratch/project/k10018/gspa-v3/out/price-stage2b-rc1/collection}
OFFGROWTH=${OFFGROWTH:-/scratch/project/k10018/gspa-v3/out/offgrowth-v1}
MODELS_TAR=${MODELS_TAR:-/scratch/project/k10018/gspa-v3/benchmark-models-v1/batches/batch_000/outputs.tar}
LAYOUT_DIR=${LAYOUT_DIR:-/project/k10018/gspa-v3/price-full-matrix-v2/layout}
REF_DIR=${REF_DIR:-/project/k10018/gspa-v3/gapsmith-reference-bacteria-v1/reference}
KEGG_MODULES=${KEGG_MODULES:-/project/k10018/gspa-v3/price-production-v1/reference/kegg_modules.tsv}
PANEL_TSV=${PANEL_TSV:-$REPO_ROOT/benchmark/v3/offgrowth/data/biolog_pm1_pm4_panel.tsv}
SWISSPROT_TSV=${SWISSPROT_TSV:-$REPO_ROOT/benchmark/v3/offgrowth/data/swissprot_experimental_ecs_price22.tsv}
PYTHON_BIN=${PYTHON_BIN:-/project/k10018/sw/envs/gspa-v3/bin/python3}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="$REPO_ROOT/benchmark/v3"
JOB=${SLURM_JOB_ID:-local}
WORK=/dev/shm/followup-targets-$JOB
OUT=$OUT_ROOT/exp2-$JOB
mkdir -p "$WORK" "$OUT_ROOT/slurm" "$OUT/updates"
trap 'rm -rf "$WORK"' EXIT
echo "job $JOB on $(hostname) repo=$REPO_ROOT"
/usr/bin/time -v -o "$OUT/time_v.txt" "$PYTHON_BIN" -m followups.exp2_targets \
  --collection "$COLLECTION" --genomes-dir "$OFFGROWTH/genomes" --layout-dir "$LAYOUT_DIR" --models-tar "$MODELS_TAR" \
  --seed-reactions "$REF_DIR/seed_reactions_corrected.tsv" --ec-aliases "$REF_DIR/seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv" \
  --kegg-modules "$KEGG_MODULES" --panel-tsv "$PANEL_TSV" --swissprot-tsv "$SWISSPROT_TSV" \
  --workers 22 --tmp-dir "$WORK" --out-dir "$OUT/updates" > "$OUT/run.log" 2>&1
echo "exit $?" > "$OUT/exit-code.txt"
