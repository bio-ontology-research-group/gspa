#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=gspa-followup-rescore
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --no-requeue
#SBATCH --exclude=nid00022,nid00497,nid01033
#SBATCH --output=/scratch/project/k10018/gspa-v3/out/stage2b-followups-v1/slurm/rescore-%j.log
#
# Experiments 1 and 3 (design tasks/v3-stage2b-followups-design-2026-10-31.md): re-score the stored off-growth
# outputs under the competitive assignment rules (exp 1) and the per-genome transporter panels (exp 3).
# No linear programs. Pass REPO_ROOT explicitly (Slurm runs a spooled copy of this script).
set -euo pipefail
: "${REPO_ROOT:?REPO_ROOT must point at the code checkout}"
OUT_ROOT=${OUT_ROOT:-/scratch/project/k10018/gspa-v3/out/stage2b-followups-v1}
OFFGROWTH=${OFFGROWTH:-/scratch/project/k10018/gspa-v3/out/offgrowth-v1}
MODELS_TAR=${MODELS_TAR:-/scratch/project/k10018/gspa-v3/benchmark-models-v1/batches/batch_000/outputs.tar}
LAYOUT_DIR=${LAYOUT_DIR:-/project/k10018/gspa-v3/price-full-matrix-v2/layout}
REF_DIR=${REF_DIR:-/project/k10018/gspa-v3/gapsmith-reference-bacteria-v1/reference}
KEGG_MODULES=${KEGG_MODULES:-/project/k10018/gspa-v3/price-production-v1/reference/kegg_modules.tsv}
PANEL_TSV=${PANEL_TSV:-$REPO_ROOT/benchmark/v3/offgrowth/data/biolog_pm1_pm4_panel.tsv}
SWISSPROT_TSV=${SWISSPROT_TSV:-$REPO_ROOT/benchmark/v3/offgrowth/data/swissprot_experimental_ecs_price22.tsv}
PYTHON_BIN=${PYTHON_BIN:-/project/k10018/sw/envs/gspa-v3/bin/python3}
EXPERIMENTS=${EXPERIMENTS:-exp1 exp3}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="$REPO_ROOT/benchmark/v3"
JOB=${SLURM_JOB_ID:-local}
WORK=/dev/shm/followup-rescore-$JOB
mkdir -p "$WORK" "$OUT_ROOT/slurm"
trap 'rm -rf "$WORK"' EXIT
echo "job $JOB on $(hostname) repo=$REPO_ROOT"
sha256sum "$REPO_ROOT/tasks/v3-stage2b-followups-design-2026-10-31.md" "$PANEL_TSV" "$KEGG_MODULES" "$SWISSPROT_TSV" > "$OUT_ROOT/slurm/rescore-$JOB.inputs.sha256"
(cd "$REPO_ROOT/benchmark/v3" && "$PYTHON_BIN" -m unittest followups.tests.test_followups 2>&1 | tail -3)
for exp in $EXPERIMENTS; do
  out=$OUT_ROOT/$exp-$JOB
  mkdir -p "$out"
  case $exp in
    exp1) mod=followups.exp1_assignment; extra=(--layout-dir "$LAYOUT_DIR") ;;
    exp3) mod=followups.exp3_panel; extra=() ;;
    *) echo "unknown experiment $exp" >&2; exit 2 ;;
  esac
  echo "$(date --iso-8601=seconds) start $exp"
  /usr/bin/time -v -o "$out/time_v.txt" "$PYTHON_BIN" -m "$mod" \
    --genomes-dir "$OFFGROWTH/genomes" --candidates-dir "$OFFGROWTH/candidates" --models-tar "$MODELS_TAR" \
    --seed-reactions "$REF_DIR/seed_reactions_corrected.tsv" --ec-aliases "$REF_DIR/seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv" \
    --kegg-modules "$KEGG_MODULES" --panel-tsv "$PANEL_TSV" --swissprot-tsv "$SWISSPROT_TSV" \
    --workers 22 --tmp-dir "$WORK" --out-dir "$out" "${extra[@]}" > "$out/run.log" 2>&1
  echo "$(date --iso-8601=seconds) done $exp (exit $?)"
done
echo "exit 0" > "$OUT_ROOT/slurm/rescore-$JOB.exit"
