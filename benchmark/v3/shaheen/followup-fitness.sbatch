#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=gspa-followup-fitness
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --exclusive
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH --exclude=nid00022,nid00497,nid01033
#SBATCH --output=/scratch/project/k10018/gspa-v3/out/stage2b-followups-v1/slurm/fitness-%j.log
#
# Experiment 6: knockout growth of the model variants on the matched panel wells for one genome (GENOME_ID), scored
# against the Fitness Browser specific phenotypes prepared locally (FITNESS_DIR/<genome>.fitness.json).
# MODE=aggregate pools the per-genome results instead.
set -euo pipefail
: "${REPO_ROOT:?REPO_ROOT must point at the code checkout}"
OUT_ROOT=${OUT_ROOT:-/scratch/project/k10018/gspa-v3/out/stage2b-followups-v1}
FITNESS_DIR=${FITNESS_DIR:-/project/k10018/gspa-v3/fitness-browser-v1/fitness}
OFFGROWTH=${OFFGROWTH:-/scratch/project/k10018/gspa-v3/out/offgrowth-v1}
MODELS_TAR=${MODELS_TAR:-/scratch/project/k10018/gspa-v3/benchmark-models-v1/batches/batch_000/outputs.tar}
LAYOUT_DIR=${LAYOUT_DIR:-/project/k10018/gspa-v3/price-full-matrix-v2/layout}
REF_DIR=${REF_DIR:-/project/k10018/gspa-v3/gapsmith-reference-bacteria-v1/reference}
KEGG_MODULES=${KEGG_MODULES:-/project/k10018/gspa-v3/price-production-v1/reference/kegg_modules.tsv}
PANEL_TSV=${PANEL_TSV:-$REPO_ROOT/benchmark/v3/offgrowth/data/biolog_pm1_pm4_panel.tsv}
PYTHON_BIN=${PYTHON_BIN:-/project/k10018/sw/envs/gspa-v3/bin/python3}
MODE=${MODE:-genome}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="$REPO_ROOT/benchmark/v3"
JOB=${SLURM_JOB_ID:-local}
OUT=$OUT_ROOT/exp6
mkdir -p "$OUT/genomes" "$OUT/logs" "$OUT_ROOT/slurm"
if [ "$MODE" = "aggregate" ]; then
  "$PYTHON_BIN" -m followups.exp6_aggregate --results-dir "$OUT/genomes" --out-dir "$OUT/report-$JOB" > "$OUT/logs/aggregate-$JOB.log" 2>&1
  exit $?
fi
: "${GENOME_ID:?GENOME_ID required}"
WORK=/dev/shm/followup-fitness-$JOB-$GENOME_ID
mkdir -p "$WORK/models"
trap 'rm -rf "$WORK"' EXIT
tar -xf "$MODELS_TAR" -C "$WORK/models" "${GENOME_ID}/"
echo "job $JOB on $(hostname) genome=$GENOME_ID"
/usr/bin/time -v -o "$OUT/logs/$GENOME_ID-$JOB.time_v.txt" "$PYTHON_BIN" -m followups.exp6_fitness \
  --genome-id "$GENOME_ID" --fitness-json "$FITNESS_DIR/$GENOME_ID.fitness.json" --genome-model-dir "$WORK/models/$GENOME_ID" \
  --genomes-dir "$OFFGROWTH/genomes" --candidates-tsv "$OFFGROWTH/candidates/$GENOME_ID.candidates.tsv" --layout-tsv "$LAYOUT_DIR/$GENOME_ID.layout.tsv" \
  --seed-reactions "$REF_DIR/seed_reactions_corrected.tsv" --ec-aliases "$REF_DIR/seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv" \
  --seed-metabolites "$REF_DIR/seed_metabolites_edited.tsv" --kegg-modules "$KEGG_MODULES" --panel-tsv "$PANEL_TSV" \
  --workers 192 --out-dir "$OUT/genomes" > "$OUT/logs/$GENOME_ID-$JOB.log" 2>&1
code=$?
echo "exit $code" > "$OUT/logs/$GENOME_ID-$JOB.exit"
exit $code
