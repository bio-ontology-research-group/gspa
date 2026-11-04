#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --exclude=nid00022,nid00497,nid01033
#SBATCH --time=06:00:00
#SBATCH --no-requeue
#SBATCH --job-name=offgrowth-full-aggregate
# Full-panel off-growth run (tasks/v3-offgrowth-full-panel-2026-11-02.md): cohort aggregation of all genome archives
# (harm against Price-148 + every Swiss-Prot experimental label of the 145 genomes), then the secondary tables.
# Submit: sbatch --output=$OUT_ROOT/slurm/aggregate-%j.log --export=ALL,REPO_ROOT=...,OUT_ROOT=...,PILOT_SUMMARY=... <this file>
set -euo pipefail
: "${REPO_ROOT:?}"; : "${OUT_ROOT:?}"; : "${PILOT_SUMMARY:?}"
REF=/project/k10018/gspa-v3/gapsmith-reference-bacteria-v1/reference
PY=/project/k10018/sw/envs/gspa-v3/bin/python3
export PYTHONPATH=$REPO_ROOT/benchmark/v3 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
R=$OUT_ROOT/report
mkdir -p "$R/tables"
sha256sum "$REPO_ROOT/benchmark/v3/offgrowth/"*.py "$OUT_ROOT/inputs/"* > "$R/code_inputs.sha256"
/usr/bin/time -v -o "$R/time_v.txt" "$PY" -m offgrowth.runner aggregate --results-dir "$OUT_ROOT/genomes" \
  --candidates-dir "$OUT_ROOT/candidates" --seed-reactions "$REF/seed_reactions_corrected.tsv" \
  --ec-aliases "$REF/seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv" \
  --unmapped-tsv "$REPO_ROOT/benchmark/v3/offgrowth/data/biolog_pm1_pm4_unmapped.tsv" \
  --swissprot-labels "$OUT_ROOT/inputs/swissprot_experimental_ecs_full145.tsv" --cases-only all --out-dir "$R" \
  > "$R/aggregate.log" 2>&1
/usr/bin/time -v -o "$R/tables/time_v.txt" "$PY" -m offgrowth.full_panel_tables --summary "$R/summary.json" \
  --genomes "$OUT_ROOT/inputs/genomes_shaheen.tsv" --candidates-dir "$OUT_ROOT/candidates" \
  --seed-reactions "$REF/seed_reactions_corrected.tsv" --ec-aliases "$REF/seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv" \
  --pilot-summary "$PILOT_SUMMARY" --out-dir "$R/tables" > "$R/tables/tables.log" 2>&1
gzip -kf "$R/summary.json"
