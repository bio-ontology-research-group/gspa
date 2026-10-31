#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --exclude=nid00022,nid00497,nid01033
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH --job-name=offgrowth-aggregate2
#SBATCH --output=/scratch/project/k10018/gspa-v3/out/offgrowth-v1/slurm/aggregate2-%j.log
# Off-growth cohort merge twice: harm against Price-148 only, and against Price-148 + Swiss-Prot experimental ECs.
set -euo pipefail
REPO=/project/k10018/code/gspa-v3
O=/scratch/project/k10018/gspa-v3/out/offgrowth-v1
REF=/project/k10018/gspa-v3/gapsmith-reference-bacteria-v1/reference
PY=/project/k10018/sw/envs/gspa-v3/bin/python3
export PYTHONPATH=$REPO/benchmark/v3 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
common=(--results-dir "$O/genomes" --candidates-dir "$O/candidates" --seed-reactions "$REF/seed_reactions_corrected.tsv"
        --ec-aliases "$REF/seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv"
        --unmapped-tsv "$REPO/benchmark/v3/offgrowth/data/biolog_pm1_pm4_unmapped.tsv")
mkdir -p "$O/report-price148" "$O/report-price148-swissprot"
sha256sum "$REPO/benchmark/v3/offgrowth/"*.py "$REPO/benchmark/v3/offgrowth/data/swissprot_experimental_ecs_price22.tsv" > "$O/report-price148-swissprot/code.sha256"
/usr/bin/time -v -o "$O/report-price148/time_v.txt" "$PY" -m offgrowth.runner aggregate "${common[@]}" --out-dir "$O/report-price148" > "$O/report-price148/aggregate.log" 2>&1 &
/usr/bin/time -v -o "$O/report-price148-swissprot/time_v.txt" "$PY" -m offgrowth.runner aggregate "${common[@]}" \
  --swissprot-labels "$REPO/benchmark/v3/offgrowth/data/swissprot_experimental_ecs_price22.tsv" \
  --out-dir "$O/report-price148-swissprot" > "$O/report-price148-swissprot/aggregate.log" 2>&1 &
wait -n; a=$?; wait -n; b=$?
echo "exit price148=$a swissprot=$b"; [ $a -eq 0 ] && [ $b -eq 0 ]
