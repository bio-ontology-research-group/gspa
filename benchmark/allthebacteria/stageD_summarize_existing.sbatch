#!/bin/bash
# Aggregation run (a): summarise the JVM outputs that Stage C already wrote
# ($S/<wave>/work/<batch>/stageC/out/<gid>.{noisy_or,gspa_full}.tsv) - read-only on the work dirs.
#   sbatch --array=1-N%30 --export=ALL,BATCH_LIST=$R/batch_lists/existing_java.txt stageD_summarize_existing.sbatch
#SBATCH --job-name=gspa_atb_sumC
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=08:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/hohndor/allthebacteria/repro_2026-10-19/logs/sumC_%A_%a.out
#SBATCH --error=/ibex/scratch/projects/c2014/hohndor/allthebacteria/repro_2026-10-19/logs/sumC_%A_%a.err
set -euo pipefail
S=/ibex/scratch/projects/c2014/hohndor/allthebacteria
R=$S/repro_2026-10-19
PY=/ibex/scratch/projects/c2014/rob/envs/gspa-neural/bin/python
: "${BATCH_LIST:?set BATCH_LIST}"
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BATCH_LIST")
WAVE=$(echo "$LINE" | cut -f1); BATCH=$(echo "$LINE" | cut -f2)
if [ "$WAVE" = wave1 ]; then MANIFEST=$S/wave1/batches/$BATCH.tsv; else MANIFEST=$S/wave2/manifests/$BATCH.tsv; fi
WORK=$S/$WAVE/work/$BATCH
$PY $R/scripts/aggregate_atb_java.py summarize \
  --work-dir "$WORK" --batch-manifest "$MANIFEST" --java-out "$WORK/stageC/out" \
  --partials "${PARTIALS:-$R/partials_stageC}" \
  --forbidden $R/ref/forbidden_by_taxon.tsv --go-dag $S/ref/go-dag.tsv \
  --complex-terms $R/ref/complex_terms.tsv --procs "${SLURM_CPUS_PER_TASK}"
