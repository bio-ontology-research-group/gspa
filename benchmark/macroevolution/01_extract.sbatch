#!/bin/bash
# One task per Stage D batch: genome-level function presence for both integrator configurations.
#   sbatch --array=1-N%40 --export=ALL,BATCH_LIST=$M/lists/batches.txt 01_extract.sbatch
# BATCH_LIST lines: <wave>\t<batch>
#SBATCH --job-name=gspa_macroevo_extract
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=/encrypted/e3008/gspa/macroevo/logs/extract_%A_%a.out
set -euo pipefail
S=/ibex/scratch/projects/c2014/hohndor/allthebacteria
R=/encrypted/e3008/allthebacteria/repro_2026-10-19
M=/encrypted/e3008/gspa/macroevo
PY=/ibex/scratch/projects/c2014/rob/envs/gspa-neural/bin/python
: "${BATCH_LIST:?set BATCH_LIST}"
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BATCH_LIST")
WAVE=$(echo "$LINE" | cut -f1); BATCH=$(echo "$LINE" | cut -f2)
if [ "$WAVE" = wave1 ]; then MANIFEST=$S/wave1/batches/$BATCH.tsv; else MANIFEST=$S/wave2/manifests/$BATCH.tsv; fi
test -f $M/presence/$BATCH.npz && { echo "$BATCH already extracted"; exit 0; }
$PY $M/code/extract_presence.py --batch $BATCH --out-root $R --manifest $MANIFEST \
  --go-dag $S/ref/go-dag.tsv --obo $S/ref/go.obo --essential-dir $M/ref/essential-functions \
  --out $M/presence --procs ${SLURM_CPUS_PER_TASK}
