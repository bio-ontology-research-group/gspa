#!/bin/bash
# Same as 01_extract.sbatch but one job loops over the batch list (used when the 2,000-job QOS cap leaves
# no room for an array): lines i with i % STRIDE == OFFSET, so STRIDE jobs can run side by side.
#   for o in 0 1 2; do sbatch --export=ALL,BATCH_LIST=$M/lists/batches.txt,STRIDE=3,OFFSET=$o 01_extract_loop.sbatch; done
#SBATCH --job-name=gspa_macroevo_extract
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=32
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=/encrypted/e3008/gspa/macroevo/logs/extract_loop_%j.out
set -euo pipefail
S=/ibex/scratch/projects/c2014/hohndor/allthebacteria
R=/encrypted/e3008/allthebacteria/repro_2026-10-19
M=/encrypted/e3008/gspa/macroevo
PY=/ibex/scratch/projects/c2014/rob/envs/gspa-neural/bin/python
: "${BATCH_LIST:?set BATCH_LIST}"; : "${STRIDE:=1}"; : "${OFFSET:=0}"
i=0
while IFS=$'\t' read -r WAVE BATCH; do
  i=$((i+1))
  [ $(( (i-1) % STRIDE )) -eq "$OFFSET" ] || continue
  if [ "$WAVE" = wave1 ]; then MANIFEST=$S/wave1/batches/$BATCH.tsv; else MANIFEST=$S/wave2/manifests/$BATCH.tsv; fi
  if [ -f $M/presence/$BATCH.npz ]; then echo "$BATCH already extracted"; continue; fi
  $PY $M/code/extract_presence.py --batch $BATCH --out-root $R --manifest $MANIFEST \
    --go-dag $S/ref/go-dag.tsv --obo $S/ref/go.obo --essential-dir $M/ref/essential-functions \
    --out $M/presence --procs ${SLURM_CPUS_PER_TASK} || echo "FAILED $BATCH"
done < "$BATCH_LIST"
echo "[$(date '+%F %T')] loop done (stride $STRIDE offset $OFFSET)"
