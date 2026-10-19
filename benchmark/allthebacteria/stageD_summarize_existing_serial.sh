#!/bin/bash
# Serial variant of stageD_summarize_existing.sbatch: one Slurm job loops over a batch list and
# summarises every batch that has no partial yet (used because per-task scheduling under the
# per-user CPU cap was the bottleneck, not compute).
#   sbatch -p batch -c 8 --mem=32G -t 06:00:00 stageD_summarize_existing_serial.sh LIST PARTIALS
set -euo pipefail
S=/ibex/scratch/projects/c2014/hohndor/allthebacteria
R=$S/repro_2026-10-19
PY=/ibex/scratch/projects/c2014/rob/envs/gspa-neural/bin/python
LIST=$1; PARTIALS=$2
mkdir -p "$PARTIALS"
while IFS=$'\t' read -r WAVE BATCH; do
  [ -f "$PARTIALS/$BATCH.genomes.tsv" ] && continue
  if [ "$WAVE" = wave1 ]; then MANIFEST=$S/wave1/batches/$BATCH.tsv; else MANIFEST=$S/wave2/manifests/$BATCH.tsv; fi
  WORK=$S/$WAVE/work/$BATCH
  $PY $R/scripts/aggregate_atb_java.py summarize --work-dir "$WORK" --batch-manifest "$MANIFEST" \
    --java-out "$WORK/stageC/out" --partials "$PARTIALS" \
    --forbidden $R/ref/forbidden_by_taxon.tsv --go-dag $S/ref/go-dag.tsv \
    --complex-terms $R/ref/complex_terms.tsv --procs "${SLURM_CPUS_PER_TASK:-4}"
done < "$LIST"
