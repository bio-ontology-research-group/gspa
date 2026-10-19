#!/bin/bash
# Stage D: GSPA JVM integrator on ALL genomes of completed AllTheBacteria batches.
# One array task per batch; the batch list is a text file with lines "<wave>\t<batch>".
#
#   sbatch --array=1-43%20 --export=ALL,BATCH_LIST=$R/batch_lists/wave1.txt stageD_java_integrate_all.sbatch
#
# Per task: (1) stageD_java_integrate_batch.py regenerates claims.jsonl / operons_chains.tsv
# exactly as Stage C does, runs `integrate --batch-manifest` (shared ontology load, two rows per
# genome: noisy_or + full) in chunks, writes gz outputs to $R/<batch>/out/;
# (2) aggregate_atb_java.py summarize writes the per-genome partials to $R/partials_stageD/.
#SBATCH --job-name=gspa_atb_stageD
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=16:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/hohndor/allthebacteria/repro_2026-10-19/logs/stageD_%A_%a.out
#SBATCH --error=/ibex/scratch/projects/c2014/hohndor/allthebacteria/repro_2026-10-19/logs/stageD_%A_%a.err
set -euo pipefail

S=/ibex/scratch/projects/c2014/hohndor/allthebacteria
R=$S/repro_2026-10-19
PY=/ibex/scratch/projects/c2014/rob/envs/gspa-neural/bin/python
export JAVA_BIN=/ibex/scratch/projects/c2014/rob/jdk-21/bin/java
: "${BATCH_LIST:?set BATCH_LIST}"
: "${CHUNK:=96}"
: "${LIMIT:=0}"
: "${NO_GZIP:=0}"
: "${OUT_ROOT:=$R}"

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BATCH_LIST")
WAVE=$(echo "$LINE" | cut -f1)
BATCH=$(echo "$LINE" | cut -f2)
if [ "$WAVE" = wave1 ]; then MANIFEST=$S/wave1/batches/$BATCH.tsv; else MANIFEST=$S/wave2/manifests/$BATCH.tsv; fi
WORK=$S/$WAVE/work/$BATCH
echo "[$(date '+%F %T')] task ${SLURM_ARRAY_TASK_ID}: $WAVE $BATCH on $(hostname), cpus=${SLURM_CPUS_PER_TASK}"
test -f "$WORK/stageC/batch_stats.tsv" || { echo "batch $BATCH has not finished Stage C; skipping"; exit 0; }

EXTRA=()
[ "$LIMIT" != 0 ] && EXTRA+=(--limit "$LIMIT")
[ "$NO_GZIP" != 0 ] && EXTRA+=(--no-gzip)
[ "${ADOPT:-1}" != 0 ] && EXTRA+=(--adopt-existing)

$PY $R/scripts/stageD_java_integrate_batch.py \
  --work-dir "$WORK" --batch-manifest "$MANIFEST" --out-root "$OUT_ROOT" \
  --ref-dir $S/ref --jar-path $S/bin/gspa-unified.jar \
  --threads "${SLURM_CPUS_PER_TASK}" --heap-gb 46 --chunk "$CHUNK" "${EXTRA[@]}"

if [ "$LIMIT" = 0 ]; then
  $PY $R/scripts/aggregate_atb_java.py summarize \
    --work-dir "$WORK" --batch-manifest "$MANIFEST" --java-out "$OUT_ROOT/$BATCH/out" \
    --partials "$OUT_ROOT/partials_stageD" \
    --forbidden $R/ref/forbidden_by_taxon.tsv --go-dag $S/ref/go-dag.tsv \
    --complex-terms $R/ref/complex_terms.tsv --procs "${SLURM_CPUS_PER_TASK}"
fi
echo "[$(date '+%F %T')] done $BATCH"
