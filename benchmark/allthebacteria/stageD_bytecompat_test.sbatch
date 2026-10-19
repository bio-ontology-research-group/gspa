#!/bin/bash
# Debug-partition test: rerun Stage D on the 20 genomes of one batch that Stage C already
# integrated and compare claims.jsonl, operons_chains.tsv, noisy_or.tsv and gspa_full.tsv byte-for-byte.
#   sbatch --export=ALL,WAVE=wave1,BATCH=batch_0001 stageD_bytecompat_test.sbatch
#SBATCH --job-name=gspa_atb_D_test
#SBATCH --account=pi-hohndor
#SBATCH --partition=debug
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=02:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/hohndor/allthebacteria/repro_2026-10-19/logs/bytecompat_%j.out
set -euo pipefail
S=/ibex/scratch/projects/c2014/hohndor/allthebacteria
R=$S/repro_2026-10-19
PY=/ibex/scratch/projects/c2014/rob/envs/gspa-neural/bin/python
: "${WAVE:=wave1}"; : "${BATCH:=batch_0001}"
if [ "$WAVE" = wave1 ]; then MANIFEST=$S/wave1/batches/$BATCH.tsv; else MANIFEST=$S/wave2/manifests/$BATCH.tsv; fi
WORK=$S/$WAVE/work/$BATCH
T=$R/test_bytecompat
rm -rf "$T/$BATCH"
$PY $R/scripts/stageD_java_integrate_batch.py --work-dir "$WORK" --batch-manifest "$MANIFEST" \
  --out-root "$T" --ref-dir $S/ref --jar-path $S/bin/gspa-unified.jar \
  --java-bin /ibex/scratch/projects/c2014/rob/jdk-21/bin/java \
  --threads "$SLURM_CPUS_PER_TASK" --heap-gb 46 --chunk 20 --limit 20 --no-gzip --keep-tmp
n_same=0; n_diff=0
for f in "$T/$BATCH"/out/*.tsv; do
  b=$(basename "$f")
  if cmp -s "$f" "$WORK/stageC/out/$b"; then n_same=$((n_same+1)); else n_diff=$((n_diff+1)); echo "DIFF $b"; fi
done
c_same=0; c_diff=0
for d in "$T/$BATCH"/tmp/chunk_0000/*/; do
  g=$(basename "$d")
  for k in claims.jsonl operons_chains.tsv; do
    if cmp -s "$d/$k" "$WORK/genomes/$g/$k"; then c_same=$((c_same+1)); else c_diff=$((c_diff+1)); echo "DIFF $g/$k"; fi
  done
done
echo "RESULT outputs identical=$n_same differ=$n_diff ; inputs identical=$c_same differ=$c_diff"
rm -rf "$T/$BATCH/tmp"
[ "$n_diff" -eq 0 ] && [ "$c_diff" -eq 0 ] && [ "$n_same" -gt 0 ]
