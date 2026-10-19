#!/bin/bash
#SBATCH --job-name=gspa_atb_wave1
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch,largemem
#SBATCH --array=1-42%42
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/hohndor/allthebacteria/logs/03_wave1_%A_%a.out
#SBATCH --error=/ibex/scratch/projects/c2014/hohndor/allthebacteria/logs/03_wave1_%A_%a.err

set -euo pipefail

export GSPA_ATB_BASE="/ibex/scratch/projects/c2014/hohndor/allthebacteria"
export GSPA_ATB_REF="$GSPA_ATB_BASE/ref"
export GSPA_ATB_CKPT="$GSPA_ATB_BASE/checkpoints"
export GSPA_ATB_BUNDLE="$GSPA_ATB_CKPT/gspa_seq_3b_bundle.pt"
export GSPA_JAR="$GSPA_ATB_BASE/bin/gspa-unified.jar"
export JAVA_BIN="/ibex/scratch/projects/c2014/rob/jdk-21/bin/java"
export PYRODIGAL_BIN="/ibex/scratch/projects/c2014/rob/miniconda3/envs/pyrodigal/bin/pyrodigal"
export MMSEQS_BIN="/ibex/scratch/projects/c2014/rob/miniconda3/envs/mmseqs2/bin/mmseqs"
export DIAMOND_BIN="/ibex/scratch/projects/c2014/rob/miniconda3/envs/bakta/bin/diamond"
export HMMSEARCH_BIN="/ibex/scratch/projects/c2014/rob/miniconda3/envs/bakta/bin/hmmsearch"
PYTHON_BIN="/ibex/scratch/projects/c2014/rob/envs/gspa-neural/bin/python"

BATCH_ID=$(printf "batch_%04d" "${SLURM_ARRAY_TASK_ID}")
BATCH_MANIFEST="$GSPA_ATB_BASE/wave1/batches/${BATCH_ID}.tsv"
WORK_DIR="$GSPA_ATB_BASE/wave1/work/${BATCH_ID}"

if [ ! -f "$BATCH_MANIFEST" ]; then
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Batch manifest $BATCH_MANIFEST not present (array index beyond Wave 1 batch count); exiting cleanly."
  exit 0
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting Wave 1 species-rep batch ${BATCH_ID} on $(hostname)..."

"$PYTHON_BIN" "$GSPA_ATB_BASE/scripts/stageA_extract_and_linclust.py" \
  --batch-manifest "$BATCH_MANIFEST" \
  --work-dir "$WORK_DIR" \
  --threads 16

"$PYTHON_BIN" "$GSPA_ATB_BASE/scripts/stageB_gspa_seq_inference.py" \
  --work-dir "$WORK_DIR" \
  --ref-dir "$GSPA_ATB_REF" \
  --ckpt-path "$GSPA_ATB_BUNDLE" \
  --threads 16

"$PYTHON_BIN" "$GSPA_ATB_BASE/scripts/stageC_whole_genome_integrate.py" \
  --batch-manifest "$BATCH_MANIFEST" \
  --work-dir "$WORK_DIR" \
  --ref-dir "$GSPA_ATB_REF" \
  --jar-path "$GSPA_JAR" \
  --threads 16 \
  --heap-gb 48

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed Wave 1 batch ${BATCH_ID}."
