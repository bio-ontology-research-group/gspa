#!/bin/bash
#SBATCH --job-name=gspa_atb_wave2
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch,largemem
#SBATCH --array=1-1027%16
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/hohndor/allthebacteria/logs/04_wave2_%A_%a.out
#SBATCH --error=/ibex/scratch/projects/c2014/hohndor/allthebacteria/logs/04_wave2_%A_%a.err

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

BATCH_ID=$(printf "w2_%04d" "${SLURM_ARRAY_TASK_ID}")
BATCH_MANIFEST="$GSPA_ATB_BASE/wave2/manifests/${BATCH_ID}.tsv"
WORK_DIR="$GSPA_ATB_BASE/wave2/work/${BATCH_ID}"

if [ ! -f "$BATCH_MANIFEST" ]; then
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Wave 2 manifest $BATCH_MANIFEST not present; exiting cleanly."
  exit 0
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting Wave 2 full 2.4M archive batch ${BATCH_ID} on $(hostname)..."

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

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed Wave 2 archive batch ${BATCH_ID}."
