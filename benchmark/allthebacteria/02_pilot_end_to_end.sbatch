#!/bin/bash
#SBATCH --job-name=gspa_atb_pilot
#SBATCH --account=pi-hohndor
#SBATCH --partition=debug,batch
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/hohndor/allthebacteria/logs/02_pilot_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/hohndor/allthebacteria/logs/02_pilot_%j.err

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

PILOT_MANIFEST="$GSPA_ATB_BASE/wave1/batches/batch_0000_pilot.tsv"
WORK_DIR="$GSPA_ATB_BASE/wave1/work/batch_0000_pilot"
OUT_RESULTS="$GSPA_ATB_BASE/results/pilot_insights"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting End-to-End Pilot Batch on $(hostname)..."

if [ ! -f "$GSPA_ATB_BUNDLE" ]; then
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Building GSPA-Seq 3B checkpoint bundle..."
  "$PYTHON_BIN" "$GSPA_ATB_BASE/scripts/00_build_gspa_seq_checkpoint.py"
fi

if [ ! -s "$WORK_DIR/dedup/reps.faa" ]; then
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Stage A: Extracting proteomes/operons and running MMseqs2 linclust..."
  "$PYTHON_BIN" "$GSPA_ATB_BASE/scripts/stageA_extract_and_linclust.py" \
    --batch-manifest "$PILOT_MANIFEST" \
    --work-dir "$WORK_DIR" \
    --threads 16
else
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Stage A already completed ($WORK_DIR/dedup/reps.faa exists), reusing..."
fi

if [ ! -s "$WORK_DIR/stageB/rep_claims.jsonl" ]; then
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Stage B: Running GSPA-Seq 3B multi-modal inference on cluster representatives..."
  "$PYTHON_BIN" "$GSPA_ATB_BASE/scripts/stageB_gspa_seq_inference.py" \
    --work-dir "$WORK_DIR" \
    --ref-dir "$GSPA_ATB_REF" \
    --ckpt-path "$GSPA_ATB_BUNDLE" \
    --threads 16
else
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Stage B already completed ($WORK_DIR/stageB/rep_claims.jsonl exists), reusing..."
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Stage C: Running whole-genome prior integration via gspa-unified.jar..."
"$PYTHON_BIN" "$GSPA_ATB_BASE/scripts/stageC_whole_genome_integrate.py" \
  --batch-manifest "$PILOT_MANIFEST" \
  --work-dir "$WORK_DIR" \
  --ref-dir "$GSPA_ATB_REF" \
  --jar-path "$GSPA_JAR" \
  --threads 16 \
  --heap-gb 48

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Stage D: Aggregating cross-phylum biological insights..."
"$PYTHON_BIN" "$GSPA_ATB_BASE/scripts/aggregate_atb_insights.py" \
  --work-dirs "$WORK_DIR" \
  --batch-manifests "$PILOT_MANIFEST" \
  --out-dir "$OUT_RESULTS"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] End-to-End Pilot Batch completed successfully!"
