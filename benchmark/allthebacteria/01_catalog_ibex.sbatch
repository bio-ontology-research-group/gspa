#!/bin/bash
#SBATCH --job-name=gspa_atb_catalog
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/hohndor/allthebacteria/logs/01_catalog_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/hohndor/allthebacteria/logs/01_catalog_%j.err

set -euo pipefail

export GSPA_ATB_BASE="/ibex/scratch/projects/c2014/hohndor/allthebacteria"
PYTHON_BIN="/ibex/scratch/projects/c2014/rob/envs/gspa-neural/bin/python"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting AllTheBacteria catalog download and Wave 1 / Wave 2 batch generation on $(hostname)..."

"$PYTHON_BIN" "$GSPA_ATB_BASE/scripts/01_download_and_catalog.py" \
  --base-dir "$GSPA_ATB_BASE" \
  --reps-per-species 3 \
  --wave1-batch-size 500

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed AllTheBacteria catalog generation."
