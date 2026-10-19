#!/bin/bash
#SBATCH --job-name=gspa_atb_agg
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/hohndor/allthebacteria/logs/05_aggregate_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/hohndor/allthebacteria/logs/05_aggregate_%j.err

set -euo pipefail

export GSPA_ATB_BASE="/ibex/scratch/projects/c2014/hohndor/allthebacteria"
PYTHON_BIN="/ibex/scratch/projects/c2014/rob/envs/gspa-neural/bin/python"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Aggregating completed AllTheBacteria batches on $(hostname)..."

shopt -s nullglob
WORK_DIRS=("$GSPA_ATB_BASE"/wave1/work/batch_* "$GSPA_ATB_BASE"/wave2/work/w2_*)
MANIFESTS=("$GSPA_ATB_BASE"/wave1/batches/batch_*.tsv "$GSPA_ATB_BASE"/wave2/manifests/w2_*.tsv)

"$PYTHON_BIN" "$GSPA_ATB_BASE/scripts/aggregate_atb_insights.py" \
  --work-dirs "${WORK_DIRS[@]}" \
  --batch-manifests "${MANIFESTS[@]}" \
  --out-dir "$GSPA_ATB_BASE/results/aggregated_insights"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed Stage D aggregation."
