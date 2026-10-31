#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=gspa-followup-rerun
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --exclusive
#SBATCH --time=06:00:00
#SBATCH --no-requeue
#SBATCH --exclude=nid00022,nid00497,nid01033
#SBATCH --output=/scratch/project/k10018/gspa-v3/out/stage2b-followups-v1/slurm/rerun-%j.log
#
# Experiments 4 and 5: re-attribute the stored Stage 2B cells (CELLS, e.g. "0-21") with the production sidecar code
# under the synteny / candidate-limit arms; writes sparse update records. Dense cells hold payloads of up to 3.4 GB,
# so WORKERS (default 6) bounds concurrent cells per node.
set -euo pipefail
: "${REPO_ROOT:?REPO_ROOT must point at the code checkout}"
: "${CELLS:?CELLS (e.g. 0-21) required}"
OUT_ROOT=${OUT_ROOT:-/scratch/project/k10018/gspa-v3/out/stage2b-followups-v1}
COLLECTION=${COLLECTION:-/scratch/project/k10018/gspa-v3/out/price-stage2b-rc1/collection}
FBA_DATA_DIR=${FBA_DATA_DIR:-/scratch/project/k10018/gspa-v3/out/price-stage2b-rc1/fba_data}
LAYOUT_DIR=${LAYOUT_DIR:-/project/k10018/gspa-v3/price-full-matrix-v2/layout}
KEGG_MODULES=${KEGG_MODULES:-/project/k10018/gspa-v3/price-production-v1/reference/kegg_modules.tsv}
SIDECAR=${SIDECAR:-$REPO_ROOT/scripts/stoichiometric_locus_refiner.py}
PYTHON_BIN=${PYTHON_BIN:-/project/k10018/sw/envs/gspa-v3/bin/python3}
WORKERS=${WORKERS:-6}
UPDATES=${UPDATES:-$OUT_ROOT/exp45/updates}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="$REPO_ROOT/benchmark/v3"
JOB=${SLURM_JOB_ID:-local}
mkdir -p "$UPDATES" "$OUT_ROOT/slurm" "$OUT_ROOT/exp45/logs"
echo "job $JOB on $(hostname) cells=$CELLS workers=$WORKERS"
sha256sum "$SIDECAR" "$KEGG_MODULES" "$REPO_ROOT/tasks/v3-stage2b-followups-design-2026-10-31.md" > "$OUT_ROOT/exp45/logs/rerun-$JOB.inputs.sha256"
/usr/bin/time -v -o "$OUT_ROOT/exp45/logs/rerun-$JOB.time_v.txt" "$PYTHON_BIN" -m followups.stage2b_rerun \
  --collection "$COLLECTION" --sidecar "$SIDECAR" --fba-data-dir "$FBA_DATA_DIR" --layout-dir "$LAYOUT_DIR" --kegg-modules "$KEGG_MODULES" \
  --cells "$CELLS" --workers "$WORKERS" --python-bin "$PYTHON_BIN" --out-dir "$UPDATES" > "$OUT_ROOT/exp45/logs/rerun-$JOB.log" 2>&1
code=$?
echo "exit $code" > "$OUT_ROOT/exp45/logs/rerun-$JOB.exit"
exit $code
