#!/bin/bash
#SBATCH --job-name=gspa-final-lp-check
#SBATCH --partition=debug
#SBATCH --nodelist=node004
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:05:00
set -euo pipefail
GSPA_DIAGNOSTIC_ROOT=/data/hohndor/gspa-v3/nutrient-phenotypes-v1
GSPA_DIAGNOSTIC_WORK=/scratch/$SLURM_JOB_ID/diagnosis
mkdir -p "$GSPA_DIAGNOSTIC_WORK"
cp -a "$GSPA_DIAGNOSTIC_ROOT/plan" "$GSPA_DIAGNOSTIC_ROOT/models" "$GSPA_DIAGNOSTIC_ROOT/scripts" "$GSPA_DIAGNOSTIC_WORK/"
export TMPDIR=/scratch/$SLURM_JOB_ID OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
/data/hohndor/gspa-v3/tools/metabolic-venv/bin/python "$GSPA_DIAGNOSTIC_ROOT/diagnose_glpk_final.py" "$GSPA_DIAGNOSTIC_WORK" "$GSPA_DIAGNOSTIC_ROOT/oracle-v2/10148" > "$GSPA_DIAGNOSTIC_ROOT/glpk-final-diagnosis.json"
/data/hohndor/gspa-v3/tools/metabolic-venv/bin/python "$GSPA_DIAGNOSTIC_ROOT/diagnose_finite.py" "$GSPA_DIAGNOSTIC_WORK" > "$GSPA_DIAGNOSTIC_ROOT/highs-finite-diagnosis.json"
