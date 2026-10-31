#!/bin/bash
#SBATCH --job-name=gspa-lp-check
#SBATCH --partition=debug
#SBATCH --exclude=node001,node005,node008,unimatrix
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:05:00
set -euo pipefail
GSPA_DIAGNOSTIC_ROOT=/data/hohndor/gspa-v3/nutrient-phenotypes-v1
GSPA_DIAGNOSTIC_WORK=/scratch/$SLURM_JOB_ID/diagnosis
mkdir -p "$GSPA_DIAGNOSTIC_WORK"
cp -a "$GSPA_DIAGNOSTIC_ROOT/plan" "$GSPA_DIAGNOSTIC_ROOT/models" "$GSPA_DIAGNOSTIC_WORK/"
export TMPDIR=/scratch/$SLURM_JOB_ID OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
/data/hohndor/gspa-v3/tools/metabolic-venv/bin/python "$GSPA_DIAGNOSTIC_ROOT/diagnose_glpk.py" \
    "$GSPA_DIAGNOSTIC_WORK" "$GSPA_DIAGNOSTIC_ROOT/jobs/10144" > "$GSPA_DIAGNOSTIC_ROOT/glpk-diagnosis.json"
