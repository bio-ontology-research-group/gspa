#!/bin/bash
#SBATCH --job-name=gspa-nutrients
#SBATCH --partition=debug
#SBATCH --exclude=node001,node005,node007,node008,unimatrix
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
set -euo pipefail
GSPA_NUTRIENT_PACKAGE=${1:?immutable package directory required}
GSPA_NUTRIENT_PY=/data/hohndor/gspa-v3/tools/metabolic-venv/bin/python
GSPA_NUTRIENT_WORK=/scratch/$SLURM_JOB_ID/nutrient-phenotypes
mkdir -p "$GSPA_NUTRIENT_WORK"
cp -a "$GSPA_NUTRIENT_PACKAGE/benchmark" "$GSPA_NUTRIENT_PACKAGE/scripts" \
      "$GSPA_NUTRIENT_PACKAGE/plan" "$GSPA_NUTRIENT_PACKAGE/models" "$GSPA_NUTRIENT_WORK/"
cd "$GSPA_NUTRIENT_WORK"
export TMPDIR=/scratch/$SLURM_JOB_ID
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
finish_nutrient_job() {
    result_code=$?
    trap - EXIT
    printf '%s\n' "$result_code" > exit-code.txt
    mkdir -p "$GSPA_NUTRIENT_PACKAGE/jobs/$SLURM_JOB_ID"
    cp -a . "$GSPA_NUTRIENT_PACKAGE/jobs/$SLURM_JOB_ID/"
    exit "$result_code"
}
trap finish_nutrient_job EXIT
"$GSPA_NUTRIENT_PY" - <<'PY' > packages.txt
from importlib.metadata import distributions
print('\n'.join(sorted(f"{d.metadata['Name']}=={d.version}" for d in distributions())))
PY
"$GSPA_NUTRIENT_PY" benchmark/v3/run_nutrient_phenotypes.py run \
    --plan plan --out predictions --ec-model models/iML1515.xml --bs-model models/iYO844.xml
"$GSPA_NUTRIENT_PY" benchmark/v3/check_nutrient_cobra.py \
    --plan plan --run predictions --ec-model models/iML1515.xml --bs-model models/iYO844.xml --out oracle.json
