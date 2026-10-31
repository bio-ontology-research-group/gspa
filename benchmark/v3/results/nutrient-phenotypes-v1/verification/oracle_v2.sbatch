#!/bin/bash
#SBATCH --job-name=gspa-oracle-v2
#SBATCH --partition=debug
#SBATCH --nodelist=node004
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:40:00
set -euo pipefail
GSPA_ORACLE_PACKAGE=/data/hohndor/gspa-v3/nutrient-phenotypes-v1
GSPA_ORACLE_WORK=/scratch/$SLURM_JOB_ID/oracle
mkdir -p "$GSPA_ORACLE_WORK"
cp -a "$GSPA_ORACLE_PACKAGE/jobs/10144/." "$GSPA_ORACLE_WORK/"
cp "$GSPA_ORACLE_PACKAGE/check_nutrient_cobra_v2.py" "$GSPA_ORACLE_WORK/benchmark/v3/check_nutrient_cobra.py"
cd "$GSPA_ORACLE_WORK"
export TMPDIR=/scratch/$SLURM_JOB_ID OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
finish_oracle() {
    result_code=$?
    trap - EXIT
    GSPA_ORACLE_OUT="$GSPA_ORACLE_PACKAGE/oracle-v2/$SLURM_JOB_ID"
    mkdir -p "$GSPA_ORACLE_OUT"
    printf '%s\n' "$result_code" > "$GSPA_ORACLE_OUT/exit-code.txt"
    cp benchmark/v3/check_nutrient_cobra.py "$GSPA_ORACLE_OUT/"
    cp oracle-v2.json "$GSPA_ORACLE_OUT/" || true
    exit "$result_code"
}
trap finish_oracle EXIT
/data/hohndor/gspa-v3/tools/metabolic-venv/bin/python benchmark/v3/check_nutrient_cobra.py \
    --plan plan --run predictions --ec-model models/iML1515.xml --bs-model models/iYO844.xml --out oracle-v2.json
