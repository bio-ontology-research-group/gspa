#!/bin/bash
#SBATCH --job-name=gspa-cocultures
#SBATCH --partition=debug
#SBATCH --exclude=node001,node002,node005,node006,node007,node008,unimatrix
#SBATCH --array=0-5%2
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --signal=B:TERM@120
set -euo pipefail
GSPA_COMM_PACKAGE=${1:?immutable package directory required}
GSPA_COMM_PY=/data/hohndor/gspa-v3/tools/metabolic-venv/bin/python
GSPA_COMM_WORK=/scratch/$SLURM_JOB_ID/nutrient-community
GSPA_COMM_PID=
mkdir -p "$GSPA_COMM_WORK"
cp -a "$GSPA_COMM_PACKAGE/benchmark" "$GSPA_COMM_PACKAGE/scripts" \
      "$GSPA_COMM_PACKAGE/plan" "$GSPA_COMM_PACKAGE/models" "$GSPA_COMM_PACKAGE/SHA256SUMS" "$GSPA_COMM_WORK/"
cd "$GSPA_COMM_WORK"
export TMPDIR=/scratch/$SLURM_JOB_ID
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
finish_community_job() {
    result_code=$?
    trap - EXIT
    printf '%s\n' "$result_code" > exit-code.txt
    mkdir -p "$GSPA_COMM_PACKAGE/jobs/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
    cp -a . "$GSPA_COMM_PACKAGE/jobs/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}/"
    exit "$result_code"
}
stop_community_job() {
    if [[ -n "$GSPA_COMM_PID" ]]; then
        kill -TERM "$GSPA_COMM_PID" 2>/dev/null || true
        wait "$GSPA_COMM_PID" 2>/dev/null || true
    fi
    exit 143
}
trap finish_community_job EXIT
trap stop_community_job TERM
sha256sum -c SHA256SUMS > input-verification.log
"$GSPA_COMM_PY" - <<'PY' > packages.txt
from importlib.metadata import distributions
print('\n'.join(sorted(d.metadata['Name'] + '==' + d.version for d in distributions())))
PY
"$GSPA_COMM_PY" benchmark/v3/run_nutrient_community.py run --plan plan --out predictions \
    --ec-model models/iML1515.xml --bs-model models/iYO844.xml \
    --shard-index "$SLURM_ARRAY_TASK_ID" --shard-count 6 > production.log 2>&1 &
GSPA_COMM_PID=$!
wait "$GSPA_COMM_PID"
GSPA_COMM_PID=
"$GSPA_COMM_PY" benchmark/v3/check_nutrient_community_cobra.py --plan plan --run predictions \
    --ec-model models/iML1515.xml --bs-model models/iYO844.xml --out oracle > oracle.log 2>&1 &
GSPA_COMM_PID=$!
wait "$GSPA_COMM_PID"
GSPA_COMM_PID=
