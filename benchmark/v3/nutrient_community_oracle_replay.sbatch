#!/bin/bash
#SBATCH --job-name=gspa-coculture-oracle
#SBATCH --partition=debug
#SBATCH --exclude=node001,node002,node005,node006,node007,node008,unimatrix
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=03:00:00
#SBATCH --signal=B:TERM@120
set -euo pipefail
GSPA_COMM_PACKAGE=${1:?immutable package directory required}
GSPA_COMM_PREVIOUS=${2:?terminal interrupted job directory required}
GSPA_COMM_PY=/data/hohndor/gspa-v3/tools/metabolic-venv/bin/python
GSPA_COMM_WORK=/scratch/$SLURM_JOB_ID/nutrient-community
GSPA_COMM_PID=
test "$(cat "$GSPA_COMM_PREVIOUS/exit-code.txt")" = 143
test ! -e "$GSPA_COMM_PREVIOUS/oracle/verification.json"
test -f "$GSPA_COMM_PREVIOUS/predictions/run.json"
mkdir -p "$GSPA_COMM_WORK"
cp -a "$GSPA_COMM_PREVIOUS/." "$GSPA_COMM_WORK/"
cd "$GSPA_COMM_WORK"
mv oracle previous-oracle
cp exit-code.txt previous-exit-code.txt
cp "$0" oracle-replay.sbatch
export TMPDIR=/scratch/$SLURM_JOB_ID
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
finish_oracle_job() {
    result_code=$?
    trap - EXIT
    printf '%s\n' "$result_code" > exit-code.txt
    mkdir -p "$GSPA_COMM_PACKAGE/jobs/oracle-replay-$SLURM_JOB_ID"
    cp -a . "$GSPA_COMM_PACKAGE/jobs/oracle-replay-$SLURM_JOB_ID/"
    exit "$result_code"
}
stop_oracle_job() {
    if [[ -n "$GSPA_COMM_PID" ]]; then
        kill -TERM "$GSPA_COMM_PID" 2>/dev/null || true
        wait "$GSPA_COMM_PID" 2>/dev/null || true
    fi
    exit 143
}
trap finish_oracle_job EXIT
trap stop_oracle_job TERM
sha256sum -c SHA256SUMS > input-verification.log
"$GSPA_COMM_PY" - "$GSPA_COMM_PREVIOUS" <<'PY'
import hashlib,json,sys
from datetime import datetime,timezone
from pathlib import Path
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
previous=Path(sys.argv[1])
files={'predictions/run.json':'predictions/run.json','predictions/predictions.jsonl':'predictions/predictions.jsonl',
       'oracle/comparisons.jsonl':'previous-oracle/comparisons.jsonl','exit-code.txt':'previous-exit-code.txt'}
identity={}
for original,copied in files.items():
    value=sha(previous/original)
    if sha(Path(copied))!=value:raise ValueError('copied terminal input differs: '+original)
    identity[copied]=value
receipt=json.loads(Path('predictions/run.json').read_text())
if identity['predictions/predictions.jsonl']!=receipt['predictions_sha256']:raise ValueError('invalid production receipt')
record=dict(schema='gspa.nutrient-community-oracle-replay/1',created_utc=datetime.now(timezone.utc).isoformat(),
    previous_job=str(previous),previous_exit_code=143,shard_index=receipt['shard_index'],shard_count=receipt['shard_count'],
    copied_inputs=identity,submission_sha256=sha(Path('oracle-replay.sbatch')),
    reason='Complete independent verification after a terminal wall-time interruption; no production rerun or outcome-based selection.',
    allowance_hours=3,production_rerun=False,oracle_prefix_reused=False)
Path('oracle-replay.json').write_text(json.dumps(record,indent=2)+'\n')
PY
"$GSPA_COMM_PY" - <<'PY' > packages.txt
from importlib.metadata import distributions
print('\n'.join(sorted(d.metadata['Name'] + '==' + d.version for d in distributions())))
PY
"$GSPA_COMM_PY" benchmark/v3/check_nutrient_community_cobra.py --plan plan --run predictions \
    --ec-model models/iML1515.xml --bs-model models/iYO844.xml --out oracle > oracle.log 2>&1 &
GSPA_COMM_PID=$!
wait "$GSPA_COMM_PID"
GSPA_COMM_PID=
