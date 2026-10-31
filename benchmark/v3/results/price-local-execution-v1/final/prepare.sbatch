#!/bin/bash
#SBATCH --job-name=gspa-price-local-pilot
#SBATCH --partition=debug
#SBATCH --nodelist=node004
#SBATCH --cpus-per-task=2
#SBATCH --mem=36G
#SBATCH --time=02:00:00
#SBATCH --output=/data/hohndor/gspa-v3/price-production-v3/logs/prepare-%j.log
set -euo pipefail
ROOT=/data/hohndor/gspa-v3/price-production-v3
export TMPDIR=/scratch/$SLURM_JOB_ID
mkdir -p "$TMPDIR"
unset JAVA_TOOL_OPTIONS JDK_JAVA_OPTIONS _JAVA_OPTIONS
PYTHON=/data/hohndor/gspa-v3/price-production-v1/venv/bin/python
export PYTHONPATH="$ROOT/test-source/benchmark/v3:$ROOT/test-source/benchmark/v3/tests"
cd "$ROOT"
finish() {
 code=$?
 printf '%s\n' "$code" > "$ROOT/logs/prepare-$SLURM_JOB_ID.exit-code.txt"
 exit "$code"
}
trap finish EXIT
"$PYTHON" - <<'CHECK'
import hashlib,json,os,platform,sys
from pathlib import Path
root=Path.cwd(); hashes={}
for manifest,base in [('code-hashes.json',root/'code'),('test-source-hashes.json',root/'test-source')]:
 for relative,expected in json.loads((root/manifest).read_text()).items():
  path=base/relative;actual=hashlib.sha256(path.read_bytes()).hexdigest()
  if actual!=expected:raise ValueError(f'source package changed: {path}')
  hashes[str(path)]=actual
for name in ('prepare.sbatch','check_pilot.py','code.tar.gz','test-source.tar.gz','code-hashes.json','test-source-hashes.json'):
 hashes[str(root/name)]=hashlib.sha256((root/name).read_bytes()).hexdigest()
path=root/'logs'/('launch-'+os.environ['SLURM_JOB_ID']+'.json')
with path.open('x') as stream:
 json.dump(dict(hashes=hashes,python=sys.version,node=platform.node(),job=os.environ['SLURM_JOB_ID']),stream,indent=2)
 stream.write('\n')
CHECK
"$PYTHON" -m unittest -v test_price_local test_price_capacity test_price_production_collection test_verified_batch test_streaming_ledger > "$ROOT/logs/tests-$SLURM_JOB_ID.log" 2>&1
export PYTHONPATH="$ROOT/code"
cp /data/hohndor/gspa-v3/compact-audit-v1/jobs/10275/gspa-3.0.0.jar "$ROOT/gspa-3.0.0.jar"
"$PYTHON" "$ROOT/code/run_price_production_collection.py" prepare --out "$ROOT/collection" --export /data/hohndor/gspa-v3/price-full-matrix-v2/plan.json --audit /data/hohndor/gspa-v3/price-full-collection-audit-v1/verification.json --jar "$ROOT/gspa-3.0.0.jar" --reference /data/hohndor/gspa-v3/price-production-v1/reference
"$PYTHON" "$ROOT/code/run_price_local.py" --out "$ROOT/collection" --index 0 --attempt "job-$SLURM_JOB_ID" --heap 24g --timeout 7200 --scratch "$TMPDIR" --container /data/hohndor/gspa-v3/observed-gap-context-modes-v1-relocation/container-root
"$PYTHON" "$ROOT/check_pilot.py" "$SLURM_JOB_ID"
