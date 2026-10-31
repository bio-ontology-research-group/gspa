#!/bin/bash
#SBATCH --job-name=gspa-price-dpz
#SBATCH --partition=debug
#SBATCH --nodelist=node002
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=08:00:00
#SBATCH --output=/data/hohndor/gspa-v3/price-production-v1/logs/dense-cell-%j.log
set -euo pipefail
ROOT=/data/hohndor/gspa-v3/price-production-v1
export TMPDIR=/scratch/$SLURM_JOB_ID
mkdir -p "$TMPDIR"
unset JAVA_TOOL_OPTIONS JDK_JAVA_OPTIONS _JAVA_OPTIONS
export PYTHONPATH="$ROOT/code"
cd "$ROOT"
"$ROOT/venv/bin/python" - <<'PY'
import json, os
from pathlib import Path
import verified_batch as V
root=Path.cwd(); job=os.environ['SLURM_JOB_ID']
with (root/('dense-cell-launch-'+job+'.json')).open('x') as stream:
 json.dump(dict(index=66,heap='64g',allocation_mib=81920,timeout_per_arm=10800,
  inputs=V.file_hashes([root/'dense-cell.sbatch',root/'collection/registration.json'])),stream,indent=2)
PY
"$ROOT/venv/bin/python" "$ROOT/code/run_price_production_collection.py" run --out "$ROOT/collection" --index 66 --attempt "job-$SLURM_JOB_ID" --heap 64g --timeout 10800 --scratch "$TMPDIR"
"$ROOT/venv/bin/python" - <<'PY'
import json, os
from pathlib import Path
import verified_batch as V
root=Path.cwd(); job=os.environ['SLURM_JOB_ID']
attempt=root/'collection/attempts/66'/('job-'+job)
new=attempt/'production/stage1/GCF_000236665.1.integrated.tsv'
old=Path('/data/hohndor/gspa-v3/price-index-capacity-v1/dpz-GCF_000236665.1-stage1/stage1/GCF_000236665.1.integrated.tsv')
hashes=V.file_hashes([old,new]);same=hashes[str(old)]==hashes[str(new)]
with (attempt/'prior-version-comparison.json').open('x') as stream:
 json.dump(dict(byte_identical=same,hashes=hashes,interpretation='Old indexed JAR versus current frozen JAR; no labels or accuracy measure.'),stream,indent=2)
if not same: raise RuntimeError('Prior-version predictions differ; preserve both and investigate')
PY
