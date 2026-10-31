#!/bin/bash
#SBATCH --job-name=gspa-price-final-audit
#SBATCH --partition=debug
#SBATCH --nodelist=node004
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=3-00:00:00
#SBATCH --no-requeue
#SBATCH --output=/data/hohndor/gspa-v3/price-finalize-v2/slurm/%j.log
set -euo pipefail
root=/data/hohndor/gspa-v3/price-finalize-v2
collection=/data/hohndor/gspa-v3/price-production-v3
jobout=$root/jobs/$SLURM_JOB_ID
export TMPDIR=/scratch/$SLURM_JOB_ID
export PYTHONPATH=$collection/code
unset JAVA_TOOL_OPTIONS JDK_JAVA_OPTIONS _JAVA_OPTIONS
mkdir -p "$jobout"
finish() {
 code=$?
 trap - EXIT
 printf '%s\n' "$code" > "$jobout/exit-code.txt"
 date --iso-8601=seconds > "$jobout/finished-at.txt"
 exit "$code"
}
trap finish EXIT
cd "$root"
sha256sum -c source.sha256 > "$jobout/source-check.log"
cp "$0" "$jobout/executed.sbatch"
python3 - <<'PY' > "$jobout/preflight.log" 2>&1
import hashlib,json
from pathlib import Path
root=Path('/data/hohndor/gspa-v3/price-production-v3')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
r=root/'collection/registration.json'
assert sha(r)=='454f466cb635b70a9bbb65cbb70d1022065568d1931ac944be95e27a754b02dd'
reg=json.loads(r.read_text())
assert sha(root/'collection/plan.json')==reg['plan_sha256']=='2407147c51dd0fe25cea0ea8193feefcac75d8cdbe46f0009d9001474339c4a4'
files=list((root/'code').glob('*.py'))
expected={Path(name) for name in reg['inputs'] if Path(name).parent==root/'code' and Path(name).suffix=='.py'}
assert expected and set(files)==expected
for p in files:assert sha(p)==reg['inputs'][str(p)],str(p)
assert not (root/'collection/collection-verification.json').exists()
print('Original registration, plan and all',len(files),'registered Python source files verified.')
PY
/usr/bin/time -v -o "$jobout/resources.txt" /data/hohndor/gspa-v3/price-production-v1/venv/bin/python \
 "$collection/code/run_price_production_collection.py" finalize --out "$collection/collection" \
 > "$jobout/finalize.log" 2>&1
cp "$collection/collection/collection-verification.json" "$jobout/"
sha256sum "$jobout/collection-verification.json" > "$jobout/receipt.sha256"
