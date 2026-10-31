#!/bin/bash
#SBATCH --job-name=gspa-price-register
#SBATCH --partition=debug
#SBATCH --nodelist=node002
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/data/hohndor/gspa-v3/price-production-v1/logs/prepare-%j.log
set -euo pipefail
ROOT=/data/hohndor/gspa-v3/price-production-v1
export TMPDIR=/scratch/$SLURM_JOB_ID
mkdir -p "$TMPDIR"
unset JAVA_TOOL_OPTIONS JDK_JAVA_OPTIONS _JAVA_OPTIONS
cd "$ROOT"
mkdir reference
cp /data/hohndor/gspa/reference/{go.owl,ec2go.txt,kegg_pathways.tsv,kegg_modules.tsv} reference/
cp -a /data/hohndor/gspa/reference/reasoner-cache reference/
cp /data/hohndor/gspa-v3/intracellular-products-v1/jobs/10221/gspa-3.0.0.jar gspa-3.0.0.jar
export PYTHONPATH="$ROOT/code"
"$ROOT/venv/bin/python" - <<'PY'
import json, os, pathlib, subprocess, sys
import verified_batch as V
import ijson
root=pathlib.Path.cwd()
files=[p for p in root.iterdir() if p.is_file()]+list((root/'code').glob('*.py'))+list((root/'code/tests').glob('*.py'))
(root/'launch-inputs.json').write_text(json.dumps(dict(hashes=V.file_hashes(files),python=sys.version,
 java=subprocess.run(['java','-version'],capture_output=True,text=True,check=True).stderr,
 parser_version=ijson.__version__,parser_backend=ijson.backend,job=os.environ['SLURM_JOB_ID'],
 jar_source='/data/hohndor/gspa-v3/intracellular-products-v1/jobs/10221/gspa-3.0.0.jar'),indent=2,sort_keys=True)+'\n')
PY
"$ROOT/venv/bin/python" "$ROOT/code/run_price_production_collection.py" prepare --out "$ROOT/collection" --export /data/hohndor/gspa-v3/price-full-matrix-v2/plan.json --audit /data/hohndor/gspa-v3/price-full-collection-audit-v1/verification.json --jar "$ROOT/gspa-3.0.0.jar" --reference "$ROOT/reference"
"$ROOT/venv/bin/python" - <<'PY'
import json
from pathlib import Path
p=json.loads(Path('collection/plan.json').read_text())
sparse=[c for c in p['cells'] if c['predictor'] in ('topec','mapred','enzbert')]
assert [c['index'] for c in sparse]==list(range(66))
assert sum(c['claims'] for c in sparse)==13682870
print(json.dumps(dict(sparse_cells=len(sparse),sparse_claims=sum(c['claims'] for c in sparse),
 largest=max(p['cells'],key=lambda c:c['claims'])),sort_keys=True))
PY
