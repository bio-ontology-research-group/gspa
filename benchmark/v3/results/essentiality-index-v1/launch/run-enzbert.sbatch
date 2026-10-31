#!/bin/bash
#SBATCH --job-name=gspa-index-enzbert
#SBATCH --partition=debug
#SBATCH --nodelist=node003
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=04:30:00
#SBATCH --output=/data/hohndor/gspa-v3/price-index-capacity-v1/enzbert-%j.log
set -euo pipefail
export TMPDIR="${SLURM_TMPDIR:-/scratch/$SLURM_JOB_ID}"
export JAVA_TOOL_OPTIONS="-Djava.io.tmpdir=$TMPDIR"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p "$TMPDIR"
cd /data/hohndor/gspa-v3/price-index-capacity-v1
python3 - <<'CHECK'
import hashlib,json
for path,expected in json.load(open('inputs.json')).items():
    h=hashlib.sha256()
    with open(path,'rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''): h.update(chunk)
    if h.hexdigest()!=expected: raise SystemExit('frozen capacity input changed: '+path)
CHECK
/data/hohndor/gspa-v3/tools/metabolic-venv/bin/python code/run_price_matrix_capacity.py --part /data/hohndor/gspa-v3/price-full-matrix-v2/parts/enzbert/GCF_000236665.1 --layout /data/hohndor/gspa-v3/price-full-matrix-v2/layout --jar /data/hohndor/gspa-v3/price-index-capacity-v1/gspa-3.0.0.jar --reference /data/hohndor/gspa/reference --out /data/hohndor/gspa-v3/price-index-capacity-v1/enzbert-GCF_000236665.1-stage1 --java-heap 8g --arms stage1 --timeout 14400
