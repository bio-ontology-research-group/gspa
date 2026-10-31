#!/bin/bash
#SBATCH --job-name=gspa-price-capacity-dpz
#SBATCH --partition=debug
#SBATCH --nodelist=node003
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=02:30:00
#SBATCH --output=/data/hohndor/gspa-v3/price-capacity-v2/run-%j.log
set -euo pipefail
export TMPDIR="${SLURM_TMPDIR:-/scratch/$SLURM_JOB_ID}"
export JAVA_TOOL_OPTIONS="-Djava.io.tmpdir=$TMPDIR"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p "$TMPDIR"
cd /data/hohndor/gspa-v3/price-capacity-v2
python3 - <<'CHECK'
import hashlib,json
for path,expected in json.load(open('inputs.json')).items():
    h=hashlib.sha256()
    with open(path,'rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''): h.update(chunk)
    if h.hexdigest()!=expected: raise SystemExit('frozen capacity input changed: '+path)
CHECK
/data/hohndor/gspa-v3/tools/metabolic-venv/bin/python code/run_price_matrix_capacity.py --part /data/hohndor/gspa-v3/price-full-matrix-v2/parts/dpz/GCF_000236665.1 --layout /data/hohndor/gspa-v3/price-full-matrix-v2/layout --jar /data/hohndor/gspa-v3/price-capacity-v2/gspa-3.0.0.jar --reference /data/hohndor/gspa/reference --out /data/hohndor/gspa-v3/price-capacity-v2/dpz-GCF_000236665.1-stage0 --java-heap 64g --arms stage0
