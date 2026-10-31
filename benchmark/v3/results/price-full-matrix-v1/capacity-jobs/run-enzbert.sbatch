#!/bin/bash
#SBATCH --job-name=gspa-price-capacity-enzbert
#SBATCH --partition=debug
#SBATCH --nodelist=node003
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=02:30:00
#SBATCH --output=/data/hohndor/gspa-v3/price-capacity-v1/run-%j.log
set -euo pipefail
export TMPDIR="${SLURM_TMPDIR:-/scratch/$SLURM_JOB_ID}"
export JAVA_TOOL_OPTIONS="-Djava.io.tmpdir=$TMPDIR"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p "$TMPDIR"
/data/hohndor/gspa-v3/tools/metabolic-venv/bin/python /data/hohndor/gspa-v3/price-capacity-v1/code/run_price_matrix_capacity.py --part /data/hohndor/gspa-v3/price-full-matrix-v2/parts/enzbert/GCF_000236665.1 --layout /data/hohndor/gspa-v3/price-full-matrix-v2/layout --jar /data/hohndor/gspa-v3/price-capacity-v1/gspa-3.0.0.jar --reference /data/hohndor/gspa/reference --out /data/hohndor/gspa-v3/price-capacity-v1/enzbert-GCF_000236665.1
