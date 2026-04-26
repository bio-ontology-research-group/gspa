#!/usr/bin/env bash
#SBATCH --job-name=gspa-spl-compile-long
#SBATCH --partition=debug
#SBATCH -c 4
#SBATCH --mem=128G
#SBATCH -t 7-00:00:00
#SBATCH -o /data/hohndor/gspa-neural/logs/spl-compile-%x-%j.out
#SBATCH -e /data/hohndor/gspa-neural/logs/spl-compile-%x-%j.err

# Long-running SDD compile (up to 7 days). Uses auto_gc_and_minimize=True
# since it produced the most compact SDDs on CC. If it finishes for MF or
# BP we gain full-vocab SPL coverage; otherwise we keep CC-only SPL as the
# deliverable.

set -euo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate gapfix
export PYTHONUNBUFFERED=1

ASPECT="${1:?usage: sbatch -J spl-long-<ASPECT> --export=ALL $0 <ASPECT>}"

python -u /data/hohndor/gspa-neural/benchmark/spl/compile_sdd.py \
    --hierarchy /data/hohndor/gspa-neural/work/go_hierarchy_${ASPECT}.npz \
    --out-dir /data/hohndor/gspa-neural/work/constraints_long \
    --stem go_true_path_${ASPECT}

echo DONE ${ASPECT}; date
