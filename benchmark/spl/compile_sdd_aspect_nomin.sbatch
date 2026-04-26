#!/usr/bin/env bash
#SBATCH --job-name=gspa-spl-compile-nomin
#SBATCH --partition=debug
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -t 08:00:00
#SBATCH -o /data/hohndor/gspa-neural/logs/spl-compile-%x-%j.out
#SBATCH -e /data/hohndor/gspa-neural/logs/spl-compile-%x-%j.err

set -euo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate gapfix
export PYTHONUNBUFFERED=1

ASPECT="${1:?usage: sbatch -J spl-$ASPECT --export=ALL $0 $ASPECT}"

python -u /data/hohndor/gspa-neural/benchmark/spl/compile_sdd.py \
    --hierarchy /data/hohndor/gspa-neural/work/go_hierarchy_${ASPECT}.npz \
    --out-dir /data/hohndor/gspa-neural/work/constraints \
    --stem go_true_path_${ASPECT} \
    --no-minimize

echo DONE ${ASPECT}; date
