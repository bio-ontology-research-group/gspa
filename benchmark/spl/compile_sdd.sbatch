#!/usr/bin/env bash
#SBATCH --job-name=gspa-spl-compile
#SBATCH --partition=debug
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -o /data/hohndor/gspa-neural/logs/spl-compile-%j.out
#SBATCH -e /data/hohndor/gspa-neural/logs/spl-compile-%j.err

set -euo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate gapfix

cd /data/hohndor/gspa-neural/work

python /data/hohndor/gspa-neural/benchmark/spl/compile_sdd.py \
    --hierarchy /data/hohndor/gspa-neural/work/go_hierarchy_5707.npz \
    --out-dir /data/hohndor/gspa-neural/work/constraints \
    --stem go_true_path_5707

echo DONE; date
ls -la /data/hohndor/gspa-neural/work/constraints/
