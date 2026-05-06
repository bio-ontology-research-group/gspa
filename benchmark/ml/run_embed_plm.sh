#!/usr/bin/env bash
#SBATCH --job-name=plm-embed
#SBATCH --partition=debug
#SBATCH --nodelist=node005
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=/data/hohndor/gspa/proteomes/bench_gtdb30/plm-%A_%a.out
#SBATCH --error=/data/hohndor/gspa/proteomes/bench_gtdb30/plm-%A_%a.err
# Usage: sbatch --array=0-29 run_embed_plm.sh

set -eo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate esmfold-v2
set -u

ROOT=/data/hohndor/gspa/proteomes/bench_gtdb30
TAGS=($(tail -n +2 $ROOT/panel_manifest.tsv | cut -f1))
tag=${TAGS[$SLURM_ARRAY_TASK_ID]}

OUT=$ROOT/plm
mkdir -p $OUT

if [[ -f $OUT/${tag}_esm2t30.npy ]]; then
    echo "[$tag] already embedded, skipping"
    exit 0
fi

echo "=== $tag: ESM-2 t30 embedding ==="
date
python3 /data/hohndor/gspa/benchmark-py/embed_proteins_plm.py \
    --proteomes-dir $ROOT/proteomes \
    --out-dir $OUT \
    --tag $tag \
    --batch-size 8
echo "DONE"
date
