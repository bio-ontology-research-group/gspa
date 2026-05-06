#!/usr/bin/env bash
#SBATCH --job-name=cult-embed
#SBATCH --partition=debug
#SBATCH --nodelist=node005
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/data/hohndor/gspa/proteomes/cultures/logs/embed-%A_%a.out
#SBATCH --error=/data/hohndor/gspa/proteomes/cultures/logs/embed-%A_%a.err

set -eo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate esmfold-v2
set -u

MANIFEST=/data/hohndor/gspa/proteomes/cultures/manifest.tsv
TAGS=($(tail -n +2 $MANIFEST | cut -f1))
tag=${TAGS[$SLURM_ARRAY_TASK_ID]}
ROOT=/data/hohndor/gspa/proteomes/cultures
OUT=$ROOT/plm
mkdir -p $OUT

# Prodigal FAA → ESM2 embeddings
# embed_proteins_plm.py expects --proteomes-dir with {tag}.faa; symlink.
mkdir -p $ROOT/faa_for_embed
ln -sf $ROOT/$tag/prodigal/$tag.faa $ROOT/faa_for_embed/$tag.faa

if [[ -f $OUT/${tag}_esm2t30.npy ]]; then
    echo "[$tag] already embedded"
    exit 0
fi

python3 /data/hohndor/gspa/benchmark-py/embed_proteins_plm.py \
    --proteomes-dir $ROOT/faa_for_embed \
    --out-dir $OUT \
    --tag $tag --batch-size 8
echo DONE; date
