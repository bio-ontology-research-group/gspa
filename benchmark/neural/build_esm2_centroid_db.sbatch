#!/usr/bin/env bash
#SBATCH --job-name=gspa-esm2c-build
#SBATCH --partition=debug
#SBATCH --gres=gpu:1
#SBATCH --exclude=node001,node002,node003,node004,node007
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH -o /data/hohndor/gspa-neural/logs/esm2c-build-%j.out
#SBATCH -e /data/hohndor/gspa-neural/logs/esm2c-build-%j.err

set -euo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate gapfix
export PYTHONUNBUFFERED=1
export TORCH_HOME=/data/hohndor/gapfix/data/esm2

NEURAL=/data/hohndor/gspa-neural/benchmark/neural
WORK=/data/hohndor/gspa-neural/sprot
OUT=/data/hohndor/gspa-neural/models/esm2_centroids_t12.npz

python $NEURAL/build_esm2_centroids.py \
    --swissprot-fasta /data/hohndor/gspa/benchmark/uniprot_sprot.fasta.gz \
    --swissprot-tsv $WORK/swissprot_go_ec.tsv \
    --model esm2_t12_35M_UR50D \
    --batch-size 64 \
    --min-class-size 10 \
    --out $OUT

ls -la $OUT
echo DONE; date
