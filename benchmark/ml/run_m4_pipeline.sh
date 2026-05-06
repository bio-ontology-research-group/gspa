#!/usr/bin/env bash
#SBATCH --job-name=m4-pipeline
#SBATCH --partition=debug
#SBATCH --exclude=node003
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -o /data/hohndor/gspa/proteomes/bench_gtdb30/m4-pipeline-%j.out
#SBATCH -e /data/hohndor/gspa/proteomes/bench_gtdb30/m4-pipeline-%j.err

set -eo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics
set -u

ROOT=/data/hohndor/gspa/proteomes/bench_gtdb30
BM=/data/hohndor/gspa/benchmark-py

cd $ROOT

echo "=== Stage 1: build PLM centroids ==="
date
python3 $BM/build_plm_centroids.py \
    --manifest panel_manifest.tsv \
    --root $ROOT \
    --plm-dir $ROOT/plm \
    --ec2go /data/hohndor/gspa/reference/ec2go.txt \
    --exclude-tag mg1655 \
    --tau 0.3 \
    --out-dir $ROOT/plm_centroids
ls -lh $ROOT/plm_centroids/

echo "=== Stage 2: augment training.tsv with PLM features ==="
date
python3 $BM/augment_training_with_plm.py \
    --in-dir $ROOT/m3_features \
    --out-dir $ROOT/m4_features \
    --plm-dir $ROOT/plm \
    --centroids-dir $ROOT/plm_centroids

echo "=== Stage 3: merge + split ==="
date
python3 $BM/merge_and_split.py \
    --in-dir $ROOT/m4_features \
    --out-train $ROOT/ml/train_m4.tsv \
    --out-valid $ROOT/ml/valid_m4.tsv \
    --valid-frac 0.2 --seed 42

echo "=== Stage 4: train LambdaMART ==="
date
python3 $BM/train_lambdamart.py \
    --train $ROOT/ml/train_m4.tsv \
    --valid $ROOT/ml/valid_m4.tsv \
    --out $ROOT/ml/model_m4.txt \
    --iters 500 --early-stop 30 \
    --feature-importance $ROOT/ml/importance_m4.json

echo "DONE"
date
