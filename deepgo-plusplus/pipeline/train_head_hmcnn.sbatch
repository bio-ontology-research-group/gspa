#!/bin/bash
#SBATCH --partition=freecycle-h100
#SBATCH --qos=freecycle
#SBATCH --requeue
#SBATCH --job-name=hmcnn
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1 --cpus-per-task=16 --mem=128G --gres=gpu:1
set -eo pipefail
cd /mnt/data/u/hohndor/cafa6-plm
export PATH=$HOME/.pixi/bin:$PATH; export CONDA_OVERRIDE_CUDA=12.0
PIXI="pixi run --manifest-path $HOME/pixi-envs/gspa-glm/pixi.toml"
mkdir -p preds/hmcnn

for MODEL in prostt5 esm2_3b; do
  for LOSS in bce mcm softreg; do
    OUT=preds/hmcnn/head_${MODEL}_${LOSS}.tsv
    echo "[$(date +%T)] === $MODEL / $LOSS (test mode) ==="
    $PIXI python train_head_hmcnn.py --loss $LOSS --mode test \
        --train-emb emb/train_${MODEL}.npz --test-emb emb/test_${MODEL}.npz \
        --labels train_ipr.tsv --dag go-dag.tsv --out $OUT --epochs 30
    echo "[$(date +%T)] $MODEL/$LOSS rows=$(wc -l < $OUT)"
  done
done
echo "[$(date +%T)] ALL DONE"
ls -l preds/hmcnn/
