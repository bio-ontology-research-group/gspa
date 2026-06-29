#!/bin/bash
#SBATCH --partition=freecycle-h100
#SBATCH --qos=freecycle
#SBATCH --requeue
#SBATCH --job-name=hmcnn_oof
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=06:00:00
#SBATCH --nodes=1 --cpus-per-task=16 --mem=128G --gres=gpu:1
set -eo pipefail
cd /mnt/data/u/hohndor/cafa6-plm
export PATH=$HOME/.pixi/bin:$PATH; export CONDA_OVERRIDE_CUDA=12.0
PIXI="pixi run --manifest-path $HOME/pixi-envs/gspa-glm/pixi.toml"
mkdir -p preds/hmcnn

# OOF mcm train predictions (integrator training features, leak-free k-fold)
for MODEL in prostt5 esm2_3b; do
  OUT=preds/hmcnn/head_${MODEL}_mcm_trainoof.tsv
  echo "[$(date +%T)] === $MODEL / mcm (oof mode) ==="
  $PIXI python train_head_hmcnn.py --loss mcm --mode oof \
      --train-emb emb/train_${MODEL}.npz \
      --labels train_ipr.tsv --dag go-dag.tsv --kfold 5 \
      --out $OUT --epochs 30
  echo "[$(date +%T)] $MODEL/mcm-oof rows=$(wc -l < $OUT)"
done
echo "[$(date +%T)] OOF DONE"
ls -l preds/hmcnn/head_*_mcm_trainoof.tsv
