#!/bin/bash
#SBATCH --partition=freecycle-h100
#SBATCH --qos=freecycle
#SBATCH --requeue
#SBATCH --signal=B:USR1@90
#SBATCH --job-name=h650m
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=06:00:00
#SBATCH --nodes=1 --cpus-per-task=16 --mem=128G --gres=gpu:1
set -eo pipefail
cd /mnt/data/u/hohndor/cafa6-plm && mkdir -p logs emb preds/hmcnn
export PATH=$HOME/.pixi/bin:$PATH; export CONDA_OVERRIDE_CUDA=12.0
export TORCH_HOME=/mnt/data/u/hohndor/esm2; export HF_HOME=/mnt/data/u/hohndor/hf_cache
PIXI="pixi run --manifest-path $HOME/pixi-envs/gspa-glm/pixi.toml"

# 1) extract ESM2-650M embeddings: train pop + nktest (covers the clean-GT proteins)
[ -f emb/train_esm2_650m.npz ] || $PIXI python extract_embeddings.py --model esm2_650m \
    --fasta seqs/uniprot_sprot_exp.fa --out emb/train_esm2_650m --batch-size 16
[ -f emb/test_esm2_650m.npz ]  || $PIXI python extract_embeddings.py --model esm2_650m \
    --fasta seqs/nktest.fasta --out emb/test_esm2_650m --batch-size 16

# 2) train bce + mcm heads (test mode -> predict nktest)
for LOSS in bce mcm; do
  OUT=preds/hmcnn/head_650m_${LOSS}.tsv
  echo "[$(date +%T)] === 650m / $LOSS ==="
  $PIXI python train_head_hmcnn.py --loss $LOSS --mode test \
      --train-emb emb/train_esm2_650m.npz --test-emb emb/test_esm2_650m.npz \
      --labels train_ipr.tsv --dag go-dag.tsv --out $OUT --epochs 30
  echo "[$(date +%T)] 650m/$LOSS rows=$(wc -l < $OUT)"
done
echo "[$(date +%T)] 650M DONE"
