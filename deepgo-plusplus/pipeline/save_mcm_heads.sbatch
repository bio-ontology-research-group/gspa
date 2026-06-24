#!/bin/bash
#SBATCH --partition=freecycle-h100
#SBATCH --qos=freecycle
#SBATCH --requeue
#SBATCH --job-name=savemcm
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1 --cpus-per-task=16 --mem=128G --gres=gpu:1
set -eo pipefail
cd /mnt/data/u/hohndor/cafa6-plm && mkdir -p logs emb preds/hmcnn models
export PATH=$HOME/.pixi/bin:$PATH; export CONDA_OVERRIDE_CUDA=12.0
export TORCH_HOME=/mnt/data/u/hohndor/esm2; export HF_HOME=/mnt/data/u/hohndor/hf_cache
PIXI="pixi run --manifest-path $HOME/pixi-envs/gspa-glm/pixi.toml"

# nktest prostt5 embeddings (covers the clean-GT; small)
[ -f emb/test_prostt5_nk.npz ] || $PIXI python extract_embeddings.py --model prostt5 \
    --fasta seqs/nktest.fasta --out emb/test_prostt5_nk --batch-size 8

# save deployable MCM heads + nktest predictions for the re-freeze
$PIXI python train_head_hmcnn.py --loss mcm --mode test \
    --train-emb emb/train_prostt5.npz --test-emb emb/test_prostt5_nk.npz \
    --labels train_ipr.tsv --dag go-dag.tsv --emb-model prostt5 \
    --save-model models/head_prostt5_mcm.pt --out preds/hmcnn/head_prostt5_mcm_nk.tsv --epochs 30
echo "[$(date +%T)] prostt5 saved"

$PIXI python train_head_hmcnn.py --loss mcm --mode test \
    --train-emb emb/train_esm2_650m.npz --test-emb emb/test_esm2_650m.npz \
    --labels train_ipr.tsv --dag go-dag.tsv --emb-model esm2_650m \
    --save-model models/head_650m_mcm.pt --out preds/hmcnn/head_650m_mcm_nk.tsv --epochs 30
echo "[$(date +%T)] 650m saved"

# self-test the apply path (load .pt, score nktest, compare row count)
$PIXI python train_head_hmcnn.py --load-model models/head_650m_mcm.pt \
    --apply-emb emb/test_esm2_650m.npz --out preds/hmcnn/head_650m_mcm_apply.tsv
echo "[$(date +%T)] apply-path rows: train=$(wc -l < preds/hmcnn/head_650m_mcm_nk.tsv) apply=$(wc -l < preds/hmcnn/head_650m_mcm_apply.tsv)"
ls -l models/head_*_mcm.pt
echo "[$(date +%T)] SAVE-HEADS DONE"
