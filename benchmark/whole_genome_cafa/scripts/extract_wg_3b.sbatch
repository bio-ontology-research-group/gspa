#!/bin/bash
#SBATCH --partition=freecycle-h200,freecycle-h100,batch-h100
#SBATCH --qos=freecycle
#SBATCH --requeue
#SBATCH --job-name=wg_esm2_3b
#SBATCH --output=/mnt/data/u/hohndor/gspa-wg-cafa/logs/%x-%A_%a.out
#SBATCH --error=/mnt/data/u/hohndor/gspa-wg-cafa/logs/%x-%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=8 --mem=64G
#SBATCH --gres=gpu:1
#SBATCH --array=0-3
set -eo pipefail
export PATH=$HOME/.pixi/bin:$PATH
export CONDA_OVERRIDE_CUDA=12.6
export TORCH_HOME=/mnt/data/u/hohndor/esm2
export HF_HOME=/mnt/data/u/hohndor/hf_cache
T=${SLURM_ARRAY_TASK_ID}
cd /home/hohndor/pixi-envs/gspa-glm
/home/hohndor/.pixi/bin/pixi run --frozen python /mnt/data/u/hohndor/cafa6-plm/extract_embeddings.py     --model esm2_3b     --fasta /mnt/data/u/hohndor/gspa-wg-cafa/benchmark/whole_genome_cafa/data/wg_missing_3b.fasta     --out /mnt/data/u/hohndor/gspa-wg-cafa/benchmark/whole_genome_cafa/data/wg_missing_3b_${T}     --shard ${T}/4     --batch-size 12
