#!/bin/bash
#SBATCH --job-name=p10retune
#SBATCH --partition=debug
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --output=/data/hohndor/gspa/proteomes/bench10/phase10_retune/slurm-%j.out

set -euo pipefail
echo "host=$(hostname) date=$(date) jobid=${SLURM_JOB_ID:-?}"
mkdir -p /data/hohndor/gspa/proteomes/bench10/phase10_retune
cd /data/hohndor/gspa/proteomes/bench10
bash /data/hohndor/gspa-v1.5.0/benchmark/phase10_retune.sh
