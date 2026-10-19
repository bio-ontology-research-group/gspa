#!/bin/bash
#SBATCH -J p03_panel21
#SBATCH -p debug
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -w node007
#SBATCH -o /home/hohndor/p03logs/p03_panel21_%j.out
#SBATCH -e /home/hohndor/p03logs/p03_panel21_%j.err
# P03 21-genome panel: build (channels + ESM-2 t30 kNN) -> integrate (GSPA D+P) -> eval
set -euo pipefail
PY=/storage/miniforge3/envs/metagenomics/bin/python3
export P03_DATA=${P03_DATA:-/srv/data}
CODE=${P03_DATA}/hohndor/gspa/proteomes/bench10/repro_2026-10-19/code
STEPS=${STEPS:-"build integrate eval"}
echo "host=$(hostname) start=$(date -Is) steps=${STEPS}"
for s in ${STEPS}; do
  t0=$(date +%s)
  case $s in
    build)     $PY $CODE/p03_panel21.py build ;;
    integrate) $PY $CODE/p03_panel21.py integrate --threads 8 --xmx 28g ;;
    eval)      $PY $CODE/p03_panel21.py eval ;;
  esac
  echo "STEP $s wall=$(( $(date +%s) - t0 ))s"
done
echo "end=$(date -Is)"
