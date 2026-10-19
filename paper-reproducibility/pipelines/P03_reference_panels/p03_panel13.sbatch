#!/bin/bash
#SBATCH -J p03_panel13
#SBATCH -p debug
#SBATCH -c 16
#SBATCH --mem=88G
#SBATCH -w node007
#SBATCH -o /home/hohndor/p03logs/p03_panel13_%j.out
#SBATCH -e /home/hohndor/p03logs/p03_panel13_%j.err
# node007 serves /data over NFS and sees it as /srv/data (no /data mount there)
# P03 13-genome panel: build -> integrate -> evaluate -> audit (t=0.5 and t=F-max threshold) -> perturbation survivors
set -euo pipefail
PY=/storage/miniforge3/envs/metagenomics/bin/python3
export P03_DATA=${P03_DATA:-/srv/data}
CODE=${P03_DATA}/hohndor/gspa/proteomes/bench10/repro_2026-10-19/code
STEPS=${STEPS:-"build integrate eval audit audit_fmaxt perturb"}
echo "host=$(hostname) start=$(date -Is) steps=${STEPS}"
for s in ${STEPS}; do
  t0=$(date +%s)
  case $s in
    build)     $PY $CODE/p03_build.py ;;
    integrate) $PY $CODE/p03_integrate.py --threads 12 --xmx 72g ;;
    audit)     $PY $CODE/p03_audit.py --xmx 64g --variant t05 ;;
    audit_fmaxt) $PY $CODE/p03_audit.py --xmx 64g --variant fmaxt ;;
    eval)      $PY $CODE/p03_eval13.py ;;
    perturb)   $PY $CODE/p03_perturb_eval.py ;;
  esac
  echo "STEP $s wall=$(( $(date +%s) - t0 ))s"
done
echo "end=$(date -Is)"
