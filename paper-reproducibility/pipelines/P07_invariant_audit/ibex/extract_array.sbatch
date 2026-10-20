#!/bin/bash
#SBATCH --job-name=p07_extract
#SBATCH --partition=batch
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --array=0-107%40
#SBATCH --output=/encrypted/e3008/gspa/meteor108_repro_2026-10-20/logs/extract_%A_%a.out
# P07: per-genome draft / METEOR / Grow-Support extraction with GSPA's own comparator code
# (scripts/meteor_threshold_repair.py, staged next to this script). Read-only on the METEOR tree.
#   W            work area (scripts/, data/, per_genome/, logs/)
#   METEOR_ROOT  archived third-party results (meteor_v8_evw_p2mu3_run/, meteor_v7_run/, paperA_2026/, psb_revision/)
#   VENV         python env with numpy/scipy/pandas (the P07 venv of 2026-10-19 has them)
set -o pipefail
W=${W:-/encrypted/e3008/gspa/meteor108_repro_2026-10-20}
METEOR_ROOT=${METEOR_ROOT:-/ibex/scratch/projects/c2014/kexin/funcarve}
VENV=${VENV:-/ibex/scratch/projects/c2014/hohndor/meteor108_repro_2026-10-19/venv}
module load python/3.11.0 >/dev/null 2>&1
source $VENV/bin/activate
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $W/scripts/panel108_gram.tsv)
GCA=$(echo "$LINE" | cut -f1); GRAM=$(echo "$LINE" | cut -f2)
echo "START $GCA $GRAM $(date)"
python3 -u $W/scripts/extract_meteor108_drafts.py --gca "$GCA" --gram "$GRAM" --outdir $W/per_genome \
  --data-dir $W/data --meteor-root $METEOR_ROOT
echo "END $GCA rc=$? $(date)"
