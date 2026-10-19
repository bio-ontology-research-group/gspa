#!/bin/bash
#SBATCH --job-name=p14-agg
#SBATCH --partition=batch
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
set -euo pipefail
OUT=${OUT:-/ibex/scratch/projects/c2014/hohndor/mags_repro_2026-10-19}
source $OUT/code/jobs/common.sh
module load python/3.11.0 >/dev/null 2>&1
LIST=${MAG_LIST:-$OUT/mag_list.tsv}
W=${WORK_DIR:-$WORK}
TAG=${TAG:-v2}
mkdir -p $RES
python3 $PY aggregate --mag-list $LIST --work $W --out $RES/pillar3_dark_matter_mags_${TAG}.tsv \
    --per-mag-out $RES/pillar3_dark_matter_mags_${TAG}_per_mag.tsv
# all fills, one file
first=1; for f in $W/*/*_fills.tsv; do if [[ $first == 1 ]]; then cat $f; first=0; else tail -n +2 $f; fi; done \
  | gzip > $RES/mag_gap_fills_${TAG}.tsv.gz
for d in $W/*/; do g=$(basename $d); if [[ -s $d/timing.tsv ]]; then sed "s/^/$g\t/" $d/timing.tsv; fi; done > $RES/timing_${TAG}.tsv
for d in $W/*/; do g=$(basename $d); if [[ -s $d/input_sha256.txt ]]; then sed "s/^/$g\t/" $d/input_sha256.txt; fi; done > $RES/input_sha256_per_mag_${TAG}.tsv
column -t -s$'\t' $RES/pillar3_dark_matter_mags_${TAG}.tsv | cut -c1-250
