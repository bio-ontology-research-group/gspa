#!/bin/bash
#SBATCH --job-name=p13-eval
#SBATCH --partition=batch
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
set -euo pipefail
OUT=${OUT:-/ibex/scratch/projects/c2014/hohndor/bgc_repro_2026-10-19}
source $OUT/code/jobs/common.sh
module load python/3.11.0 >/dev/null 2>&1
LIST=${GENOME_LIST:-$OUT/genomes.tsv}
INTEG_DIR=${INTEG_DIR:-$INTEG}
TAG=${TAG:-v2}
python3 $BGC evaluate --genomes $LIST --work $WORK --integrated $INTEG_DIR --obo $REF/go.obo \
  --threshold 0.50 --n-boot ${NBOOT:-1000} \
  --out $RES/pillar2_3_antismash_bgc_specificity_${TAG}.tsv --per-genome $RES/bgc_per_genome_counts_${TAG}.tsv
# inventory: clusters / genomes / classes
cat $WORK/*/*_regions.tsv | awk -F'\t' '$1!="genome"' > $RES/bgc_regions_all_${TAG}.tsv
cat $WORK/*/*_timing.tsv | awk 'NR==1 || $1!="genome"' > $RES/bgc_timing_${TAG}.tsv
cat $RES/pillar2_3_antismash_bgc_specificity_${TAG}.tsv
