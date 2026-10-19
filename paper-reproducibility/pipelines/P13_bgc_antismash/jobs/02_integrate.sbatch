#!/bin/bash
#SBATCH --job-name=p13-integrate
#SBATCH --partition=batch
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
# GSPA integrate in --batch-manifest mode for one configuration (env CFG = stage0_noisy_or | operon_only | gspa_A9).
set -euo pipefail
OUT=${OUT:-/ibex/scratch/projects/c2014/hohndor/bgc_repro_2026-10-19}
source $OUT/code/jobs/common.sh
module load python/3.11.0 >/dev/null 2>&1
LIST=${GENOME_LIST:-$OUT/genomes.tsv}
MDIR=${MANIFEST_DIR:-$OUT/manifests}
INTEG_DIR=${INTEG_DIR:-$INTEG}
mkdir -p $MDIR $INTEG_DIR/$CFG
python3 $BGC manifest --genomes $LIST --work $WORK --integrated $INTEG_DIR --out-dir $MDIR --taxon bacteria
$JAVA -Xmx56g -jar $JAR integrate --batch-manifest $MDIR/manifest_${CFG}.tsv --batch-threads 12 \
  --batch-stats-out $MDIR/batch_stats_${CFG}.tsv \
  --go-owl $REF/go.owl --pathways $REF/kegg_pathways.tsv --modules $REF/kegg_modules.tsv --ec2go $REF/ec2go.txt \
  --essential-profile bacteria --taxon bacteria --lite
echo "DONE $CFG $(date -Is)"
