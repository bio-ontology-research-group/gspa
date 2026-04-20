#!/usr/bin/env bash
# Phase 1 orchestrator: given genome_inventory.tsv from
# enumerate_genomes.sh, stage the FASTAs and submit CheckM2, GTDBtk,
# and skani as parallel SLURM jobs.  Each writes its own output into
# $ROOT/{checkm2,gtdbtk,skani}.  Once all three complete, run
# build_genome_manifest.py to emit genome_manifest.tsv.

set -euo pipefail

ROOT=/data/hohndor/gspa/proteomes/culture_panel/phase1
BM=/data/hohndor/gspa/benchmark-py/panel

mkdir -p $ROOT/logs
cd $ROOT

# 1. Stage FASTAs as symlinks under $ROOT/staged/<genome_id>.fna
echo "=== Stage ==="
bash $BM/stage_fastas.sh \
    $ROOT/genome_inventory.tsv \
    $ROOT/staged \
    $ROOT/genome_list.tsv

N=$(ls $ROOT/staged | wc -l)
echo "staged $N genomes."

if [[ $N -eq 0 ]]; then
    echo "ERROR: no genomes staged. Aborting."
    exit 1
fi

# 2. Submit parallel QC jobs
echo "=== Submit ==="
JCK=$(sbatch --parsable $BM/run_checkm2.sh)
JGT=$(sbatch --parsable $BM/run_gtdbtk.sh)
JSK=$(sbatch --parsable $BM/run_skani_derep.sh)
echo "checkm2:  $JCK"
echo "gtdbtk:   $JGT"
echo "skani:    $JSK"

# 3. Build manifest after all three finish
JMAN=$(sbatch --parsable \
    --dependency=afterok:$JCK:$JGT:$JSK \
    --job-name=panel-manifest \
    --partition=debug \
    --exclude=node003 \
    -c 2 --mem=8G -t 01:00:00 \
    -o $ROOT/logs/manifest-%j.out \
    -e $ROOT/logs/manifest-%j.err \
    --wrap="source /storage/miniforge3/etc/profile.d/conda.sh && \
conda activate metagenomics && \
python3 $BM/build_genome_manifest.py \
    --inventory $ROOT/genome_inventory.tsv \
    --checkm2-out $ROOT/checkm2/quality_report.tsv \
    --gtdbtk-bac120 $ROOT/gtdbtk/gtdbtk.bac120.summary.tsv \
    --gtdbtk-ar53 $ROOT/gtdbtk/gtdbtk.ar53.summary.tsv \
    --skani-clusters $ROOT/skani/clusters.tsv \
    --count-contigs \
    --out $ROOT/genome_manifest.tsv && \
echo 'manifest quality tiers:' && \
awk -F'\\t' 'NR>1 {t[\$NF]++} END {for (k in t) print k, t[k]}' \
    $ROOT/genome_manifest.tsv")

echo "manifest job: $JMAN  (runs after $JCK $JGT $JSK)"
echo
echo "watch with:"
echo "  squeue -j $JCK,$JGT,$JSK,$JMAN -u hohndor"
