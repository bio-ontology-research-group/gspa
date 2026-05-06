#!/usr/bin/env bash
#SBATCH --job-name=lro-stratify
#SBATCH --partition=debug
#SBATCH --exclude=node003
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH -o /data/hohndor/gspa/proteomes/bench_gtdb30/stratify-%j.out
#SBATCH -e /data/hohndor/gspa/proteomes/bench_gtdb30/stratify-%j.err

set -eo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics
set -u

ROOT=/data/hohndor/gspa/proteomes/bench_gtdb30
DIAMOND=/storage/miniforge3/envs/metagenomics/bin/diamond
WORK=$ROOT/homology_strat
mkdir -p $WORK

echo "=== Prepare panel (no mg1655) + mg1655 query FASTA ==="
date
awk 'BEGIN{RS=">"; ORS=""} NR>1 && $1 !~ /^mg1655:/ {print ">"$0}' \
    $ROOT/ortho/all_panel_proteins.faa > $WORK/panel_nomg.faa
grep -c "^>" $WORK/panel_nomg.faa

echo "=== Build DIAMOND DB ==="
date
$DIAMOND makedb --in $WORK/panel_nomg.faa -d $WORK/panel_nomg --threads 8

echo "=== Query mg1655 vs panel ==="
date
$DIAMOND blastp \
    -q $ROOT/proteomes/mg1655.faa \
    -d $WORK/panel_nomg \
    -o $WORK/mg1655_vs_panel.tsv \
    --outfmt 6 qseqid sseqid pident length evalue bitscore \
    --max-target-seqs 200 --evalue 1e-5 --threads 8 --quiet

wc -l $WORK/mg1655_vs_panel.tsv
echo DONE
date
