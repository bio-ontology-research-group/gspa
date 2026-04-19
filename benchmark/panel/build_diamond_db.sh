#!/bin/bash
#SBATCH --job-name=panel-dmnd-db
#SBATCH --partition=debug
#SBATCH --output=/data/hohndor/gspa/proteomes/bench_gtdb30/dmnd-db-%j.out
#SBATCH --error=/data/hohndor/gspa/proteomes/bench_gtdb30/dmnd-db-%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8

set -eo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics

ROOT=/data/hohndor/gspa/proteomes/bench_gtdb30
BP=/data/hohndor/gspa/benchmark-py
cd $ROOT

echo "[1/3] Building panel exclusion list"
python3 panel/build_panel_excludes.py \
    --manifest panel_manifest.tsv \
    --sprot-fasta /data/hohndor/gspa/benchmark/uniprot_sprot.fasta \
    --out panel_excluded.txt
wc -l panel_excluded.txt

echo "[2/3] Filter Swiss-Prot minus panel"
python3 $BP/filter_fasta_by_exclude.py \
    --fasta /data/hohndor/gspa/benchmark/uniprot_sprot.fasta \
    --exclude panel_excluded.txt \
    --out reference_panel.fasta
echo "  kept: $(grep -c '^>' reference_panel.fasta)"

echo "[3/3] makedb"
diamond makedb --in reference_panel.fasta -d reference_panel --threads 8

echo "[verify] 0 panel-taxon accs should remain in DB headers"
grep "^>" reference_panel.fasta | awk '{sub(/^>[^|]*\|/, ""); sub(/\|.*/, ""); print}' | sort -u > /tmp/db_accs.txt
comm -12 <(sort panel_excluded.txt) /tmp/db_accs.txt | wc -l

echo DONE
date
