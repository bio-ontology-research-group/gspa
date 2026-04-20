#!/usr/bin/env bash
#SBATCH --job-name=panel-p4b
#SBATCH --partition=debug
#SBATCH --exclude=node003
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 03:00:00
#SBATCH -o /data/hohndor/gspa/proteomes/culture_panel/phase4/logs/augall-%j.out
#SBATCH -e /data/hohndor/gspa/proteomes/culture_panel/phase4/logs/augall-%j.err

# Phase 4B: augment all 97 dark_matter.tsv files in one job. Loads the
# filtered catalog once (log_lr >= 0.3 only, ~127K entries) into RAM.

set -eo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics
set -u

BM=/data/hohndor/gspa/benchmark-py
PANEL=/data/hohndor/gspa/proteomes/culture_panel

python3 $BM/panel/augment_all.py \
    --manifest $PANEL/phase2_manifest.tsv \
    --phase4-root $PANEL/phase4 \
    --catalog $PANEL/nonanchor_catalog_panel.tsv \
    --ortho-map $PANEL/phase3/ortho/orthogroup_map_50.tsv \
    --ec-aliases /data/hohndor/gspa/bin/gapsmith/data_merged/seed_Enzyme_Class_Reactions_Aliases_unique.tsv \
    --min-log-lr 0.3

echo DONE; date
