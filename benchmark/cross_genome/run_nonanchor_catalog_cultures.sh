#!/usr/bin/env bash
#SBATCH --job-name=nonanc-cult
#SBATCH --partition=debug
#SBATCH --exclude=node003
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 08:00:00
#SBATCH -o /data/hohndor/gspa/proteomes/bench_gtdb30/nonanc-cult-%j.out
#SBATCH -e /data/hohndor/gspa/proteomes/bench_gtdb30/nonanc-cult-%j.err

set -eo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics
set -u

ROOT=/data/hohndor/gspa/proteomes/bench_gtdb30
BM=/data/hohndor/gspa/benchmark-py
python3 $BM/build_nonanchor_catalog.py \
    --manifest $ROOT/panel_manifest.tsv \
    --root $ROOT \
    --orthogroup-map $ROOT/ortho/orthogroup_map_50.tsv \
    --reactions-tsv /data/hohndor/gspa/bin/gapsmith/data_merged/seed_reactions.tsv \
    --diffusion-tsv /data/hohndor/gspa/bin/gapsmith/data_merged/diffusion_mets.tsv \
    --ec-aliases-tsv /data/hohndor/gspa/bin/gapsmith/data_merged/seed_Enzyme_Class_Reactions_Aliases_unique.tsv \
    --ec2go /data/hohndor/gspa/reference/ec2go.txt \
    --tau 0.3 \
    --restrict-to-gaps-file /data/hohndor/gspa/proteomes/cultures/all_culture_gaps.jsonl \
    --out $ROOT/nonanchor_catalog_cultures.tsv
wc -l $ROOT/nonanchor_catalog_cultures.tsv
echo DONE; date
