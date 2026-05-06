#!/bin/bash
#SBATCH --job-name=rxn-catalog-int
#SBATCH --partition=debug
#SBATCH --exclude=node003
#SBATCH --output=/data/hohndor/gspa/proteomes/bench_gtdb30/catalog-int-%j.out
#SBATCH --error=/data/hohndor/gspa/proteomes/bench_gtdb30/catalog-int-%j.err
#SBATCH --time=06:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

set -eo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics

ROOT=/data/hohndor/gspa/proteomes/bench_gtdb30
cd $ROOT

# Catalog from *integrated posteriors* instead of raw DIAMOND claims.
# Held-out: mg1655. Reaction universe: restricted to EC-equivalents of mg1655 gaps.
python3 /data/hohndor/gspa/benchmark-py/build_catalog.py \
    --manifest panel_manifest.tsv \
    --root $ROOT \
    --orthogroup-map ortho/orthogroup_map_50.tsv \
    --reactions-tsv /data/hohndor/gspa/bin/gapsmith/data_merged/seed_reactions.tsv \
    --diffusion-tsv /data/hohndor/gspa/bin/gapsmith/data_merged/diffusion_mets.tsv \
    --ec-aliases-tsv /data/hohndor/gspa/bin/gapsmith/data_merged/seed_Enzyme_Class_Reactions_Aliases_unique.tsv \
    --ec2go /data/hohndor/gspa/reference/ec2go.txt \
    --exclude-tag mg1655 \
    --restrict-to-gaps-file /data/hohndor/gspa/proteomes/bench_ecoli/gaps/mg1655_gaps.jsonl \
    --tau 0.3 \
    --out catalog_panel_excl_mg1655_integrated.tsv

wc -l catalog_panel_excl_mg1655_integrated.tsv
echo DONE
date
