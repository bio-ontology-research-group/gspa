#!/bin/bash
#SBATCH --job-name=rxn-catalog
#SBATCH --partition=debug
#SBATCH --exclude=node003
#SBATCH --output=/data/hohndor/gspa/proteomes/bench_gtdb30/catalog-%j.out
#SBATCH --error=/data/hohndor/gspa/proteomes/bench_gtdb30/catalog-%j.err
#SBATCH --time=06:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

set -eo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics

ROOT=/data/hohndor/gspa/proteomes/bench_gtdb30
cd $ROOT

# Use claims as posteriors (faster than a full integrate pass per genome)
# Catalog holds out MG1655 so M2 evaluation on MG1655 is an out-of-panel test.
python3 /home/leechuck/Public/software/gspa/benchmark/cross_genome/build_catalog.py \
    --manifest panel_manifest.tsv \
    --root $ROOT \
    --orthogroup-map ortho/orthogroup_map_50.tsv \
    --reactions-tsv /data/hohndor/gspa/bin/gapsmith/data_merged/seed_reactions.tsv \
    --diffusion-tsv /data/hohndor/gspa/bin/gapsmith/data_merged/diffusion_mets.tsv \
    --ec-aliases-tsv /data/hohndor/gspa/bin/gapsmith/data_merged/seed_Enzyme_Class_Reactions_Aliases_unique.tsv \
    --ec2go /data/hohndor/gspa/reference/ec2go.txt \
    --exclude-tag mg1655 \
    --use-claims-as-posteriors \
    --restrict-to-gaps-file /data/hohndor/gspa/proteomes/bench_ecoli/gaps/mg1655_gaps.jsonl \
    --out catalog_panel_excl_mg1655.tsv

wc -l catalog_panel_excl_mg1655.tsv
echo DONE
date
