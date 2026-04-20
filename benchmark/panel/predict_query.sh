#!/usr/bin/env bash
#SBATCH --job-name=panel-p4
#SBATCH --partition=debug
#SBATCH --exclude=node003
#SBATCH --output=/data/hohndor/gspa/proteomes/culture_panel/phase4/logs/pq-%A_%a.out
#SBATCH --error=/data/hohndor/gspa/proteomes/culture_panel/phase4/logs/pq-%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=4

# Phase 4: per-query dark-matter prediction + cross-genome augmentation.
# Array task per genome:
#   1. predict_dark_matter.py over gapsmith gaps
#   2. augment with nonanchor_catalog (no self-ref filter — panel is
#      already ANI-95 dereplicated, so each genome contributes at most
#      1 to n_base counts; the ~1% optimism is negligible on 97 genomes)
#   3. Add reaction/EC names via augment_names.py

set -eo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics
set -u

MANIFEST=/data/hohndor/gspa/proteomes/culture_panel/phase2_manifest.tsv
PH2=/data/hohndor/gspa/proteomes/culture_panel/phase2
PH4=/data/hohndor/gspa/proteomes/culture_panel/phase4
BM=/data/hohndor/gspa/benchmark-py
REF=/data/hohndor/gspa/reference
GSD=/data/hohndor/gspa/bin/gapsmith/data_merged
CATALOG=/data/hohndor/gspa/proteomes/culture_panel/nonanchor_catalog_panel.tsv
ORTHO=/data/hohndor/gspa/proteomes/culture_panel/phase3/ortho/orthogroup_map_50.tsv

TAGS=($(tail -n +2 $MANIFEST | cut -f1))
tag=${TAGS[$SLURM_ARRAY_TASK_ID]}

WD=$PH4/$tag
mkdir -p $WD $PH4/logs

echo "=== [$tag] $(date) predict start ==="
cd $WD

# 1. predict dark matter
if [[ ! -s $WD/dark_matter.tsv ]]; then
    python3 $BM/cultures/predict_dark_matter.py \
        --tag $tag \
        --reactions-tbl $PH2/$tag/gapsmith/$tag-all-Reactions.tbl \
        --layout $PH2/$tag/layout/${tag}_layout.tsv \
        --integrated $PH2/$tag/integrated/${tag}_integrated.tsv \
        --reactions-tsv $GSD/seed_reactions.tsv \
        --diffusion-tsv $GSD/diffusion_mets.tsv \
        --ec-aliases-tsv $GSD/seed_Enzyme_Class_Reactions_Aliases_unique.tsv \
        --ec2go $REF/ec2go.txt \
        --top-k 5 --min-anchors 3 \
        --out $WD/dark_matter.tsv
fi
wc -l $WD/dark_matter.tsv

# 2. cross-genome augmentation.
# augment_cross_genome.py's mmseqs m8 input is an inter-genome search.
# Here the panel IS made of the same genomes we're augmenting, so the
# protein's orthogroup is a direct lookup in the ortho map (no mmseqs
# needed). We massage inputs to match the expected format by emitting
# a minimal m8 file where each predicted candidate maps to itself
# with the tag: prefix.
if [[ ! -s $WD/cand_self.m8 ]]; then
    awk -F"\t" -v tag=$tag "NR>1 {print tag\":\"\$5\"\t\"tag\":\"\$5\"\t100.0\t100\t100\t100\t1e-100\t300\"}" \
        $WD/dark_matter.tsv | sort -u > $WD/cand_self.m8
fi

if [[ ! -s $WD/dark_matter_augmented.tsv ]]; then
    python3 $BM/cultures/augment_cross_genome.py \
        --predictions $WD/dark_matter.tsv \
        --mmseqs-m8 $WD/cand_self.m8 \
        --panel-ortho-map $ORTHO \
        --catalog $CATALOG \
        --ec-aliases $GSD/seed_Enzyme_Class_Reactions_Aliases_unique.tsv \
        --out $WD/dark_matter_augmented.tsv
fi
wc -l $WD/dark_matter_augmented.tsv

# Naming is deferred to Phase 5 (single pass across all 97 outputs).

echo "=== [$tag] DONE $(date) ==="
