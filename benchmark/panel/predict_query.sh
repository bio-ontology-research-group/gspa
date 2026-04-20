#!/usr/bin/env bash
#SBATCH --job-name=panel-p4a
#SBATCH --partition=debug
#SBATCH --exclude=node003
#SBATCH --output=/data/hohndor/gspa/proteomes/culture_panel/phase4/logs/pq-%A_%a.out
#SBATCH --error=/data/hohndor/gspa/proteomes/culture_panel/phase4/logs/pq-%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2

# Phase 4A: run predict_dark_matter.py per genome. Fast & small.
# Augmentation is deferred to Phase 4B (single job, single catalog
# load), since the 128M-row catalog is 15GB in RAM and can't be loaded
# by every array task.

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

TAGS=($(tail -n +2 $MANIFEST | cut -f1))
tag=${TAGS[$SLURM_ARRAY_TASK_ID]}

WD=$PH4/$tag
mkdir -p $WD

echo "=== [$tag] $(date) predict start ==="
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
echo "=== [$tag] DONE $(date) ==="
