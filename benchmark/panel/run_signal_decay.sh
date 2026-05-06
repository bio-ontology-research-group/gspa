#!/usr/bin/env bash
#SBATCH --job-name=panel-p6b
#SBATCH --partition=debug
#SBATCH --exclude=node003
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH -o /data/hohndor/gspa/proteomes/culture_panel/phase6/logs/p6b-%A_%a.out
#SBATCH -e /data/hohndor/gspa/proteomes/culture_panel/phase6/logs/p6b-%A_%a.err

# Phase 6b: rebuild the non-anchor catalog at smaller panel sizes so
# we can check that log_lr stabilizes (rather than inflating) as the
# panel grows.
#
# Array index  0 → N=29 subpanel
#              1 → N=60 subpanel
#              2 → N=97 reference (copy of the existing full catalog —
#                    written as signal_decay_n97.tsv for symmetric
#                    analysis)
#
# For each subpanel we take the first N tags from the full
# panel_manifest (sorted for determinism), keeping the full 97-genome
# ortho map (orthogroup definitions are the same across sizes — only
# the per-genome presence counts change).

set -eo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics
set -u

SIZES=(29 60 97)
N=${SIZES[$SLURM_ARRAY_TASK_ID]}

PANEL=/data/hohndor/gspa/proteomes/culture_panel
PH3=$PANEL/phase3
PH6=$PANEL/phase6
BM=/data/hohndor/gspa/benchmark-py
GSD=/data/hohndor/gspa/bin/gapsmith/data_merged

mkdir -p $PH6/logs $PH6/decay

# Deterministic sample: sorted by tag, first N
SUBMAN=$PH6/decay/panel_manifest_n${N}.tsv
head -1 $PH3/panel_manifest.tsv > $SUBMAN
tail -n +2 $PH3/panel_manifest.tsv | sort | head -n $N >> $SUBMAN
echo "subpanel size: $(($(wc -l < $SUBMAN) - 1))"

OUT=$PH6/decay/nonanchor_catalog_n${N}.tsv
if [[ ! -s $OUT ]]; then
    python3 $BM/build_nonanchor_catalog.py \
        --manifest $SUBMAN \
        --root $PH3 \
        --orthogroup-map $PH3/ortho/orthogroup_map_50.tsv \
        --reactions-tsv $GSD/seed_reactions.tsv \
        --diffusion-tsv $GSD/diffusion_mets.tsv \
        --ec-aliases-tsv $GSD/seed_Enzyme_Class_Reactions_Aliases_unique.tsv \
        --ec2go /data/hohndor/gspa/reference/ec2go.txt \
        --tau 0.3 \
        --restrict-to-gaps-file $PH3/all_panel_gaps.jsonl \
        --out $OUT 2>&1 | tail -15
fi

wc -l $OUT
echo DONE; date
