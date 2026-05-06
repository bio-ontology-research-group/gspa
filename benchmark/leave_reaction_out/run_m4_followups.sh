#!/usr/bin/env bash
#SBATCH --job-name=m4-setup
#SBATCH --partition=debug
#SBATCH --exclude=node003
#SBATCH -c 4
#SBATCH --mem=20G
#SBATCH -t 01:00:00
#SBATCH -o /data/hohndor/gspa/proteomes/bench_gtdb30/m4-followup-%j.out
#SBATCH -e /data/hohndor/gspa/proteomes/bench_gtdb30/m4-followup-%j.err

set -eo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics
set -u

ROOT=/data/hohndor/gspa/proteomes/bench_gtdb30
LRO=/data/hohndor/gspa/proteomes/bench_ecoli/leave_reaction_out

echo "=== (a) build ec_catalysts panel for strict Track A ==="
date
python3 /data/hohndor/gspa/benchmark-py/build_ec_catalysts.py \
    --manifest $ROOT/panel_manifest.tsv \
    --root $ROOT \
    --plm-dir $ROOT/plm \
    --ec2go /data/hohndor/gspa/reference/ec2go.txt \
    --exclude-tag mg1655 \
    --tau 0.3 \
    --out-dir $ROOT/plm_catalysts
ls -lh $ROOT/plm_catalysts/

echo "=== (b) build bsubtilis-held-out centroids ==="
date
python3 /data/hohndor/gspa/benchmark-py/build_plm_centroids.py \
    --manifest $ROOT/panel_manifest.tsv \
    --root $ROOT \
    --plm-dir $ROOT/plm \
    --ec2go /data/hohndor/gspa/reference/ec2go.txt \
    --exclude-tag bsubtilis \
    --tau 0.3 \
    --out-dir $ROOT/plm_centroids_nobsubt

echo "=== (c) build bsubtilis LRO cases ==="
date
python3 $LRO/build_testset.py \
    --reactions-tbl $ROOT/gapsmith/bsubtilis/bsubtilis-all-Reactions.tbl \
    --gff $ROOT/genomes/bsubtilis_genomic.gff \
    --map-tsv $ROOT/maps/bsubtilis_map.tsv \
    --operons $ROOT/operons/bsubtilis_operons.tsv \
    --ortho-map $ROOT/ortho/orthogroup_map_50.tsv \
    --truth $ROOT/truth/bsubtilis_truth_all.tsv \
    --ec2go /data/hohndor/gspa/reference/ec2go.txt \
    --genome-tag bsubtilis \
    --out $LRO/cases_bsubtilis.tsv
wc -l $LRO/cases_bsubtilis.tsv

echo DONE; date
