#!/bin/bash
#SBATCH --job-name=panel-gapsmith
#SBATCH --partition=debug
#SBATCH --output=/data/hohndor/gspa/proteomes/bench_gtdb30/gapsmith-%A_%a.out
#SBATCH --error=/data/hohndor/gspa/proteomes/bench_gtdb30/gapsmith-%A_%a.err
#SBATCH --time=08:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=8
#SBATCH --array=0-29

set -eo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
# metagenomics env has diamond; gapsmith binary runs on glibc 2.35 as shipped
conda activate metagenomics

ROOT=/data/hohndor/gspa/proteomes/bench_gtdb30
cd $ROOT
mkdir -p gapsmith

TAGS=($(tail -n +2 panel_manifest.tsv | cut -f1))
tag=${TAGS[$SLURM_ARRAY_TASK_ID]}
faa=$ROOT/proteomes/${tag}.faa
[[ -s $faa ]] || { echo "missing $faa; skip"; exit 0; }

outdir=gapsmith/${tag}
mkdir -p $outdir
[[ -s $outdir/${tag}-all-Reactions.tbl ]] && { echo "$tag already done"; exit 0; }

echo "=== $tag: gapsmith find -p all ==="
date
GAPSMITH=/data/hohndor/gspa/bin/gapsmith/gapsmith-linux
DATA=/data/hohndor/gspa/bin/gapsmith/data_merged

$GAPSMITH --data-dir $DATA find -p all -t Bacteria -o $outdir -u all -A diamond -K 8 $faa

ls -la $outdir
date
echo DONE_$tag
