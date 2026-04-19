#!/usr/bin/env bash
#SBATCH --job-name=cult-ortho
#SBATCH --partition=debug
#SBATCH --exclude=node003
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 02:00:00
#SBATCH -o /data/hohndor/gspa/proteomes/cultures/logs/ortho-%j.out
#SBATCH -e /data/hohndor/gspa/proteomes/cultures/logs/ortho-%j.err

set -eo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics
set -u

ROOT=/data/hohndor/gspa/proteomes/cultures
PANEL_FAA=/data/hohndor/gspa/proteomes/bench_gtdb30/ortho/all_panel_proteins.faa
PANEL_MAP=/data/hohndor/gspa/proteomes/bench_gtdb30/ortho/orthogroup_map_50.tsv
WD=$ROOT/ortho
mkdir -p $WD/tmp

# Concatenate culture proteomes, prefixed by tag
cat > $WD/all_cultures.faa <<EOF
EOF
for tag in MR59-1 MR60-1 C-1.1 C-1.3; do
    sed "s/^>/>${tag}:/" $ROOT/$tag/prodigal/$tag.faa >> $WD/all_cultures.faa
done
echo "culture proteins: $(grep -c '^>' $WD/all_cultures.faa)"
echo "panel proteins:   $(grep -c '^>' $PANEL_FAA)"

# MMseqs2 easy-search: best hits from cultures vs panel
# 50% seq id / 80% query cov to align with panel cluster cutoffs
mmseqs easy-search \
    $WD/all_cultures.faa \
    $PANEL_FAA \
    $WD/cult_vs_panel.m8 \
    $WD/tmp \
    --min-seq-id 0.5 -c 0.8 --cov-mode 0 \
    --max-seqs 1 --threads 8 \
    --format-output "query,target,pident,qcov,tcov,evalue,bits" 2>&1 | tail -20

wc -l $WD/cult_vs_panel.m8
head -3 $WD/cult_vs_panel.m8
echo DONE
date
