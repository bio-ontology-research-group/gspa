#!/bin/bash
#SBATCH --job-name=panel-cluster
#SBATCH --partition=debug
# node003 has an older CPU that SIGILLs on the metagenomics-env mmseqs
# binary. Exclude until a portable binary is vendored.
#SBATCH --exclude=node003
#SBATCH --output=/data/hohndor/gspa/proteomes/bench_gtdb30/cluster-%j.out
#SBATCH --error=/data/hohndor/gspa/proteomes/bench_gtdb30/cluster-%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8

set -eo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics

ROOT=/data/hohndor/gspa/proteomes/bench_gtdb30
cd $ROOT
mkdir -p ortho/work50 ortho/work90

# Combine all panel proteomes with tag-prefixed headers
echo "[1/3] combine proteomes"
ALL=ortho/all_panel_proteins.faa
if [[ ! -s $ALL ]]; then
    > $ALL
    for tag in $(tail -n +2 panel_manifest.tsv | cut -f1); do
        faa=proteomes/${tag}.faa
        [[ -s $faa ]] || continue
        awk -v t=$tag 'BEGIN{RS=">"; ORS=""} NR>1 {n=$1; sub(/^[^\n]*\n/,""); print ">"t":"n"\n"$0}' $faa >> $ALL
    done
fi
echo "  proteins: $(grep -c '^>' $ALL)"

# Stage to node-local to avoid GlusterFS tmp-copy errors
STAGE=/tmp/panel_cluster_$$
mkdir -p $STAGE
cp $ALL $STAGE/all.faa

echo "[2/3] MMseqs2 easy-cluster @ 50% id / 80% cov"
cd $STAGE
mmseqs easy-cluster all.faa cl50 tmp50 \
    --min-seq-id 0.5 -c 0.8 --cov-mode 0 \
    --threads 8 --remove-tmp-files 1

echo "[3/3] MMseqs2 easy-cluster @ 90% id / 80% cov"
mmseqs easy-cluster all.faa cl90 tmp90 \
    --min-seq-id 0.9 -c 0.8 --cov-mode 0 \
    --threads 8 --remove-tmp-files 1

cp cl50_cluster.tsv $ROOT/ortho/cluster_rep_member_50.tsv
cp cl90_cluster.tsv $ROOT/ortho/cluster_rep_member_90.tsv

# Build orthogroup_map.tsv (member <TAB> rep) used by downstream builds
awk -F'\t' '{sub(/^[^:]+:/, "", $2); print $2"\t"$1}' $ROOT/ortho/cluster_rep_member_50.tsv > $ROOT/ortho/orthogroup_map_50.tsv
awk -F'\t' '{sub(/^[^:]+:/, "", $2); print $2"\t"$1}' $ROOT/ortho/cluster_rep_member_90.tsv > $ROOT/ortho/orthogroup_map_90.tsv

cd $ROOT
rm -rf $STAGE
echo "[done] 50%: $(wc -l < ortho/orthogroup_map_50.tsv) rows, 90%: $(wc -l < ortho/orthogroup_map_90.tsv) rows"
date
