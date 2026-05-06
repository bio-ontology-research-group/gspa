#!/bin/bash
#SBATCH --job-name=panel-ops-truth
#SBATCH --partition=debug
#SBATCH --output=/data/hohndor/gspa/proteomes/bench_gtdb30/ot-%A_%a.out
#SBATCH --error=/data/hohndor/gspa/proteomes/bench_gtdb30/ot-%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --array=0-29

set -eo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics

ROOT=/data/hohndor/gspa/proteomes/bench_gtdb30
BP=/data/hohndor/gspa/benchmark-py
cd $ROOT
mkdir -p operons truth gaps

TAGS=($(tail -n +2 panel_manifest.tsv | cut -f1))
tag=${TAGS[$SLURM_ARRAY_TASK_ID]}
gff=$ROOT/genomes/${tag}_genomic.gff
map=$ROOT/maps/${tag}_map.tsv

[[ -s $gff && -s $map ]] || { echo "missing $gff or $map; skip"; exit 0; }

echo "=== $tag: operons (from GFF, then remap via $map) ==="
python3 $BP/make_operons.py $gff operons/${tag}_operons_refseq.tsv 100
# Remap refseq -> uniprot
python3 - <<PY
import csv, sys
m = {}
with open("${map}") as f:
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p) >= 2:
            m[p[0]] = p[1]
with open("operons/${tag}_operons_refseq.tsv") as fin, open("operons/${tag}_operons.tsv", 'w') as fout:
    for line in fin:
        members = line.rstrip('\n').split('\t')
        up = [m[x] for x in members if x in m]
        if len(up) >= 2:
            fout.write('\t'.join(up) + '\n')
PY
echo "  operons: $(wc -l < operons/${tag}_operons.tsv)"

echo "=== $tag: GOA truth ==="
# extract_goa_dual writes out-dir/{tag}_truth_all.tsv + out-dir/{tag}_truth_exp.tsv
cut -f2 $map > truth/${tag}_accs.txt
python3 $BP/extract_goa_dual.py \
    --goa /data/hohndor/gspa/benchmark/goa_uniprot_all.gaf.gz \
    --accessions ${tag}:truth/${tag}_accs.txt \
    --out-dir truth/ || true
echo "  truth: $(wc -l < truth/${tag}_truth_all.tsv 2>/dev/null) all, $(wc -l < truth/${tag}_truth_exp.tsv 2>/dev/null) exp"

echo "=== $tag: layout ==="
python3 $BP/make_layout.py --gff $gff --map $map --out layout/${tag}_layout.tsv || true

date
echo DONE_$tag
