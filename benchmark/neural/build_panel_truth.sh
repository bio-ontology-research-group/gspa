#!/usr/bin/env bash
#SBATCH --job-name=gspa-panel-truth
#SBATCH --partition=debug
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 02:00:00
#SBATCH -o /data/hohndor/gspa-neural/logs/panel-truth-%j.out
#SBATCH -e /data/hohndor/gspa-neural/logs/panel-truth-%j.err

# Build per-genome SwissProt GO + EC truth for the 21-genome panel.
# Artifacts under /data/hohndor/gspa-neural/panel/:
#   maps/<tag>_map.tsv          RefSeq_id <TAB> UniProt_acc
#   truth/<tag>_truth_exp.tsv   UniProt_acc <TAB> aspect <TAB> GO:NNNNNNN  (exp evidence)
#   truth/<tag>_truth_all.tsv   same, all evidence (incl IEA), exclude NOT
#   truth/<tag>_truth_refseq.tsv  RefSeq_id <TAB> aspect <TAB> GO:NNNNNNN  (joined)
#   truth/<tag>_ec_refseq.tsv     RefSeq_id <TAB> '' <TAB> EC:X.X.X.X       (from uniprot_sprot.dat)

set -euo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate gapfix

PANEL=/data/hohndor/gspa-neural/panel
REF=/data/hohndor/gspa/reference
BENCH=/data/hohndor/gspa/benchmark
SCRIPTS=/data/hohndor/gspa-neural/benchmark/panel

cd "$PANEL"
mkdir -p maps truth

# 1) RefSeq -> UniProt per genome
python "$SCRIPTS/build_refseq_uniprot_map.py" \
    --manifest panel_manifest.tsv \
    --genomes-dir fastas \
    --collab "$REF/gene_refseq_uniprotkb_collab.gz" \
    --out-dir maps

# 2) For each tag, extract its UniProt accessions (col 2 of the map)
ACC_ARGS=()
while IFS=$'\t' read tag rest; do
    [[ "$tag" == "tag" ]] && continue
    map="maps/${tag}_map.tsv"
    [[ -s $map ]] || { echo "skip $tag (no map)"; continue; }
    awk -F'\t' 'NR>0{print $2}' $map | sort -u > truth/${tag}_uniprot_accs.txt
    n=$(wc -l < truth/${tag}_uniprot_accs.txt)
    echo "  $tag: $n uniprot accs"
    ACC_ARGS+=(--accessions "${tag}:truth/${tag}_uniprot_accs.txt")
done < panel_manifest.tsv

# 3) Scan GOA once, emit all truth TSVs (UniProt-keyed, aspect as MF/BP/CC)
python "$SCRIPTS/extract_goa_dual.py" \
    --goa "$BENCH/goa_uniprot_all.gaf.gz" \
    "${ACC_ARGS[@]}" \
    --out-dir truth

# 4) Rekey truth by RefSeq (so it matches what predictors emit)
while IFS=$'\t' read tag rest; do
    [[ "$tag" == "tag" ]] && continue
    map="maps/${tag}_map.tsv"
    truth_exp="truth/${tag}_truth_exp.tsv"
    truth_all="truth/${tag}_truth_all.tsv"
    [[ -s $map && -s $truth_exp ]] || continue
    # join on UniProt acc; emit refseq-keyed
    for src in $truth_exp $truth_all; do
        base=$(basename $src .tsv)
        out="truth/${base}_refseq.tsv"
        python -c "
import sys
m = {}
for line in open('$map'):
    p = line.rstrip('\n').split('\t')
    if len(p) >= 2:
        m.setdefault(p[1], []).append(p[0])
out = set()
hdr = None
with open('$src') as fh:
    for line in fh:
        p = line.rstrip('\n').split('\t')
        if hdr is None:
            hdr = p; continue
        if len(p) < 3: continue
        acc, asp, go = p[0], p[1], p[2]
        for rs in m.get(acc, []):
            out.add((rs, asp, go))
with open('$out','w') as fh:
    fh.write('protein_id\taspect\tfunction_id\n')
    for r in sorted(out):
        fh.write('\t'.join(r)+'\n')
print('wrote', '$out', len(out), 'rows')
"
    done
done < panel_manifest.tsv

echo DONE; date
