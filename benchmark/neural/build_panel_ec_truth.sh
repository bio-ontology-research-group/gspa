#!/usr/bin/env bash
#SBATCH --job-name=gspa-panel-ec-truth
#SBATCH --partition=debug
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 01:00:00
#SBATCH -o /data/hohndor/gspa-neural/logs/panel-ec-truth-%j.out
#SBATCH -e /data/hohndor/gspa-neural/logs/panel-ec-truth-%j.err

set -euo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate gapfix
export PYTHONUNBUFFERED=1

PANEL=/data/hohndor/gspa-neural/panel
NEURAL=/data/hohndor/gspa-neural/benchmark/neural
SPROT=/data/hohndor/gspa/reference/uniprot_sprot.dat.gz

# 1) extract all EC rows for the union of panel uniprot accessions
cat $PANEL/truth/*_uniprot_accs.txt | sort -u > $PANEL/truth/all_panel_uniprot_accs.txt

python $NEURAL/extract_ec_truth.py \
    --dat $SPROT \
    --accessions $PANEL/truth/all_panel_uniprot_accs.txt \
    --out $PANEL/truth/all_panel_ec.tsv

# 2) per-genome EC truth, keyed on UniProt (then on RefSeq via map)
while IFS=$'\t' read tag rest; do
    [[ "$tag" == "tag" ]] && continue
    accs_file="$PANEL/truth/${tag}_uniprot_accs.txt"
    map_file="$PANEL/maps/${tag}_map.tsv"
    [[ -s "$accs_file" && -s "$map_file" ]] || continue

    # filter EC rows to this genome's UniProt accs
    awk -F'\t' 'NR==FNR{a[$0]=1; next} FNR==1 || $1 in a' \
        $accs_file $PANEL/truth/all_panel_ec.tsv > $PANEL/truth/${tag}_ec_uniprot.tsv

    # rekey to RefSeq
    python -c "
import sys
m = {}
for line in open('$map_file'):
    p = line.rstrip('\n').split('\t')
    if len(p) >= 2:
        m.setdefault(p[1], []).append(p[0])
out = set()
hdr = None
with open('$PANEL/truth/${tag}_ec_uniprot.tsv') as fh:
    for line in fh:
        p = line.rstrip('\n').split('\t')
        if hdr is None:
            hdr = p; continue
        if len(p) < 3: continue
        acc, asp, ec = p[0], p[1], p[2]
        for rs in m.get(acc, []):
            out.add((rs, asp, ec))
with open('$PANEL/truth/${tag}_ec_refseq.tsv','w') as fh:
    fh.write('protein_id\taspect\tfunction_id\n')
    for r in sorted(out):
        fh.write('\t'.join(r)+'\n')
print('wrote ${tag}_ec_refseq.tsv:', len(out), 'rows')
"
done < $PANEL/panel_manifest.tsv

echo DONE; date
wc -l $PANEL/truth/*_ec_refseq.tsv | tail -25
