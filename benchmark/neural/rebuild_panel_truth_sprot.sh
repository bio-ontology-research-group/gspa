#!/usr/bin/env bash
#SBATCH --job-name=gspa-retruth-sprot
#SBATCH --partition=debug
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH -o /data/hohndor/gspa-neural/logs/retruth-sprot-%j.out
#SBATCH -e /data/hohndor/gspa-neural/logs/retruth-sprot-%j.err

# Rebuild panel RefSeq-keyed truth using SwissProt-filtered RefSeq↔UniProt
# map, so that curated (exp + IEA-from-sprot) annotations land on our
# panel proteins instead of TrEMBL-only IEA. Emits _sprot variants
# alongside the originals so existing truth isn't clobbered.

set -euo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate gapfix
export PYTHONUNBUFFERED=1

PANEL=/data/hohndor/gspa-neural/panel
SP=/data/hohndor/gspa-neural/sprot
REF=/data/hohndor/gspa/reference

[[ -s $SP/sprot_refseq_map.tsv ]] || { echo "missing $SP/sprot_refseq_map.tsv"; exit 2; }
[[ -s $SP/sprot_go.tsv ]]         || { echo "missing $SP/sprot_go.tsv"; exit 2; }
[[ -s $SP/sprot_ec.tsv ]]         || { echo "missing $SP/sprot_ec.tsv"; exit 2; }

# Build per-genome SwissProt-only RefSeq→UniProt map
while IFS=$'\t' read tag rest; do
    [[ "$tag" == "tag" ]] && continue
    gff=$PANEL/fastas/${tag}_genomic.gff
    [[ -s $gff ]] || { echo "skip $tag (no gff)"; continue; }
    python - <<EOF
import re, gzip
gff = "$gff"
sprot_map = "$SP/sprot_refseq_map.tsv"
out_map = "$PANEL/maps/${tag}_sprot_map.tsv"

refseq_ids = set()
opener = gzip.open if gff.endswith(".gz") else open
with opener(gff, "rt") as fh:
    for line in fh:
        if line.startswith("#"): continue
        parts = line.split("\t")
        if len(parts) < 9 or parts[2] != "CDS": continue
        m = re.search(r"protein_id=([^;]+)", parts[8])
        if m: refseq_ids.add(m.group(1))

keep = {}
with open(sprot_map) as fh:
    next(fh, None)  # header (might be absent in filtered file but safe)
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) < 2: continue
        rs, acc = p[0], p[1]
        if rs in refseq_ids and rs not in keep:
            keep[rs] = acc

with open(out_map, "w") as fh:
    fh.write("refseq_id\tuniprot_acc\n")
    for rs, acc in sorted(keep.items()):
        fh.write(f"{rs}\t{acc}\n")
print(f"$tag: {len(keep)} sprot-mapped refseqs (of {len(refseq_ids)} cds proteins)")
EOF
done < $PANEL/panel_manifest.tsv

# Rekey GO + EC truth via SwissProt-only maps
while IFS=$'\t' read tag rest; do
    [[ "$tag" == "tag" ]] && continue
    map="$PANEL/maps/${tag}_sprot_map.tsv"
    [[ -s $map ]] || continue
    python - <<EOF
from collections import defaultdict
tag = "$tag"
PANEL = "$PANEL"
SP = "$SP"

# refseq -> [uniprot accs]
m = defaultdict(list)
with open(f"{PANEL}/maps/{tag}_sprot_map.tsv") as fh:
    next(fh, None)
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2:
            m[p[0]].append(p[1])
# inverse: uniprot -> [refseqs]
inv = defaultdict(list)
for rs, accs in m.items():
    for a in accs:
        inv[a].append(rs)

def rekey(src, out, drop_header=True):
    rows = set()
    with open(src) as fh:
        if drop_header:
            next(fh, None)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3: continue
            acc, asp, term = p[0], p[1], p[2]
            for rs in inv.get(acc, []):
                rows.add((rs, asp, term))
    with open(out, "w") as fh:
        fh.write("protein_id\taspect\tfunction_id\n")
        for r in sorted(rows):
            fh.write("\t".join(r) + "\n")
    return len(rows)

n_go = rekey(f"{SP}/sprot_go.tsv", f"{PANEL}/truth/{tag}_truth_sprot_refseq.tsv")
# EC has format "acc\t\tEC:X" already (aspect blank)
n_ec = rekey(f"{SP}/sprot_ec.tsv", f"{PANEL}/truth/{tag}_ec_sprot_refseq.tsv")
print(f"{tag}: sprot GO={n_go}, sprot EC={n_ec}")
EOF
done < $PANEL/panel_manifest.tsv

# Propagate GO truth (ancestor-closed)
for f in $PANEL/truth/*_truth_sprot_refseq.tsv; do
    prop=${f%.tsv}_prop.tsv
    [[ -s $prop ]] && continue
    python /data/hohndor/gspa-neural/benchmark/neural/propagate_truth.py \
        --go-obo $REF/go.obo --in $f --out $prop
done

echo DONE; date
wc -l $PANEL/truth/*_sprot_*.tsv | tail -50
