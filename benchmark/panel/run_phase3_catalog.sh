#!/usr/bin/env bash
#SBATCH --job-name=panel-p3
#SBATCH --partition=debug
#SBATCH --exclude=node003
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH -o /data/hohndor/gspa/proteomes/culture_panel/phase3/logs/p3-%j.out
#SBATCH -e /data/hohndor/gspa/proteomes/culture_panel/phase3/logs/p3-%j.err

# Phase 3 driver: rebuild the non-anchor cross-genome catalog over the
# expanded culture panel (97 HQ+rep genomes from culture_panel/phase2/).
#
# Steps
# -----
# 1. Concat per-genome proteomes into all_panel_proteins.faa (prefixed
#    with "<tag>:") and cluster with MMseqs2 easy-cluster at 50% id /
#    80% cov. Produce orthogroup_map_50.tsv in the format
#    <protein_id>\t<tag:rep>.
# 2. Flatten phase2/<tag>/integrated/ and phase2/<tag>/layout/ into
#    phase3/{integrated,layout}/ as symlinks, and build
#    phase3/panel_manifest.tsv (one tag per row).
# 3. Collect gap (rxn, ec) pairs across the panel to build
#    phase3/all_panel_gaps.jsonl — used as --restrict-to-gaps-file.
# 4. Invoke build_nonanchor_catalog.py over the panel.

set -eo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics
set -u

PANEL=/data/hohndor/gspa/proteomes/culture_panel
PH2=$PANEL/phase2
PH3=$PANEL/phase3
BM=/data/hohndor/gspa/benchmark-py
MANIFEST=$PANEL/phase2_manifest.tsv   # tag\tfasta_path

mkdir -p $PH3/{logs,ortho,integrated,layout,tmp}

TAGS=$(tail -n +2 $MANIFEST | cut -f1)
N=$(echo "$TAGS" | wc -l)
echo "panel tags: $N"
date

echo "=== step 1: concat proteomes + MMseqs2 cluster ==="
FAA=$PH3/ortho/all_panel_proteins.faa
if [[ ! -s $FAA ]]; then
    > $FAA
    for tag in $TAGS; do
        src=$PH2/$tag/prodigal/$tag.faa
        if [[ ! -s $src ]]; then
            echo "  [warn] missing $src — skipped"
            continue
        fi
        sed "s/^>/>${tag}:/" $src >> $FAA
    done
    echo "  total proteins: $(grep -c '^>' $FAA)"
fi

MAP=$PH3/ortho/orthogroup_map_50.tsv
if [[ ! -s $MAP ]]; then
    cd $PH3/ortho
    # Use node-local /tmp for MMseqs scratch — GlusterFS's sticky-T
    # split-brain entries trip up MMseqs' copy-after-merge step.
    MM_TMP=$(mktemp -d /tmp/mmseqs_panel.XXXXXX)
    mmseqs easy-cluster $FAA cluster50 $MM_TMP \
        --min-seq-id 0.5 -c 0.8 --cov-mode 0 \
        --threads 16 2>&1 | tail -20
    rm -rf $MM_TMP
    # cluster50_cluster.tsv: rep<TAB>member
    awk '{print $2 "\t" $1}' cluster50_cluster.tsv > $MAP
    echo "  orthogroup_map rows: $(wc -l < $MAP)"
    echo "  unique clusters:    $(cut -f2 $MAP | sort -u | wc -l)"
fi

echo
echo "=== step 2: flatten integrated/ + layout/ as symlinks ==="
for tag in $TAGS; do
    isrc=$PH2/$tag/integrated/${tag}_integrated.tsv
    idst=$PH3/integrated/${tag}_integrated.tsv
    lsrc=$PH2/$tag/layout/${tag}_layout.tsv
    ldst=$PH3/layout/${tag}_layout.tsv
    [[ -s $isrc && ! -e $idst ]] && ln -sf $isrc $idst
    [[ -s $lsrc && ! -e $ldst ]] && ln -sf $lsrc $ldst
done
echo "  integrated: $(ls $PH3/integrated/ | wc -l)"
echo "  layout:     $(ls $PH3/layout/ | wc -l)"

cat > $PH3/panel_manifest.tsv <<EOF
tag
EOF
for tag in $TAGS; do
    if [[ -s $PH3/integrated/${tag}_integrated.tsv && \
          -s $PH3/layout/${tag}_layout.tsv ]]; then
        echo "$tag"
    fi
done >> $PH3/panel_manifest.tsv
echo "  panel_manifest: $(wc -l < $PH3/panel_manifest.tsv)"

echo
echo "=== step 3: collect (rxn, ec) gaps across panel ==="
GAPS=$PH3/all_panel_gaps.jsonl
if [[ ! -s $GAPS ]]; then
    python3 <<PY
import csv
import json
from pathlib import Path

gaps = set()
with open("$MANIFEST") as f:
    f.readline()
    for line in f:
        tag = line.split("\t")[0].strip()
        tbl = Path(f"$PH2/{tag}/gapsmith/{tag}-all-Reactions.tbl")
        if not tbl.exists():
            continue
        with open(tbl) as g:
            rdr = csv.DictReader(g, delimiter="\t")
            for r in rdr:
                if r.get("status") not in ("bad_blast", "no_blast"):
                    continue
                ec = (r.get("ec") or "").strip()
                rxn = (r.get("rxn") or "").strip()
                if ec:
                    gaps.add((rxn, ec))

print(f"{len(gaps)} unique (rxn, ec) gaps across panel")
with open("$GAPS", "w") as f:
    for rxn, ec in sorted(gaps):
        f.write(json.dumps({"reaction_id": rxn, "ec_number": ec}) + "\n")
PY
fi
wc -l $GAPS

echo
echo "=== step 4: build non-anchor catalog ==="
OUT=$PANEL/nonanchor_catalog_panel.tsv
python3 $BM/build_nonanchor_catalog.py \
    --manifest $PH3/panel_manifest.tsv \
    --root $PH3 \
    --orthogroup-map $MAP \
    --reactions-tsv /data/hohndor/gspa/bin/gapsmith/data_merged/seed_reactions.tsv \
    --diffusion-tsv /data/hohndor/gspa/bin/gapsmith/data_merged/diffusion_mets.tsv \
    --ec-aliases-tsv /data/hohndor/gspa/bin/gapsmith/data_merged/seed_Enzyme_Class_Reactions_Aliases_unique.tsv \
    --ec2go /data/hohndor/gspa/reference/ec2go.txt \
    --tau 0.3 \
    --restrict-to-gaps-file $GAPS \
    --out $OUT 2>&1 | tail -20

wc -l $OUT
echo
echo "=== summary ==="
echo "panel size:       $N"
echo "orthogroup_map:   $(wc -l < $MAP)"
echo "panel_manifest:   $(wc -l < $PH3/panel_manifest.tsv)"
echo "catalog rows:     $(wc -l < $OUT)"
echo DONE; date
