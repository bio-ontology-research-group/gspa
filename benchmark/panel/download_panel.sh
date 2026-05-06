#!/bin/bash
# Download RefSeq FASTA + GFF + proteome for all panel genomes.
# Usage: bash download_panel.sh <manifest.tsv> <outdir>
set -eo pipefail
MANIFEST=${1:-panel_manifest.tsv}
OUTDIR=${2:-genomes}
mkdir -p $OUTDIR
cd $OUTDIR

tail -n +2 ../$MANIFEST | while IFS=$'\t' read tag species phylum taxid acc_asm; do
    [[ -z "$tag" ]] && continue
    acc_dir=$(echo $acc_asm | awk -F_ '{print $1"_"$2}')
    p1=${acc_dir:4:3}; p2=${acc_dir:7:3}; p3=${acc_dir:10:3}
    base=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/${p1}/${p2}/${p3}/${acc_asm}
    echo "=== $tag: $acc_asm ==="
    for ext in genomic.fna.gz genomic.gff.gz protein.faa.gz; do
        out=${tag}_${ext%.gz}
        if [[ -s $out ]]; then
            echo "  $out: exists"
            continue
        fi
        curl -sSL --retry 3 "${base}/${acc_asm}_${ext}" -o ${tag}_${ext}
        if [[ -s ${tag}_${ext} ]]; then
            gunzip -f ${tag}_${ext}
            echo "  $out: $(wc -c < $out) bytes"
        else
            echo "  $out: DOWNLOAD FAILED"
            rm -f ${tag}_${ext}
        fi
    done
done
