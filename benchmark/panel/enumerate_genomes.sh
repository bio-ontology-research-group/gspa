#!/usr/bin/env bash
# Scoped enumeration of genome FASTAs across isolates/, enrichment/,
# rh_sequencing/.  v3: every loop appends to $OUT.

set -euo pipefail

ROOT=/data/emptyquarter/sequencing-results
OUT=/data/hohndor/gspa/proteomes/culture_panel/phase1/genome_inventory.tsv
mkdir -p "$(dirname "$OUT")"

{
    echo -e "source_dir\tsample_id\tfasta_path\tsize_bytes"

    # ========== isolates ==========
    # <site>/<SAMPLE>/assembly/<NN>_<SAMPLE>*.{fasta,fa}
    # (depth-3: site/sample/assembly/file, so path goes to depth 4)
    # skip assembly/mags/... by requiring the file parent to be assembly/
    for site in site59 site60; do
        find -L $ROOT/isolates/$site -mindepth 3 -maxdepth 3 \
            -type f \( -name '*.fasta' -o -name '*.fa' -o -name '*.fna' \) \
            -size +100k 2>/dev/null \
        | while IFS= read -r p; do
            sample=$(echo "$p" | awk -F/ '{print $(NF-2)}')
            size=$(stat -c%s "$p" 2>/dev/null || echo 0)
            echo -e "isolates\t${site}/${sample}\t${p}\t${size}"
        done
    done

    # ========== enrichment ==========
    # A: enrichment/<C-N>/mags/<BIN>.fa  (short-read bins, depth 3)
    find -L $ROOT/enrichment -mindepth 3 -maxdepth 3 -type f \
        -path '*/mags/*.fa' -size +100k 2>/dev/null \
    | while IFS= read -r p; do
        cname=$(echo "$p" | awk -F/ '{print $(NF-2)}')
        bin=$(basename "$p" .fa)
        size=$(stat -c%s "$p" 2>/dev/null || echo 0)
        echo -e "enrichment\t${cname}/${bin}\t${p}\t${size}"
    done

    # B: enrichment/<C-N>/<tech>/assembly/mags/<BIN>/<BIN>_assembly.fa
    #    long-read curated MAGs  (depth 6 from enrichment/)
    find -L $ROOT/enrichment -mindepth 6 -maxdepth 6 -type f \
        -path '*/assembly/mags/*' -name '*_assembly.fa' \
        -size +100k 2>/dev/null \
    | while IFS= read -r p; do
        cname=$(echo "$p" \
            | awk -F"/enrichment/" '{print $2}' \
            | cut -d/ -f1)
        bin=$(echo "$p" | awk -F/ '{print $(NF-1)}')
        size=$(stat -c%s "$p" 2>/dev/null || echo 0)
        echo -e "enrichment\t${cname}/${bin}\t${p}\t${size}"
    done

    # ========== rh_sequencing ==========
    # pick assembly files, skip mags/ duplicates
    find -L $ROOT/rh_sequencing -type f -name '*_assembly.fa' \
        -size +100k -not -path '*/assembly/mags/*' 2>/dev/null \
    | while IFS= read -r p; do
        rhdir=$(echo "$p" \
            | awk -F"/rh_sequencing/" '{print $2}' \
            | cut -d/ -f1)
        fname=$(basename "$p" _assembly.fa)
        tech=$(echo "$p" | awk -F"/rh_sequencing/$rhdir/" '{print $2}' \
            | awk -F/ '{print $1}')
        if [[ "$tech" == "assembly" ]]; then
            sample="${rhdir}/${fname}"
        else
            sample="${rhdir}/${tech}/${fname}"
        fi
        size=$(stat -c%s "$p" 2>/dev/null || echo 0)
        echo -e "rh_sequencing\t${sample}\t${p}\t${size}"
    done

    # Other rh_sequencing FASTAs (not _assembly.fa named, not fastq/qc)
    find -L $ROOT/rh_sequencing -type f \
        \( -name '*.fa' -o -name '*.fasta' -o -name '*.fna' \) \
        -not -name '*_assembly.fa' \
        -not -path '*/fastq/*' \
        -not -path '*/qc/*' \
        -not -path '*/assembly/mags/*' \
        -size +100k 2>/dev/null \
    | while IFS= read -r p; do
        rhdir=$(echo "$p" \
            | awk -F"/rh_sequencing/" '{print $2}' \
            | cut -d/ -f1)
        fname=$(basename "$p" | sed 's/\.[^.]*$//')
        size=$(stat -c%s "$p" 2>/dev/null || echo 0)
        echo -e "rh_sequencing\t${rhdir}/${fname}\t${p}\t${size}"
    done
} > "$OUT"

# dedup by fasta_path
awk -F'\t' 'NR==1 || !seen[$3]++' "$OUT" > "${OUT}.tmp"
mv "${OUT}.tmp" "$OUT"

echo "=== total rows (incl. header) ==="
wc -l "$OUT"
echo "=== per source ==="
awk -F'\t' 'NR>1 {print $1}' "$OUT" | sort | uniq -c | sort -rn
echo "=== unique samples per source ==="
awk -F'\t' 'NR>1 {print $1"\t"$2}' "$OUT" | sort -u | cut -f1 | sort | uniq -c
