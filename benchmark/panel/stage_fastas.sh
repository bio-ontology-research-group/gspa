#!/usr/bin/env bash
# Stage all candidate genome FASTAs into a single flat directory with
# stable genome_id filenames so CheckM2 / GTDBtk / skani can all share
# the same input set.
#
# Reads genome_inventory.tsv (source, sample, path, size, depth).
# Writes symlinks $STAGE/<genome_id>.fna.
#
# Also writes genome_list.tsv mapping genome_id → original path.

set -euo pipefail

INV=${1:-/data/hohndor/gspa/proteomes/culture_panel/phase1/genome_inventory.tsv}
STAGE=${2:-/data/hohndor/gspa/proteomes/culture_panel/phase1/staged}
MAP=${3:-/data/hohndor/gspa/proteomes/culture_panel/phase1/genome_list.tsv}

mkdir -p "$STAGE"
echo -e "genome_id\tsource\tsample\toriginal_path" > "$MAP"

# Deduplicate by identical paths just in case find hit the same file
# twice (symlinks).  Keep the first occurrence per genome_id so we
# don't clobber.

awk -F'\t' -v stage="$STAGE" -v map="$MAP" '
NR == 1 { next }
{
    src=$1; sample=$2; path=$3
    # basename without .fna/.fa/.fasta
    n=split(path, p, "/")
    base=p[n]
    sub(/\.(fna|fa|fasta)$/, "", base)
    gsub(/\//, "_", sample)
    gid = src "__" sample "__" base
    if (seen[gid]++) next
    print gid "\t" src "\t" sample "\t" path >> map
    # symlink so we don'\''t duplicate the (potentially large) data
    cmd = "ln -sf \"" path "\" \"" stage "/" gid ".fna\""
    system(cmd)
}
' "$INV"

echo "staged: $(ls $STAGE | wc -l) FASTAs → $STAGE"
echo "map:    $(wc -l < $MAP) lines → $MAP"
