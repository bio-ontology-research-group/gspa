#!/usr/bin/env bash
# Stage candidate genome FASTAs as symlinks, validating each for:
#   - target exists + readable
#   - first byte is '>'
#   - size >= 100k (resolving symlinks)
# Writes:
#   $STAGE/<genome_id>.fna  (symlinks)
#   $MAP   tsv: genome_id, source, sample, original_path
#   $REJ   tsv: genome_id, original_path, reason

set -euo pipefail

INV=${1:-/data/hohndor/gspa/proteomes/culture_panel/phase1/genome_inventory.tsv}
STAGE=${2:-/data/hohndor/gspa/proteomes/culture_panel/phase1/staged}
MAP=${3:-/data/hohndor/gspa/proteomes/culture_panel/phase1/genome_list.tsv}
REJ=${4:-/data/hohndor/gspa/proteomes/culture_panel/phase1/rejected.tsv}

mkdir -p "$STAGE"
find "$STAGE" -maxdepth 1 -type l -delete 2>/dev/null || true

python3 - "$INV" "$STAGE" "$MAP" "$REJ" <<'PY'
import csv
import os
import re
import sys
from pathlib import Path

inv, stage, map_out, rej_out = sys.argv[1:]
seen = set()

with open(inv) as fin, \
     open(map_out, "w") as fmap, \
     open(rej_out, "w") as frej:
    fmap.write("genome_id\tsource\tsample\toriginal_path\n")
    frej.write("genome_id\toriginal_path\treason\n")
    rdr = csv.DictReader(fin, delimiter="\t")
    for r in rdr:
        src = r["source_dir"]
        sample = r["sample_id"].replace("/", "_")
        path = r["fasta_path"]
        base = re.sub(r"\.(fna|fa|fasta)$", "", os.path.basename(path))
        gid = f"{src}__{sample}__{base}"
        if gid in seen:
            continue
        seen.add(gid)

        # Resolve through symlinks for size
        try:
            size = os.stat(path).st_size
        except (OSError, FileNotFoundError):
            frej.write(f"{gid}\t{path}\ttarget_missing\n")
            continue
        if size < 100_000:
            frej.write(f"{gid}\t{path}\tsize_below_100k({size})\n")
            continue

        try:
            with open(path, "rb") as f:
                first = f.read(1)
        except (OSError, FileNotFoundError):
            frej.write(f"{gid}\t{path}\tunreadable\n")
            continue
        if first != b">":
            frej.write(f"{gid}\t{path}\tinvalid_header({first.decode('latin1')!r})\n")
            continue

        fmap.write(f"{gid}\t{src}\t{sample}\t{path}\n")
        dst = os.path.join(stage, f"{gid}.fna")
        if os.path.lexists(dst):
            os.unlink(dst)
        os.symlink(path, dst)
PY

NREJ=$(( $(wc -l < $REJ) - 1 ))
NOK=$(( $(wc -l < $MAP) - 1 ))
echo "staged: $NOK valid FASTAs → $STAGE"
echo "rejected: $NREJ → see $REJ"
if [[ $NREJ -gt 0 ]]; then
    echo "rejection reasons:"
    awk -F'\t' 'NR>1 {gsub(/\(.*/, "", $3); print $3}' $REJ \
        | sort | uniq -c | sort -rn
fi
