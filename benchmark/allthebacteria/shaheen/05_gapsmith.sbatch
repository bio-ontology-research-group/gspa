#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --exclude=nid00022
#SBATCH --time=08:00:00
#SBATCH --job-name=atbp-S5
#SBATCH --no-requeue
# S5 metabolic model per genome for the shards in SHARDS (after S1): gapsmith doall (find, find-transport, draft,
# medium, fill; DIAMOND aligner; -t Bacteria, -t Archaea for archaeal genomes) on the genome's proteome, as EQ
# 05_gapsmith.sbatch, including its medium-CSV repair (gapsmith 0.1.0 writes compound names with commas unquoted)
# and the draft-only fallback when the fill has no solution. S5_PARALLEL genomes at a time (default 48) with
# S5_THREADS threads each (default 4). Native binaries (static gapsmith, conda DIAMOND); reference read directly
# from /project (read-only; gapsmith builds its DIAMOND databases in the run directory).
# Per shard: shards/<s>/gapsmith.tar.gz (<G>/ run directory; <G>/<G>.xml -> <G>-filled.xml), gapsmith.status.tsv
# (genome, source = computed_gapsmith | computed_gapsmith_refill | gapsmith_fill_failed | gapsmith_failed), S5.done.
set -euo pipefail
source "${HERE_CODE:?}/common.sh"
make_local; install_finish S5; inode_guard
P=${S5_PARALLEL:-48}; TH=${S5_THREADS:-4}
export TH
version gapsmith "$GAPSMITH" --version
version diamond "$TOOLS/bin/diamond" version
todo=()
for s in $(shard_list); do
  [[ -s $SHARDS_DIR/$s/S1.done ]] || { echo "$s: S1 not done" >&2; exit 2; }
  [[ -s $SHARDS_DIR/$s/S5.done ]] && echo "[skip] $s" || todo+=("$s")
done
(( ${#todo[@]} )) || { echo "nothing to do"; exit 0; }
phase unpack start
for s in "${todo[@]}"; do
  mkdir -p "$LOCAL/$s/annot" "$LOCAL/$s/gs" "$LOCAL/$s/status"
  tar -xf "$SHARDS_DIR/$s/annot.tar" -C "$LOCAL/$s/annot" --wildcards '*.faa'
done
phase unpack end
model() {   # model <shard> <genome> <domain>
  local s=$1 g=$2 W=$LOCAL/$1 tax=Bacteria src=computed_gapsmith
  [[ "$3" == *rchaea* ]] && tax=Archaea
  local O=$W/gs/$g; mkdir -p "$O"
  export PATH=$TOOLS/bin:$PATH
  if ! timed model gapsmith "$s.$g" "$GAPSMITH" doall "$W/annot/$g.faa" -A diamond -K "$TH" -t "$tax" \
         --data-dir "$GAPSMITH_DATA" -f "$O" > "$O/gapsmith.log" 2>&1; then
    if [[ -s $O/$g-draft.gmod.cbor && -s $O/$g-medium.csv ]] && grep -q "malformed media row" "$O/gapsmith.log"; then
      "$PY_HOST" - "$O/$g-medium.csv" "$O/$g-medium.fixed.csv" <<'PY'
import sys
src, dst = sys.argv[1:3]
with open(src) as fi, open(dst, "w") as fo:
    for i, line in enumerate(fi):
        f = line.rstrip("\n").split(",")
        if i == 0 or len(f) <= 3:
            fo.write(line)
        else:
            fo.write(f[0] + "," + "; ".join(x.strip() for x in f[1:-1]) + "," + f[-1] + "\n")
PY
      if timed model gapsmith_fill "$s.$g" "$GAPSMITH" fill --media "$O/$g-medium.fixed.csv" -r "$O/$g-all-Reactions.tbl" \
           -K "$TH" --data-dir "$GAPSMITH_DATA" -o "$O" "$O/$g-draft.gmod.cbor" >> "$O/gapsmith.log" 2>&1; then
        src=computed_gapsmith_refill
      else src=gapsmith_fill_failed; fi
    elif [[ -s $O/$g-draft.gmod.cbor ]]; then src=gapsmith_fill_failed
    else src=gapsmith_failed; fi
  fi
  [[ "$src" == computed_* && -s $O/$g-filled.xml ]] && ln -sfn "$g-filled.xml" "$O/$g.xml"
  [[ "$src" == computed_* && ! -s $O/$g-filled.xml ]] && src=gapsmith_outputs_incomplete
  printf '%s\t%s\n' "$g" "$src" > "$W/status/$g"
  return 0
}
export -f model
phase model start
for s in "${todo[@]}"; do shard_cols "$s" genome_id domain | while IFS=$'\t' read -r g dom; do
  if [[ -s $LOCAL/$s/annot/$g.faa ]]; then printf '%s\t%s\t%s\n' "$s" "$g" "$dom"; fi; done; done \
  | xargs -P "$P" -L 1 bash -c 'model "$@"' _
phase model end
phase pack start
ng=0; nok=0
for s in "${todo[@]}"; do
  W=$LOCAL/$s D=$SHARDS_DIR/$s
  cat "$W"/status/* | sort > "$D/gapsmith.status.tsv"
  ( cd "$W/gs" && tar -czf "$D/gapsmith.tar.gz.part" . ) && mv "$D/gapsmith.tar.gz.part" "$D/gapsmith.tar.gz"
  g1=$(wc -l < "$D/gapsmith.status.tsv"); k1=$(grep -c $'\tcomputed_' "$D/gapsmith.status.tsv" || true)
  printf '%s\t%s\t%s\t%s\n' "$s" "$g1" "$k1" "$(date -Is)" > "$D/S5.done"
  ng=$((ng + g1)); nok=$((nok + k1)); echo "$s: $k1 of $g1 genomes with a gap-filled model"
done
phase pack end
RARGS+=(--count "n_genomes=$ng" --count "n_models_filled=$nok" --count "n_shards=${#todo[@]}" --param "parallel=$P" --param "threads=$TH")
