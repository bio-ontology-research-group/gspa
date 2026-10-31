#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --exclude=nid00022
#SBATCH --time=04:00:00
#SBATCH --job-name=atbp-S4
#SBATCH --no-requeue
# S4 scatter for the shards in SHARDS (after every S3 chunk of those shards is done): claims of each representative
# copied to its members (atb_scatter.py, one process per shard, S4_PARALLEL at a time), one JSONL per genome in
# node-local tmpfs, packed as shards/<s>/claims.tar (<G>.claims.jsonl) with shards/<s>/coverage.tsv and S4.done.
set -euo pipefail
source "${HERE_CODE:?}/common.sh"
make_local; install_finish S4; inode_guard
version python "$PY_HOST" --version
todo=()
for s in $(shard_list); do
  [[ -s $SHARDS_DIR/$s/S2.done ]] || { echo "$s: S2 not done" >&2; exit 2; }
  [[ -s $SHARDS_DIR/$s/S4.done ]] && echo "[skip] $s" || todo+=("$s")
done
(( ${#todo[@]} )) || { echo "nothing to do"; exit 0; }
PAR=${S4_PARALLEL:-${#todo[@]}}
scatter_shard() {
  local s=$1 W=$LOCAL/$1 D=$SHARDS_DIR/$1
  mkdir -p "$W/claims"
  timed scatter atb_scatter "$s" "$PY_HOST" "$HERE_CODE/atb_scatter.py" --shard-dir "$D" --out-dir "$W/claims" --coverage "$D/coverage.tsv"
  ( cd "$W/claims" && tar -cf "$D/claims.tar.part" . ) && mv "$D/claims.tar.part" "$D/claims.tar"
  printf '%s\t%s\t%s\n' "$s" "$(awk 'NR>1{s+=$4} END{print s+0}' "$D/coverage.tsv")" "$(date -Is)" > "$D/S4.done"
  rm -rf "$W"
}
phase scatter start
pids=(); n=0
for s in "${todo[@]}"; do
  scatter_shard "$s" & pids+=($!); n=$((n + 1))
  if (( n % PAR == 0 )); then for p in "${pids[@]}"; do wait "$p"; done; pids=(); fi
done
for p in "${pids[@]}"; do wait "$p"; done
phase scatter end
ng=0; nc=0
for s in "${todo[@]}"; do ng=$((ng + $(tail -n +2 "$SHARDS_DIR/$s/coverage.tsv" | wc -l))); nc=$((nc + $(cut -f2 "$SHARDS_DIR/$s/S4.done"))); done
RARGS+=(--count "n_genomes=$ng" --count "n_claims=$nc" --count "n_shards=${#todo[@]}")
