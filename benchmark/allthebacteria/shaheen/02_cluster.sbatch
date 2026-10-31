#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --exclude=nid00022
#SBATCH --time=06:00:00
#SBATCH --job-name=atbp-S2
#SBATCH --no-requeue
# S2 cluster for the shards in SHARDS: mmseqs2 18.8cc5c easy-linclust of the shard's pooled proteome at 90 %
# identity and 90 % bidirectional coverage (--min-seq-id 0.9 -c 0.9 --cov-mode 0, as EQ), S2_PARALLEL shards at a
# time with CORES_PER_NODE/S2_PARALLEL threads each, mmseqs temp in node-local tmpfs. The pilot clusters each shard
# on its own (the production plan adds a cross-shard pass on shard representatives; not measured here).
# Per shard: clu_cluster.tsv.gz (<rep>\t<member>), chunks/cNNNN.faa (CHUNK_SIZE representatives each), S2.done.
set -euo pipefail
source "${HERE_CODE:?}/common.sh"
make_local; install_finish S2; inode_guard
version mmseqs "$TOOLS/bin/mmseqs" version
todo=()
for s in $(shard_list); do
  [[ -s $SHARDS_DIR/$s/S1.done ]] || { echo "$s: S1 not done" >&2; exit 2; }
  [[ -s $SHARDS_DIR/$s/S2.done ]] && echo "[skip] $s" || todo+=("$s")
done
(( ${#todo[@]} )) || { echo "nothing to do"; exit 0; }
PAR=${S2_PARALLEL:-${#todo[@]}}; (( PAR > ${#todo[@]} )) && PAR=${#todo[@]}
TH=$(( CORES_PER_NODE / PAR ))
MEM=${S2_SPLIT_MEMORY:-$(( 300 / PAR ))G}
cluster() {   # cluster <shard>
  local s=$1 W=$LOCAL/$1 D=$SHARDS_DIR/$1
  mkdir -p "$W/tmp" "$D/chunks"
  timed cluster linclust "$s" "$TOOLS/bin/mmseqs" easy-linclust "$D/pooled.faa.gz" "$W/clu" "$W/tmp" \
    --min-seq-id 0.9 -c 0.9 --cov-mode 0 --threads "$TH" --split-memory-limit "$MEM" > "$W/linclust.log" 2>&1
  rm -rf "$W/tmp" "$W/clu_all_seqs.fasta"
  gzip -1 -c "$W/clu_cluster.tsv" > "$D/clu_cluster.tsv.gz.part" && mv "$D/clu_cluster.tsv.gz.part" "$D/clu_cluster.tsv.gz"
  rm -f "$D"/chunks/c*.faa
  awk -v n="$CHUNK_SIZE" -v dir="$D/chunks" 'BEGIN{c=0;k=n} /^>/{if (k>=n) {if (f) close(f); f=sprintf("%s/c%04d.faa", dir, c); c++; k=0} k++} {print > f}' "$W/clu_rep_seq.fasta"
  local nrep nprot nch
  nrep=$(grep -c '^>' "$W/clu_rep_seq.fasta"); nprot=$(wc -l < "$W/clu_cluster.tsv"); nch=$(ls "$D"/chunks/c*.faa | wc -l)
  printf '%s\t%s\t%s\t%s\t%s\n' "$s" "$nprot" "$nrep" "$nch" "$(date -Is)" > "$D/S2.done"
  gzip -c "$W/linclust.log" > "$D/S2.linclust.log.gz"
  echo "$s: $nprot proteins -> $nrep representatives in $nch chunks"
}
phase cluster start
pids=(); n=0
for s in "${todo[@]}"; do
  cluster "$s" & pids+=($!); n=$((n + 1))
  if (( n % PAR == 0 )); then for p in "${pids[@]}"; do wait "$p"; done; pids=(); fi
done
for p in "${pids[@]}"; do wait "$p"; done
phase cluster end
ng=0; np=0; nr=0; nc=0
for s in "${todo[@]}"; do
  read -r _ p r c _ < "$SHARDS_DIR/$s/S2.done"; np=$((np + p)); nr=$((nr + r)); nc=$((nc + c))
  ng=$((ng + $(wc -l < "$SHARDS_DIR/$s/counts.tsv")))
done
RARGS+=(--count "n_genomes=$ng" --count "n_proteins=$np" --count "n_representatives=$nr" --count "n_chunks=$nc"
        --count "n_shards=${#todo[@]}" --param "threads_per_shard=$TH" --param "chunk_size=$CHUNK_SIZE")
