#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --exclude=nid00022
#SBATCH --time=03:00:00
#SBATCH --job-name=atbp-S1
#SBATCH --no-requeue
# S1 proteins for the shards in SHARDS (s000:s001:...): extract the shard's genomes from its ATB archive
# (md5 checked against the catalogue) into node-local tmpfs, call genes with pyrodigal 3.7.1 (-p single; ATB has
# no Bakta annotation on Shaheen, which is the EQ rule's pyrodigal branch), S1_PARALLEL genomes at a time
# (default 192, one per physical core). Per shard on /scratch (no per-genome files):
#   shards/<s>/annot.tar        <G>.faa, <G>.gff3, <G>.fna of every genome
#   shards/<s>/pooled.faa.gz    every protein as <G>|<protein_id> (input of S2)
#   shards/<s>/counts.tsv       genome_id, proteins;   shards/<s>/S1.failed.tsv genomes without proteins
#   shards/<s>/S1.done
set -euo pipefail
source "${HERE_CODE:?}/common.sh"
make_local; install_finish S1; inode_guard
P=${S1_PARALLEL:-$CORES_PER_NODE}
version pyrodigal "$TOOLS/bin/pyrodigal" --version
version xz xz --version
todo=()
for s in $(shard_list); do [[ -s $SHARDS_DIR/$s/S1.done ]] && echo "[skip] $s" || todo+=("$s"); done
(( ${#todo[@]} )) || { echo "nothing to do"; exit 0; }

phase extract start
extract() {   # extract <shard>
  local s=$1 W=$LOCAL/$1 tar_xz md5 want
  mkdir -p "$W/raw" "$W/annot"
  shard_cols "$s" genome_id member > "$W/genomes.tsv"
  cut -f2 "$W/genomes.tsv" > "$W/members.txt"
  read -r tar_xz want < <(awk -F'\t' -v s="$s" 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; next} $c["shard"]==s{print $c["tar_xz"], $c["tar_xz_md5"]}' "$LISTS/shards.tsv")
  timed extract md5 "$s" md5sum "$ATB_ARCHIVES/$tar_xz" > "$W/md5.txt" &
  timed extract tar "$s" tar -xJf "$ATB_ARCHIVES/$tar_xz" -C "$W/raw" -T "$W/members.txt"
  wait
  md5=$(cut -d' ' -f1 "$W/md5.txt")
  [[ -z "$want" || "$md5" == "$want" ]] || { echo "$s: md5 $md5 != catalogue $want for $tar_xz" >&2; return 1; }
}
pids=(); for s in "${todo[@]}"; do extract "$s" & pids+=($!); done
fail=0; for p in "${pids[@]}"; do wait "$p" || fail=1; done
(( fail == 0 )) || { echo "archive extraction failed"; exit 3; }
phase extract end

phase genecall start
genecall() {  # genecall <shard> <genome> <member>
  local s=$1 g=$2 W=$LOCAL/$1
  timed genecall pyrodigal "$s.$g" "$TOOLS/bin/pyrodigal" -i "$W/raw/$3" -a "$W/annot/$g.faa.part" -f gff \
    -o "$W/annot/$g.gff3.part" -p single -j 1 > /dev/null 2> "$W/annot/$g.log" \
    && mv "$W/annot/$g.faa.part" "$W/annot/$g.faa" && mv "$W/annot/$g.gff3.part" "$W/annot/$g.gff3" \
    && mv "$W/raw/$3" "$W/annot/$g.fna" && rm -f "$W/annot/$g.log"
  return 0
}
export -f genecall
for s in "${todo[@]}"; do awk -F'\t' -v s="$s" '{print s "\t" $1 "\t" $2}' "$LOCAL/$s/genomes.tsv"; done \
  | xargs -P "$P" -L 1 bash -c 'genecall "$@"' _ 
phase genecall end

phase pack start
ng=0; np=0; nf=0
for s in "${todo[@]}"; do
  W=$LOCAL/$s D=$SHARDS_DIR/$s; mkdir -p "$D"
  : > "$W/counts.tsv"; : > "$W/failed.tsv"
  while IFS=$'\t' read -r g _m; do
    if [[ -s $W/annot/$g.faa ]]; then
      awk -v g="$g" '/^>/{sub(/^>/,""); split($0,a,/[ \t]/); print ">" g "|" a[1]; n++; next} {print} END{print g "\t" n+0 > "/dev/stderr"}' \
        "$W/annot/$g.faa" 2>> "$W/counts.tsv"
    else
      printf '%s\tno_proteins\t%s\n' "$g" "$(tr '\n' ' ' < "$W/annot/$g.log" 2>/dev/null | cut -c1-200)" >> "$W/failed.tsv"
      rm -f "$W/annot/$g".*
    fi
  done < "$W/genomes.tsv" | gzip -1 > "$D/pooled.faa.gz.part"
  ( cd "$W/annot" && tar -cf "$D/annot.tar.part" . )
  mv "$D/pooled.faa.gz.part" "$D/pooled.faa.gz"; mv "$D/annot.tar.part" "$D/annot.tar"
  cp "$W/counts.tsv" "$D/counts.tsv"; cp "$W/failed.tsv" "$D/S1.failed.tsv"
  g1=$(wc -l < "$W/counts.tsv"); p1=$(awk '{s+=$2} END{print s+0}' "$W/counts.tsv"); f1=$(wc -l < "$W/failed.tsv")
  printf '%s\t%s\t%s\t%s\n' "$s" "$g1" "$p1" "$(date -Is)" > "$D/S1.done"
  echo "$s: $g1 genomes, $p1 proteins, $f1 without proteins"
  ng=$((ng + g1)); np=$((np + p1)); nf=$((nf + f1))
done
phase pack end
RARGS+=(--count "n_genomes=$ng" --count "n_proteins=$np" --count "n_failed=$nf" --count "n_shards=${#todo[@]}" --param "parallel=$P")
