#!/bin/bash
# Driver of the AllTheBacteria pilot on Shaheen III (gate G0). Run on a Shaheen login node from the mirrored
# code tree ($CODE/benchmark/allthebacteria/shaheen). Nothing is submitted without --submit: by default every
# command prints the node plan and the exact sbatch lines (dependencies shown as ${J_<job>} placeholders).
#
#   bash driver.sh plan                 # node plan + sbatch lines of `all` (default)
#   bash driver.sh all [--submit]       # S0 -> S1 -> S2 -> S3 workers -> S4;  S1 -> S5;  S4 + S5 -> S6
#   bash driver.sh smoke [--submit]     # `all` on 1 shard x 24 genomes, 2,000-representative chunks, under $OUT-smoke
#   bash driver.sh s3 [--submit]        # top-up S3 workers for the chunks not done yet (after S2)
#   bash driver.sh s4|s5|s6 [--submit]  # resubmit one step for every shard (done shards are skipped by the job)
#   bash driver.sh status               # queue, done markers, chunks left, scratch file count (read-only)
#   bash driver.sh accounting           # sacct of every job in $OUT/jobs.tsv (stdout; for summarise --sacct)
#   bash driver.sh summarise [ARGS]     # summarise_receipts.py on $OUT/receipts (ARGS passed through)
#
# Shape (environment, defaults in common.sh): N_SHARDS=20 GENOMES_PER_SHARD=500 CHUNK_SIZE=50000; node grouping
# GROUP=4 shards per S1/S2/S5 node job, S6_GROUP=2 shards per S6 node job, S3_NODES=4 worker nodes, S4 one node.
# Paths: JAR, CODE, REF, IPS_HOME, OUT, PKG_ROOT ... (common.sh). Every job runs on one exclusive workq node
# (192 cores, account k10018, --exclude=nid00022); S0 runs on the shared partition.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
export HERE_CODE=${HERE_CODE:-$here}
cmd=${1:-plan}; shift || true
SUBMIT=0; PASS=()
for a in "$@"; do if [[ "$a" == --submit ]]; then SUBMIT=1; else PASS+=("$a"); fi; done
if [[ "$cmd" == smoke ]]; then
  export OUT=${SMOKE_OUT:-${OUT:-/scratch/project/k10018/gspa-v3/atb-pilot}-smoke} N_SHARDS=1 GENOMES_PER_SHARD=${SMOKE_GENOMES:-24}
  export CHUNK_SIZE=2000 S3_NODES=1 S3_SLOTS=4 GENOMES_PER_JVM=12 PARALLEL_JVM=2 GROUP=1 S6_GROUP=1
  cmd=all
fi
source "$here/common.sh"
GROUP=${GROUP:-4}; S6_GROUP=${S6_GROUP:-2}; S3_NODES=${S3_NODES:-4}
JOBS_TSV=$OUT/jobs.tsv

shards_of_group() {   # shards_of_group <first> <size> -> s000:s001:...
  local i out=()
  for (( i = $1; i < $1 + $2 && i < N_SHARDS; i++ )); do out+=("$(shard_id "$i")"); done
  (IFS=:; echo "${out[*]}")
}
# sb <name> <dependency or ""> <shards or ""> <script> [extra sbatch args]: prints, or submits with --submit;
# sets J=<job id> (a ${J_<name>} placeholder when printing)
sb() {
  local name=$1 dep=$2 shards=$3 script=$4; shift 4
  local c=(sbatch --parsable --export=ALL --job-name="atbp-$name" --output="$LOGS/%x-%j.log")
  [[ -n "$dep" ]] && c+=(--dependency="$dep")
  c+=("$@" "$HERE_CODE/$script")
  if (( SUBMIT )); then
    J=$(SHARDS=$shards "${c[@]}")
    printf '%s\t%s\t%s\t%s\n' "$J" "atbp-$name" "$shards" "$(date -Is)" >> "$JOBS_TSV"
    echo "submitted atbp-$name: $J (shards ${shards:--}${dep:+, $dep})"
  else
    printf 'SHARDS=%s %s\n' "${shards:-}" "${c[*]}"
    J="\${J_${name//-/_}}"
  fi
}
deps() { local IFS=:; echo "$*"; }

plan() {
  local ng=$(( (N_SHARDS + GROUP - 1) / GROUP )) n6=$(( (N_SHARDS + S6_GROUP - 1) / S6_GROUP ))
  cat <<EOF
AllTheBacteria pilot: $N_SHARDS shards x $GENOMES_PER_SHARD genomes = $(( N_SHARDS * GENOMES_PER_SHARD )) genomes (seed $SEED)
  out $OUT   jar $JAR   code $CODE   ips $IPS_HOME
  S0 select            1 shared job (2 CPUs)
  S1 proteins          $ng node jobs x $GROUP shards (pyrodigal, ${S1_PARALLEL:-192} genomes at a time)
  S2 cluster           $ng node jobs x $GROUP shards (linclust 90/90 per shard)       after its S1 job
  S3 channels          $S3_NODES worker nodes x ${S3_SLOTS:-24} slots x ${S3_CPUS:-8} CPUs, chunks of $CHUNK_SIZE   after every S2 job
  S4 scatter           1 node job, all shards                                       after every S3 worker
  S5 gapsmith          $ng node jobs x $GROUP shards (${S5_PARALLEL:-48} x ${S5_THREADS:-4} CPUs)             after its S1 job
  S6 integrate         $n6 node jobs x $S6_GROUP shards (${PARALLEL_JVM:-12} JVMs x ${JVM_THREADS:-16} threads, ${GENOMES_PER_JVM:-96} genomes per JVM)   after S4 and its S5 job
EOF
}

case "$cmd" in
  plan|all)
    plan
    [[ "$cmd" == plan ]] && echo && echo "sbatch lines of 'all' (nothing submitted):"
    [[ "$cmd" == all && $SUBMIT == 0 ]] && echo && echo "dry run (add --submit to submit):"
    (( SUBMIT )) && mkdir -p "$LOGS" "$LISTS" "$RECEIPTS" "$SHARDS_DIR"
    sb S0 "" "" 00_select.sbatch; J0=$J
    J1=(); J2=(); J5=()
    for (( i = 0, g = 0; i < N_SHARDS; i += GROUP, g++ )); do
      sh=$(shards_of_group "$i" "$GROUP")
      sb "S1-g$g" "afterok:$J0" "$sh" 01_proteins.sbatch; J1[g]=$J
      sb "S2-g$g" "afterok:${J1[g]}" "$sh" 02_cluster.sbatch; J2[g]=$J
      sb "S5-g$g" "afterok:${J1[g]}" "$sh" 05_gapsmith.sbatch; J5[g]=$J
    done
    J3=()
    for (( k = 0; k < S3_NODES; k++ )); do sb "S3-w$k" "afterok:$(deps "${J2[@]}")" "" 03_channels.sbatch; J3+=("$J"); done
    sb S4 "afterok:$(deps "${J3[@]}")" "$(shards_of_group 0 "$N_SHARDS")" 04_scatter.sbatch; J4=$J
    for (( i = 0, h = 0; i < N_SHARDS; i += S6_GROUP, h++ )); do
      sh=$(shards_of_group "$i" "$S6_GROUP")
      # the S5 job(s) holding these shards
      d5=(); for (( k = i; k < i + S6_GROUP && k < N_SHARDS; k++ )); do d5+=("${J5[k / GROUP]}"); done
      d5=($(printf '%s\n' "${d5[@]}" | sort -u))
      sb "S6-g$h" "afterok:$(deps "$J4" "${d5[@]}")" "$sh" 06_integrate.sbatch
    done ;;
  s3)
    total=0; left=0
    for f in "$SHARDS_DIR"/*/chunks/c*.faa; do [[ -e $f ]] || continue; total=$((total + 1)); [[ -e ${f%.faa}.done ]] || left=$((left + 1)); done
    n=${S3_NODES_TOPUP:-$(( (left + ${S3_SLOTS:-24} - 1) / ${S3_SLOTS:-24} ))}
    echo "chunks: $total, not done: $left -> $n worker nodes"
    for (( k = 0; k < n; k++ )); do sb "S3-t$k" "" "" 03_channels.sbatch; done ;;
  s4) sb S4 "" "$(shards_of_group 0 "$N_SHARDS")" 04_scatter.sbatch ;;
  s5) for (( i = 0, g = 0; i < N_SHARDS; i += GROUP, g++ )); do sb "S5-g$g" "" "$(shards_of_group "$i" "$GROUP")" 05_gapsmith.sbatch; done ;;
  s6) for (( i = 0, h = 0; i < N_SHARDS; i += S6_GROUP, h++ )); do sb "S6-g$h" "" "$(shards_of_group "$i" "$S6_GROUP")" 06_integrate.sbatch; done ;;
  status)
    squeue -u "$USER" -o '%.12i %.16j %.9T %.10M %.6D %R' | awk 'NR==1 || /atbp-/'
    for m in S1 S2 S4 S5 S6; do printf '%s.done %s of %s shards\n' "$m" "$(ls "$SHARDS_DIR"/*/"$m".done 2>/dev/null | wc -l)" "$N_SHARDS"; done
    total=0; left=0
    for f in "$SHARDS_DIR"/*/chunks/c*.faa; do [[ -e $f ]] || continue; total=$((total + 1)); [[ -e ${f%.faa}.done ]] || left=$((left + 1)); done
    echo "chunks: $total, not done: $left"
    echo "files under $OUT: $(find "$OUT" -xdev | wc -l)"
    (command -v kuq >/dev/null && kuq) || lfs quota -u "$USER" /scratch || true ;;
  accounting)
    ids=$(cut -f1 "$JOBS_TSV" | grep -E '^[0-9]+$' | paste -sd,)
    sacct -X -P -n -j "$ids" -o JobID,JobName,NNodes,ElapsedRaw,State ;;
  summarise)
    "$PY_HOST" "$HERE_CODE/summarise_receipts.py" "$RECEIPTS" "${PASS[@]}" ;;
  *) sed -n '2,15p' "$0"; exit 1 ;;
esac
