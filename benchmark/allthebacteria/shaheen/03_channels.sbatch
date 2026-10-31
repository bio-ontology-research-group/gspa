#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --exclude=nid00022
#SBATCH --time=24:00:00
#SBATCH --job-name=atbp-S3
#SBATCH --no-requeue
#SBATCH --signal=B:USR1@900
# S3 channels on representatives: a pull worker. S3_SLOTS chunk slots (default 24) of S3_CPUS CPUs (default 8)
# claim chunks (shards/*/chunks/cNNNN.faa) by an atomic mkdir of cNNNN.lock and run the EQ channel recipe on them:
#   DIAMOND 2.1.16 blastp vs experimental-GO Swiss-Prot (--evalue 1e-5 --max-target-seqs 5)
#   HMMER 3.3.2 hmmsearch Pfam-A --cut_ga
#   InterProScan (IPS_APPL, -goterms -iprlookup -dp) in the container (Java 21 of the package tree)
#   claims: bgc_benchmark.py build-claims + add_ec_claims.py (container python)
# A slot stops claiming once the job has run S3_STOP_AFTER seconds (default 64800 = 18 h) so that a running
# chunk can finish inside the 24 h limit; locks of jobs that are no longer in squeue are cleared at start, so a
# killed worker's chunk is picked up again by the next worker. Outputs next to the chunk: cNNNN.{diamond.tsv.gz,
# pfam.domtbl.gz,ips.tsv.gz,claims.jsonl.gz}, cNNNN.done; one unit receipt per chunk in receipts/S3/chunks/.
set -euo pipefail
source "${HERE_CODE:?}/common.sh"
make_local; install_finish S3; inode_guard
SLOTS=${S3_SLOTS:-24}; CPUS=${S3_CPUS:-8}; STOP_AFTER=${S3_STOP_AFTER:-64800}
export CPUS STOP_AFTER
version diamond "$TOOLS/bin/diamond" version
version hmmsearch bash -c "'$TOOLS/bin/hmmsearch' -h | sed -n 2p"
CTR_PWD=$LOCAL version interproscan ctr "$IPS_HOME/interproscan.sh" --version
mkdir -p "$RECEIPTS/S3/chunks"
# clear stale locks (owner job not queued or running)
if running=$(squeue -h -u "$USER" -o %A 2>/dev/null); then
  for l in "$SHARDS_DIR"/*/chunks/c*.lock; do
    [[ -d $l ]] || continue
    c=${l%.lock}; [[ -e $c.done ]] && continue
    [[ -n $(find "$l" -maxdepth 0 -mmin +10) ]] || continue   # a worker that is just starting may not have written owner yet
    owner=$(cat "$l/owner" 2>/dev/null || echo none)
    grep -qx "${owner%%:*}" <<< "$running" || { echo "stale lock $l (owner $owner)"; rm -rf "$l"; }
  done
fi

# Every command is checked explicitly: run_chunk is called from an if, where set -e does not apply.
run_chunk() {   # run_chunk <path/cNNNN.faa>
  local faa=$1 c=${1%.faa} s id W n ncl
  s=$(basename "$(dirname "$(dirname "$faa")")"); id=${s}_$(basename "$c")
  W=$LOCAL/chunks/$id; mkdir -p "$W/ipstmp"
  local TIMES_ROOT=$W/times   # this chunk's processes only (unit receipt); merged into the job's times below
  n=$(grep -c '^>' "$faa")
  echo "[$id] $n representatives start $(date -Is)"
  local t0; t0=$(date +%s.%N)
  timed channels diamond "$id" "$TOOLS/bin/diamond" blastp --query "$faa" --db "$SPROT_DMND" --threads "$CPUS" \
    --outfmt 6 --evalue 1e-5 --max-target-seqs 5 --quiet --tmpdir "$W" --out "$W/diamond.tsv" || return 1
  timed channels hmmsearch "$id" "$TOOLS/bin/hmmsearch" --domtblout "$W/pfam.domtbl" --noali --cpu "$CPUS" --cut_ga \
    "$PFAM_HMM" "$faa" > /dev/null || return 1
  sed '/^>/!s/\*//g' "$faa" > "$W/in.faa" || return 1
  CTR_PWD=$W timed channels interproscan "$id" ctr "$IPS_HOME/interproscan.sh" -i "$W/in.faa" -f tsv -o "$W/ips.tsv" \
    -appl "$IPS_APPL" -goterms -iprlookup -dp -cpu "$CPUS" -T "$W/ipstmp" > "$W/ips.log" 2>&1 \
    || { echo "[$id] interproscan failed"; tail -5 "$W/ips.log"; return 1; }
  rm -rf "$W/ipstmp" "$W/in.faa"
  CTR_PWD=$W timed claims build_claims "$id" ctr "$CTR_PY" "$CODE/benchmark/unified/pillar2_fba_lro_price149/bgc_benchmark.py" \
    build-claims --obo "$GO_OBO" --diamond "$W/diamond.tsv" --sprot-go "$SPROT_GO" --domtbl "$W/pfam.domtbl" \
    --pfam2go "$PFAM2GO" --ips "$W/ips.tsv" --out "$W/claims.jsonl" > "$W/claims.log" 2>&1 || return 1
  CTR_PWD=$W timed claims ec_claims "$id" ctr "$CTR_PY" "$CODE/benchmark/emptyquarter_mags/add_ec_claims.py" \
    --claims "$W/claims.jsonl" --diamond "$W/diamond.tsv" --ec2go "$EC2GO" --sprot-ec "$SPROT_EC" >> "$W/claims.log" 2>&1 || return 1
  ncl=$(wc -l < "$W/claims.jsonl")
  for f in diamond.tsv pfam.domtbl ips.tsv claims.jsonl; do gzip -1 -c "$W/$f" > "$c.$f.gz.part" || return 1; done
  for f in diamond.tsv pfam.domtbl ips.tsv claims.jsonl; do mv "$c.$f.gz.part" "$c.$f.gz"; done
  receipt write --step S3 --unit "$id" --kind unit --out "$RECEIPTS/S3/chunks/$id.json" \
    --times-root "$W/times" --job-start "$t0" --nodes 0 --count "n_representatives=$n" --count "n_claims=$ncl" \
    --param "cpus=$CPUS" --param "shard=$s" --input-sha "$faa" > /dev/null || return 1
  printf '%s\t%s\t%s\t%s\t%s\n' "$id" "$n" "$ncl" "${SLURM_JOB_ID:-local}" "$(date -Is)" > "$c.done"
  mkdir -p "$LOCAL/times" && cp -r "$W/times/." "$LOCAL/times/"
  rm -rf "$W"
  echo "[$id] done: $ncl claims $(date -Is)"
}
slot() {   # slot <k>: claim and run chunks until none is left or the stop time is reached
  local k=$1 faa c ndone=0
  while :; do
    if [[ -e $LOCAL/stop ]] || (( $(date +%s) - ${JOB_START%.*} > STOP_AFTER )); then break; fi
    faa=""
    for f in "$SHARDS_DIR"/*/chunks/c*.faa; do
      c=${f%.faa}
      [[ -e $c.done ]] && continue
      if mkdir "$c.lock" 2>/dev/null; then echo "${SLURM_JOB_ID:-local}:$k" > "$c.lock/owner"; faa=$f; break; fi
    done
    [[ -n $faa ]] || break
    if run_chunk "$faa"; then ndone=$((ndone + 1)); else echo "slot $k: chunk $faa failed" >&2; echo "$faa" >> "$LOCAL/failed_chunks"; fi
    [[ -e ${faa%.faa}.done ]] || rmdir "${faa%.faa}.lock" 2>/dev/null || rm -rf "${faa%.faa}.lock"
  done
  echo "slot $k: $ndone chunks"
}
trap 'touch "$LOCAL/stop"' USR1
phase channels start
pids=()
for (( k = 0; k < SLOTS; k++ )); do slot "$k" & pids+=($!); done
for p in "${pids[@]}"; do wait "$p" || true; done
wait || true   # a USR1 (15 min before the limit) interrupts wait; collect every slot before the receipt
phase channels end
nd=0; nr=0
for d in "$SHARDS_DIR"/*/chunks/c*.done; do
  [[ -e $d ]] || continue
  read -r _ n _ job _ < "$d"; [[ "$job" == "${SLURM_JOB_ID:-local}" ]] && { nd=$((nd + 1)); nr=$((nr + n)); }
done
left=0; for f in "$SHARDS_DIR"/*/chunks/c*.faa; do [[ -e ${f%.faa}.done ]] || left=$((left + 1)); done
nfail=$(nlines "$LOCAL/failed_chunks")
RARGS+=(--count "n_chunks=$nd" --count "n_representatives=$nr" --count "n_chunks_failed=$nfail" --count "n_chunks_left=$left"
        --param "slots=$SLOTS" --param "cpus_per_slot=$CPUS" --param "ips_appl=$IPS_APPL")
echo "worker done: $nd chunks ($nr representatives) here, $left chunks not done overall, $nfail failed here"
