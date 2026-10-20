#!/bin/bash
#SBATCH --job-name=eq-chunk
#SBATCH --partition=batch
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
# One chunk of cluster representatives per array task (line SLURM_ARRAY_TASK_ID+1 of lists/chunks.txt):
# DIAMOND vs experimental-GO Swiss-Prot, HMMER Pfam-A, InterProScan -> claims JSONL keyed by representative id.
# Idempotent: finished channels are skipped.
set -euo pipefail
# start from a clean module environment: the submitting shell may export LD_LIBRARY_PATH/PYTHONPATH of
# loaded modules (python/gcc), which breaks the conda tools and the static-linked gapsmith binary
module purge >/dev/null 2>&1 || true; unset LD_LIBRARY_PATH PYTHONPATH PYTHONHOME
source ${OUT:-/encrypted/e3008/gspa/eq_mags}/code/common.sh
module load python/3.11.0 diamond/2.1.16 hmmer/3.3.2 >/dev/null 2>&1
FAA=$(sed -n "$(( SLURM_ARRAY_TASK_ID + 1 ))p" $LISTS/chunks.txt); [[ -n "$FAA" ]] || { echo "no chunk"; exit 0; }
C=$(basename $FAA .faa); D=$CHUNKS/$C; mkdir -p $D/times; cd $D
T=${SLURM_CPUS_PER_TASK:-8}
[[ -s DONE ]] && { echo "[skip] $C done"; exit 0; }
N=$(grep -c '^>' $FAA)
echo "== $C: $N representatives $(date -Is) host=$(hostname) cpus=$T"
TV="/usr/bin/time -v"; run() { local n=$1; shift; $TV -o times/$n.time "$@"; }
[[ -s diamond.tsv || -s diamond.tsv.gz ]] || { run diamond diamond blastp --query $FAA --db $REF/sprot_exp.dmnd --threads $T \
    --outfmt 6 --evalue 1e-5 --max-target-seqs 5 --quiet --out diamond.tsv.part && mv diamond.tsv.part diamond.tsv; }
[[ -s pfam.domtbl || -s pfam.domtbl.gz ]] || { run hmmsearch hmmsearch --domtblout pfam.domtbl.part --noali --cpu $T --cut_ga \
    $PFAM_HMM $FAA > /dev/null && mv pfam.domtbl.part pfam.domtbl; }
if [[ ! -s ips.tsv && ! -s ips.tsv.gz ]]; then
  module load $IPS_MODULE >/dev/null 2>&1
  command -v interproscan.sh >/dev/null || { echo "interproscan.sh not on PATH after module load $IPS_MODULE"; exit 1; }
  TMPD=${TMPDIR:-/tmp}/ips_${C}_$$; mkdir -p $TMPD
  sed '/^>/!s/\*//g' $FAA > $TMPD/in.faa
  run interproscan interproscan.sh -i $TMPD/in.faa -f tsv -o ips.tsv.part -appl $IPS_APPL -goterms -iprlookup -dp -cpu $T -T $TMPD > ips.log 2>&1
  mv ips.tsv.part ips.tsv; rm -rf $TMPD
fi
for f in diamond.tsv pfam.domtbl ips.tsv; do [[ -s $f.gz && ! -s $f ]] && gunzip -k $f.gz; done
run claims python3 $BGC build-claims --obo $REF/go.obo --diamond diamond.tsv --sprot-go $REF/sprot_exp_go.tsv \
    --domtbl pfam.domtbl --pfam2go $REF/pfam2go.txt --ips ips.tsv --out claims.jsonl
# EC claims for Stage 2B and the gap prior: GO terms with a unique complete EC (ec2go) and Swiss-Prot EC transfer
run ec_claims python3 $CODE/add_ec_claims.py --claims claims.jsonl --diamond diamond.tsv --ec2go $REF/ec2go.txt --sprot-ec $REF/sprot_ec.tsv
NCL=$(wc -l < claims.jsonl)
gzip -f diamond.tsv pfam.domtbl ips.tsv
{ diamond version; hmmsearch -h | sed -n 2p; interproscan.sh --version 2>/dev/null | head -2; } > tool_versions.txt 2>&1 || true
for f in times/*.time; do printf '%s\t%s\t%s\n' $(basename $f .time) \
  "$(grep 'Elapsed (wall' $f | awk '{print $NF}')" "$(grep 'Maximum resident' $f | awk '{print $NF}')"; done > timing.tsv
printf '%s\t%s\t%s\t%s\n' "$C" "$N" "$NCL" "$(date -Is)" > DONE
echo "DONE $C representatives=$N claims=$NCL $(date -Is)"
