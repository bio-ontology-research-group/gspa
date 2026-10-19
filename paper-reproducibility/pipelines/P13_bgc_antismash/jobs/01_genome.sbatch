#!/bin/bash
#SBATCH --job-name=p13-genome
#SBATCH --partition=batch
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
# One genome per array task: antiSMASH gbk -> proteins/GFF/labels -> operons -> DIAMOND + HMMER + InterProScan -> claims.
# Env: GENOME_LIST (default $OUT/genomes.tsv), IPS_MAX_SEQS (0 = all; >0 = test subset), SKIP_IPS=1 to drop channel 3.
set -euo pipefail
OUT=${OUT:-/ibex/scratch/projects/c2014/hohndor/bgc_repro_2026-10-19}
source $OUT/code/jobs/common.sh
module load python/3.11.0 diamond/2.1.16 hmmer/3.3.2 >/dev/null 2>&1
LIST=${GENOME_LIST:-$OUT/genomes.tsv}
# PER_TASK=K: array task i processes genome indices i*K .. i*K+K-1 sequentially
# (keeps the number of queued jobs below the per-user submit limit).
K=${PER_TASK:-1}
if [[ $K -gt 1 && -z "${P13_INNER:-}" ]]; then
  rc=0
  for j in $(seq 0 $((K-1))); do
    P13_INNER=1 SLURM_ARRAY_TASK_ID=$((SLURM_ARRAY_TASK_ID*K+j)) bash "$0" || { echo "FAILED index $((SLURM_ARRAY_TASK_ID*K+j))"; rc=1; }
  done
  exit $rc
fi
G=$(awk -F'\t' -v i=$SLURM_ARRAY_TASK_ID 'NR>1 && NR-2==i{print $1}' $LIST)
[[ -n "$G" ]] || { echo "no genome for index $SLURM_ARRAY_TASK_ID"; exit 0; }
C=${SLURM_CPUS_PER_TASK:-8}
D=$WORK/$G; mkdir -p $D; cd $D
echo "== $G  $(date -Is)  host=$(hostname)"
t0=$(date +%s)
python3 $BGC prep-genome --gbk $ANTISMASH/$G/$G.gbk --genome-id $G --out-dir $D
python3 $CODE/make_operons.py $G.gff ${G}_operons.tsv 300
t1=$(date +%s)
[[ -s ${G}_diamond.tsv ]] || { diamond blastp --query $G.faa --db $REF/sprot_exp.dmnd --threads $C \
    --outfmt 6 --evalue 1e-5 --max-target-seqs 5 --quiet --out ${G}_diamond.tsv.part \
    && mv ${G}_diamond.tsv.part ${G}_diamond.tsv; }
t2=$(date +%s)
[[ -s ${G}_pfam.domtbl ]] || { hmmsearch --domtblout ${G}_pfam.domtbl.part --noali --cpu $C --cut_ga $PFAM_HMM $G.faa > /dev/null \
    && mv ${G}_pfam.domtbl.part ${G}_pfam.domtbl; }
t3=$(date +%s)
IPS_ARG=""
if [[ "${SKIP_IPS:-0}" != "1" ]]; then
  if [[ ! -s ${G}_ips.tsv ]]; then
    module load interproscan/5.61-93.0 >/dev/null 2>&1
    TMPD=${TMPDIR:-/tmp}/ips_${G}_$$; mkdir -p $TMPD
    if [[ "${IPS_MAX_SEQS:-0}" -gt 0 ]]; then
      awk -v n=$IPS_MAX_SEQS '/^>/{k++} k<=n' $G.faa | sed '/^>/!s/\*//g' > $TMPD/in.faa
    else
      sed '/^>/!s/\*//g' $G.faa > $TMPD/in.faa
    fi
    interproscan.sh -i $TMPD/in.faa -f tsv -o ${G}_ips.tsv.part -appl $IPS_APPL -goterms -iprlookup -dp \
        -cpu $C -T $TMPD > ${G}_ips.log 2>&1
    mv ${G}_ips.tsv.part ${G}_ips.tsv
    rm -rf $TMPD
  fi
  IPS_ARG="--ips ${G}_ips.tsv"
fi
t4=$(date +%s)
python3 $BGC build-claims --obo $REF/go.obo --diamond ${G}_diamond.tsv --sprot-go $REF/sprot_exp_go.tsv \
    --domtbl ${G}_pfam.domtbl --pfam2go $REF/pfam2go.txt $IPS_ARG --out ${G}_claims.jsonl
t5=$(date +%s)
nprot=$(grep -c '^>' $G.faa)
printf 'genome\tn_proteins\tprep_s\tdiamond_s\thmmer_s\tips_s\tclaims_s\tcpus\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  $G $nprot $((t1-t0)) $((t2-t1)) $((t3-t2)) $((t4-t3)) $((t5-t4)) $C > ${G}_timing.tsv
cat ${G}_timing.tsv
echo "DONE $G $(date -Is)"
