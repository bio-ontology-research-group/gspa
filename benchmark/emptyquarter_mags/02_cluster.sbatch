#!/bin/bash
#SBATCH --job-name=eq-cluster
#SBATCH --partition=batch
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=180G
# Pool the proteomes of every bin in $BINS (ids <mag>|<protein_id>), cluster them with mmseqs linclust at
# 90 % identity and 90 % bidirectional coverage, split the representatives into chunks of $CHUNK_SIZE,
# write per-bin member maps, and submit the downstream jobs: the chunk channel array, the scatter job
# (afterok the array) and the integration array (afterok scatter, aftercorr the gapsmith array $GAPSMITH_JOB).
# Idempotent per stage (pool, cluster, chunks); re-run with RESUBMIT_ONLY=1 to only resubmit downstream jobs.
set -euo pipefail
# start from a clean module environment: the submitting shell may export LD_LIBRARY_PATH/PYTHONPATH of
# loaded modules (python/gcc), which breaks the conda tools and the static-linked gapsmith binary
module purge >/dev/null 2>&1 || true; unset LD_LIBRARY_PATH PYTHONPATH PYTHONHOME
source ${OUT:-/encrypted/e3008/gspa/eq_mags}/code/common.sh
T=${SLURM_CPUS_PER_TASK:-32}
mkdir -p $CL $CHUNKS $CBS $LOGS; cd $CL
NB=$(( $(wc -l < $BINS) - 1 ))
echo "== cluster step: $NB bins, chunk size $CHUNK_SIZE, $(date -Is) host=$(hostname) cpus=$T"

if [[ "${RESUBMIT_ONLY:-0}" != "1" ]]; then
  # --- 1. pool ---------------------------------------------------------------------------------
  if [[ ! -s pooled.faa.gz || ! -s pooled.counts.tsv ]]; then
    missing=0; : > pooled.counts.tsv.part
    tail -n +2 $BINS | cut -f1 | while read -r G; do
      f=$WORK/$G/annot/$G.faa
      if [[ -s $f ]]; then
        awk -v g="$G" '/^>/{sub(/^>/,""); split($0,a,/[ \t]/); print ">" g "|" a[1]; n++; next} {print} END{print g "\t" n > "/dev/stderr"}' $f 2>> pooled.counts.tsv.part
      else
        echo "MISSING proteome $G" >&2; missing=$((missing+1))
      fi
    done | gzip -1 > pooled.faa.gz.part
    mv pooled.faa.gz.part pooled.faa.gz; mv pooled.counts.tsv.part pooled.counts.tsv
    echo "pooled $(awk '{s+=$2} END{print s}' pooled.counts.tsv) proteins from $(wc -l < pooled.counts.tsv) bins"
    [[ $(wc -l < pooled.counts.tsv) -eq $NB ]] || { echo "proteome count $(wc -l < pooled.counts.tsv) != bins $NB"; exit 1; }
  fi
  # --- 2. linclust ---------------------------------------------------------------------------------
  if [[ ! -s clu_rep_seq.fasta || ! -s clu_cluster.tsv ]]; then
    rm -rf tmp_mmseqs; mkdir -p tmp_mmseqs
    /usr/bin/time -v -o linclust.time $MMSEQS easy-linclust pooled.faa.gz clu tmp_mmseqs \
      --min-seq-id 0.9 -c 0.9 --cov-mode 0 --threads $T --split-memory-limit 150G > linclust.log 2>&1
    rm -rf tmp_mmseqs clu_all_seqs.fasta
    echo "linclust: $(grep -c '^>' clu_rep_seq.fasta) representatives for $(wc -l < clu_cluster.tsv) proteins; $(grep 'Elapsed' linclust.time)"
  fi
  # --- 3. per-bin member maps (member -> representative) and chunks of representatives -----------
  if [[ ! -s members_by_bin.DONE ]]; then
    mkdir -p members_by_bin
    # clu_cluster.tsv: <rep>\t<member>; bin = text before the first '|' of the member id
    awk -F'\t' '{split($2,a,"|"); print a[1] "\t" $2 "\t" $1}' clu_cluster.tsv | sort -k1,1 -S 40G --parallel=$T \
      | awk -F'\t' 'BEGIN{OFS="\t"} {if ($1!=b) {if (b!="") close(f); b=$1; f="members_by_bin/" b ".tsv"} print $2, $3 > f}'
    [[ $(ls members_by_bin | wc -l) -eq $NB ]] || echo "WARNING: $(ls members_by_bin | wc -l) member maps for $NB bins"
    touch members_by_bin.DONE
  fi
  if [[ ! -s $LISTS/chunks.txt ]]; then
    rm -rf $CHUNKS/chunk_*; mkdir -p $CHUNKS
    awk -v n=$CHUNK_SIZE -v dir=$CHUNKS 'BEGIN{c=0;k=n} /^>/{if (k>=n) {if (f) close(f); f=sprintf("%s/chunk_%05d.faa", dir, c); c++; k=0} k++} {print > f}' clu_rep_seq.fasta
    ls $CHUNKS/chunk_*.faa | sort > $LISTS/chunks.txt.part; mv $LISTS/chunks.txt.part $LISTS/chunks.txt
    echo "$(wc -l < $LISTS/chunks.txt) chunks of <= $CHUNK_SIZE representatives"
  fi
fi

# --- 4. submit downstream jobs ---------------------------------------------------------------------
# sbatch is retried while the user's submitted-job cap (2,000, array tasks included) is reached
submit() { local id; for (( i = 0; i < 144; i++ )); do
  if id=$(sbatch --parsable "$@" 2> sbatch.err); then echo "$id"; return 0; fi
  grep -q "QOSMaxSubmitJobPerUserLimit\|job submit limit" sbatch.err || { cat sbatch.err >&2; return 1; }
  echo "submit limit reached, retrying in 10 min ($(date -Is))" >&2; sleep 600; done; return 1; }
NC=$(wc -l < $LISTS/chunks.txt)
EXP=ALL,OUT=$OUT,BINS=$BINS,WORK=$WORK,CL=$CL,CHUNKS=$CHUNKS,CBS=$CBS,LISTS=$LISTS,LOGS=$LOGS,GENOMES_PER_TASK=$GENOMES_PER_TASK
JC=$(submit --array=0-$(( NC - 1 ))%${THROTTLE_CH:-140} --export=$EXP -o $LOGS/chunk_%A_%a.out $CODE/03_chunk_channels.sbatch)
JS=$(submit --dependency=afterok:$JC -p batch -t 12:00:00 -c 4 --mem=64G -J eq-scatter -o $LOGS/scatter_%j.out --export=$EXP \
      --wrap="module load python/3.11.0; python3 $CODE/04_scatter.py --bins $BINS --cluster $CL --chunks $CHUNKS --out $CBS")
DEP="afterok:$JS"; [[ -n "${GAPSMITH_JOB:-}" ]] && DEP="$DEP,aftercorr:$GAPSMITH_JOB"
NT=$(n_tasks); JI=$(submit --array=0-$(( NT - 1 ))%${THROTTLE_IN:-60} --dependency=$DEP --export=$EXP -o $LOGS/in_%A_%a.out $CODE/06_integrate.sbatch)
printf 'chunks\t%s\nscatter\t%s\nintegrate\t%s\n' $JC $JS $JI >> ${JOBS_TSV:-$OUT/jobs.tsv}
echo "submitted: chunk array $JC ($NC chunks), scatter $JS, integrate array $JI ($NT tasks of $GENOMES_PER_TASK genomes; dependency $DEP)"
