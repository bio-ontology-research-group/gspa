#!/bin/bash
# Stage D follow-up: when a Stage D array has ended, resubmit every wave-2 batch that still has no partial.
# Batches whose task hit the 24 h wall (the 4,000-genome batches run 20-29 h at 96 genomes per 40 min on
# 16 CPUs) keep the genomes they integrated: stageD_java_integrate_batch.py skips genomes whose two .gz
# outputs exist, so the resubmitted task only runs the remainder and then writes the partial.
#
#   sbatch --dependency=afterany:<STAGE_D_ARRAY> --export=ALL,CONC=60 stageD_autoresume.sbatch
#
# Repeats itself: the new array gets its own afterany-chained autoresume, until no batch is left.
#SBATCH --job-name=gspa_atb_autoresume
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --output=/encrypted/e3008/allthebacteria/repro_2026-10-19/logs/autoresume_%j.out
set -euo pipefail
R=/encrypted/e3008/allthebacteria/repro_2026-10-19
: "${CONC:=60}"
: "${SOURCE_LIST:=$R/batch_lists/wave2_done_2026-10-19a.txt}"
list=$R/batch_lists/wave2_remaining_$(date +%F_%H%M).txt
awk -F'\t' -v P=$R/partials_stageD '{ if (system("test -f " P "/" $2 ".genomes.tsv") != 0) print }' "$SOURCE_LIST" > "$list"
n=$(wc -l < "$list")
echo "$(date -Is) $n batches without a partial -> $list"
if [[ $n -eq 0 ]]; then echo "Stage D complete: every batch of $SOURCE_LIST has a partial"; exit 0; fi
# The per-user QOS allows 2,000 submitted jobs and every array task counts, so submit only what fits
# (keeping a margin for other submissions); the chained autoresume picks up the rest when this array ends.
free=$(( 2000 - $(squeue -u hohndor -h -r | wc -l) - 60 ))
if [[ $free -lt 10 ]]; then
  echo "only $free job slots free; retrying in 2 h"
  sbatch --parsable --begin=now+2hours --export=ALL,CONC=$CONC,SOURCE_LIST=$SOURCE_LIST $R/scripts/stageD_autoresume.sbatch
  exit 0
fi
m=$(( n < free ? n : free ))
jid=$(sbatch --parsable --array=1-$m%$CONC --cpus-per-task=16 --mem=64G --time=1-00:00:00 \
  --export=ALL,BATCH_LIST=$list $R/scripts/stageD_java_integrate_all.sbatch)
echo "resubmitted $m of $n batches as array $jid"
next=$(sbatch --parsable --dependency=afterany:$jid --export=ALL,CONC=$CONC,SOURCE_LIST=$SOURCE_LIST $R/scripts/stageD_autoresume.sbatch)
echo "next autoresume: $next (afterany:$jid)"
printf '%s\t%s\t%s/%s\t%s\n' "$(date -Is)" "$jid" "$m" "$n" "$next" >> $R/batch_lists/autoresume_chain.tsv
