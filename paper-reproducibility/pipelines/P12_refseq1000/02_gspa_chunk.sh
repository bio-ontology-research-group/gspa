#!/bin/bash
# P12 step 02 wrapper: run 02_gspa_genome.sh for CHUNK consecutive genomes
# (array task K covers genomes.txt lines K*CHUNK+1 .. (K+1)*CHUNK).  Needed
# because the IBEX per-user submit limit (2,000 queued tasks) does not allow
# a 1,000-task array next to other work.  Failures of one genome do not stop
# the chunk; re-submitting skips genomes with work/<g>/DONE.
R=${R:-/ibex/scratch/projects/c2014/hohndor/refseq1000_repro_2026-10-19}
CHUNK=${CHUNK:-10}
K=${SLURM_ARRAY_TASK_ID:-0}
rc=0
for g in $(sed -n "$(( K*CHUNK + 1 )),$(( (K+1)*CHUNK ))p" $R/genomes.txt); do
  bash $R/scripts/02_gspa_genome.sh "$g" || { echo "[FAIL] $g"; rc=1; }
done
exit $rc
