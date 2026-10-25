#!/bin/bash
# Macro-evolution of bacterial function from the GSPA annotation of AllTheBacteria (Stage D outputs).
#   bash run.sh deploy            # copy code and the gspa-core essential-function profiles to ibex
#   bash run.sh extract [N]       # ibex: N looping jobs (default 3) over lists/batches.txt (per-batch presence files)
#   bash run.sh list              # ibex: rebuild lists/batches.txt from the Stage D partials
#   bash run.sh profiles [STEPS]  # ibex job: 02 species profiles, 03 tree + parsimony, 04 tables (default "02 03 04")
#   bash run.sh status            # ibex: job states and counts
#   bash run.sh fetch             # local: results/, tree summaries and species table into benchmark/macroevolution/results/
set -euo pipefail
M=/encrypted/e3008/gspa/macroevo
R=/encrypted/e3008/allthebacteria/repro_2026-10-19
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$HERE/../.." && pwd)
case "${1:-help}" in
  deploy)
    ssh ibex "mkdir -p $M/code $M/ref $M/logs $M/presence $M/lists $M/results"
    scp "$HERE"/*.py "$HERE"/*.sbatch ibex:$M/code/
    scp -r "$REPO/gspa-core/src/main/resources/essential-functions" ibex:$M/ref/ ;;
  list)
    ssh ibex "for b in \$(ls $R/partials_stageD/*.genomes.tsv | xargs -n1 basename | sed s/.genomes.tsv//); do \
      if [[ \$b == batch_* ]]; then echo -e \"wave1\t\$b\"; else echo -e \"wave2\t\$b\"; fi; done > $M/lists/batches.txt; wc -l $M/lists/batches.txt" ;;
  extract)
    n=${2:-3}
    ssh ibex "for o in \$(seq 0 $((n-1))); do sbatch --export=ALL,BATCH_LIST=$M/lists/batches.txt,STRIDE=$n,OFFSET=\$o $M/code/01_extract_loop.sbatch; done" ;;
  profiles)
    ssh ibex "sbatch --export=ALL,STEPS='${2:-02 03 04}' $M/code/02_profiles.sbatch" ;;
  status)
    ssh ibex "squeue -u hohndor -h -o '%i %j %T %M' | grep macroevo; echo presence: \$(ls $M/presence/*.npz 2>/dev/null | wc -l) of \$(wc -l < $M/lists/batches.txt); ls $M/results 2>/dev/null | head -20" ;;
  fetch)
    mkdir -p "$HERE/results"
    scp ibex:$M/results/\*.tsv "$HERE/results/"
    scp ibex:$M/profiles/species.tsv ibex:$M/tree/tips.tsv ibex:$M/tree/nodes.tsv ibex:$M/tree/unmapped_species.tsv ibex:$M/tree/tree_pruned.nwk "$HERE/results/"
    scp ibex:$M/tree/parsimony_\*.tsv "$HERE/results/" ;;
  *) sed -n 2,9p "$0" ;;
esac
