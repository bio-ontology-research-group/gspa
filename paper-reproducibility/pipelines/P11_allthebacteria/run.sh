#!/bin/bash
# P11 -- AllTheBacteria pan-phylum analysis from the real GSPA integrator (JVM) output.
# Host: KAUST Ibex (Slurm, account pi-hohndor). Never run compute on the login node:
# every step below is an sbatch submission; the login node only copies files.
#
# Usage (from a local checkout; ssh alias `ibex`):
#   bash run.sh deploy                 # copy Stage D scripts + constraint tables to $R
#   bash run.sh forbidden              # (local, optional) re-dump the jar's forbidden-term tables
#   bash run.sh lists                  # (ibex job) inventory: completed batches, Java runtimes
#   bash run.sh test                   # (ibex, debug) byte-compat test on batch_0001 (20 genomes)
#   bash run.sh stageD-wave1           # (ibex) Java integrate on all genomes of the 43 wave-1 batches
#   bash run.sh stageD-wave2 FROM TO   # (ibex) same for lines FROM..TO of the completed wave-2 list
#   bash run.sh summarize-existing     # (ibex) run (a): summarise Java outputs Stage C already wrote
#   bash run.sh combine PARTIALS TAG   # (ibex job) merge partials -> tables; TAG = file suffix after _java ("" for run (a), _legacy, _stageD)
#   bash run.sh fetch TAG              # copy tables back into results/allthebacteria/
#
# Stages A-C are the original pipeline (benchmark/allthebacteria/0{3,4}_*.sbatch, wave-1 job
# 53958xxx / wave-2 array 53958302); this script only adds Stage D and the aggregation.
set -euo pipefail
S=/ibex/scratch/projects/c2014/hohndor/allthebacteria
R=$S/repro_2026-10-19
PY=/ibex/scratch/projects/c2014/rob/envs/gspa-neural/bin/python
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$HERE/../../.." && pwd)
B=$REPO/benchmark/allthebacteria

mode=${1:-help}
case "$mode" in
  deploy)
    ssh ibex "mkdir -p $R/scripts $R/ref $R/logs $R/batch_lists $R/partials_stageC $R/partials_stageC_legacy $R/partials_stageD"
    scp "$B"/stageD_*.py "$B"/stageD_*.sbatch "$B"/aggregate_atb_java.py "$B"/dump_forbidden_terms.groovy ibex:$R/scripts/
    ;;
  forbidden)
    # Uses the same jar (sha256 e415d9c6...72a7) locally; input = all taxids in the batch manifests.
    out=${2:-$PWD/forbidden_by_taxon.tsv}
    ssh ibex "cat $S/wave1/batches/*.tsv $S/wave2/manifests/*.tsv" | awk -F'\t' '$1!="sample"{print $5}' | sort -u > taxa_all.txt
    java -cp "$REPO/gspa-unified.jar" groovy.ui.GroovyMain "$B/dump_forbidden_terms.groovy" "$out" $(cat taxa_all.txt)
    unzip -p "$REPO/gspa-unified.jar" complex-terms/complex_terms.tsv > complex_terms.tsv
    scp "$out" "${out%.tsv}.sets.tsv" complex_terms.tsv ibex:$R/ref/
    ;;
  lists)
    ssh ibex "sbatch -p batch -t 02:00:00 -c 2 --mem=8G -J atb_survey -o $R/logs/survey_%j.out \
      --wrap='$PY $R/scripts/stageD_survey_batches.py --out-dir $R/survey'"
    # batch lists (wave<TAB>batch) are then derived from $R/survey/batch_inventory.tsv:
    #   wave1.txt               all 43 wave-1 batches (stageC_done=1)
    #   wave2_done_<date>.txt   wave-2 batches with stageC/batch_stats.tsv
    #   existing_current.txt    batches whose Stage C Java run used the final config (java_batch_stats.tsv present)
    #   existing_legacy.txt     w2_0013 w2_0018 w2_0019 w2_0020 w2_0021 w2_0024 (older Stage C: max_iter 5, pairwise operons)
    ;;
  test)
    ssh ibex "sbatch --export=ALL,WAVE=wave1,BATCH=batch_0001 $R/scripts/stageD_bytecompat_test.sbatch"
    ;;
  stageD-wave1)
    ssh ibex "sbatch --array=1-43%20 --export=ALL,BATCH_LIST=$R/batch_lists/wave1.txt $R/scripts/stageD_java_integrate_all.sbatch"
    ;;
  stageD-wave2)
    from=${2:?FROM}; to=${3:?TO}; list=${4:-$R/batch_lists/wave2_done_2026-10-19a.txt}
    ssh ibex "sbatch --array=${from}-${to}%12 --cpus-per-task=16 --mem=64G --time=1-00:00:00 \
      --export=ALL,BATCH_LIST=$list $R/scripts/stageD_java_integrate_all.sbatch"
    ;;
  summarize-existing)
    ssh ibex "n=\$(wc -l < $R/batch_lists/existing_current.txt); \
      sbatch --array=1-\$n%30 --export=ALL,BATCH_LIST=$R/batch_lists/existing_current.txt $R/scripts/stageD_summarize_existing.sbatch; \
      sbatch --array=1-6 --export=ALL,BATCH_LIST=$R/batch_lists/existing_legacy.txt,PARTIALS=$R/partials_stageC_legacy $R/scripts/stageD_summarize_existing.sbatch"
    ;;
  combine)
    parts=${2:?PARTIALS dir, e.g. $R/partials_stageD}; tag=${3-}
    ssh ibex "sbatch -p batch -t 02:00:00 -c 2 --mem=16G -J atb_combine -o $R/logs/combine_%j.out \
      --wrap='$PY $R/scripts/aggregate_atb_java.py combine --partials $parts --out-dir $R/results$tag --suffix _java$tag --summary-tag $tag'"
    ;;
  fetch)
    tag=${2-}
    mkdir -p "$REPO/results/allthebacteria"
    scp "ibex:$R/results$tag/*" "$REPO/results/allthebacteria/"
    ;;
  *)
    sed -n '2,20p' "$0"
    ;;
esac
