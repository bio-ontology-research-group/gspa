#!/bin/bash
# P13 — GSPA BGC benchmark on antiSMASH 8 regions of the 108 METEOR panel genomes (ibex, SLURM).
#
# Run from the laptop/repo root:
#   bash paper-reproducibility/pipelines/P13_bgc_antismash/run.sh sync      # rsync code to ibex
# then on ibex (sbatch only; nothing heavy runs on the login node):
#   bash $OUT/code/run.sh ref         # references (job 00)
#   bash $OUT/code/run.sh test        # 3 genomes, debug partition, IPS on 200 proteins, all 3 configs + eval
#   bash $OUT/code/run.sh all         # 108-genome array -> 3 integrate jobs -> evaluate (dependency chain)
#   bash $OUT/code/run.sh fetch       # (laptop) copy result TSVs back into results/
set -euo pipefail
export OUT=${OUT:-/ibex/scratch/projects/c2014/hohndor/bgc_repro_2026-10-19}
HERE=$(cd "$(dirname "$0")" && pwd)
J=$OUT/code/jobs
SB="sbatch --parsable --account=pi-hohndor --export=ALL"
TEST_GENOMES="GCA_003598085.1 GCA_022305115.1 GCA_046479585.1"
cfgs="stage0_noisy_or operon_only gspa_A9"

case "${1:-}" in
  sync)
    REPO=$(cd "$HERE/../../.." && pwd)
    ssh ibex "mkdir -p $OUT/code/jobs"
    rsync -a "$HERE/run.sh" "$HERE/README.md" "$HERE/OUTPUTS.md" ibex:$OUT/code/
    rsync -a "$HERE/jobs/" ibex:$OUT/code/jobs/
    rsync -a "$REPO/benchmark/unified/pillar2_fba_lro_price149/bgc_benchmark.py" "$REPO/benchmark/make_operons.py" ibex:$OUT/code/
    ;;
  ref)
    mkdir -p $OUT/logs
    $SB -o $OUT/logs/ref-%j.out $J/00_ref.sbatch
    ;;
  test)
    mkdir -p $OUT/logs $OUT/test
    L=$OUT/test/genomes_test.tsv
    printf 'genome\ttaxid\torganism\n' > $L
    for g in $TEST_GENOMES; do printf '%s\t-\t-\n' $g >> $L; done
    DEPA=${DEP:+-d afterok:$DEP}   # optional: wait for the ref job
    a=$($SB $DEPA -p debug -t 02:00:00 --array=0-2 -o $OUT/logs/test-genome-%A_%a.out \
         --export=ALL,GENOME_LIST=$L,IPS_MAX_SEQS=200 $J/01_genome.sbatch)
    ids=""
    for c in $cfgs; do
      i=$($SB -p debug -t 01:00:00 --mem=64G -d afterok:$a -o $OUT/logs/test-integ-$c-%j.out \
          --export=ALL,GENOME_LIST=$L,CFG=$c,MANIFEST_DIR=$OUT/test/manifests,INTEG_DIR=$OUT/test/integrated $J/02_integrate.sbatch)
      ids="$ids:$i"
    done
    e=$($SB -p debug -t 00:30:00 -d afterok$ids -o $OUT/logs/test-eval-%j.out \
        --export=ALL,GENOME_LIST=$L,INTEG_DIR=$OUT/test/integrated,TAG=test,NBOOT=200 $J/03_evaluate.sbatch)
    echo "test jobs: genome=$a integrate=${ids#:} eval=$e"
    ;;
  all)
    mkdir -p $OUT/logs
    K=${PER_TASK:-4}     # genomes per array task (ibex per-user submit limit is 2000 queued jobs)
    ng=$(( $(wc -l < $OUT/genomes.tsv) - 1 ))
    n=$(( (ng + K - 1) / K - 1 ))
    a=$($SB --array=0-$n -o $OUT/logs/genome-%A_%a.out --export=ALL,PER_TASK=$K $J/01_genome.sbatch)
    ids=""
    for c in $cfgs; do
      i=$($SB -d afterany:$a -o $OUT/logs/integ-$c-%j.out --export=ALL,CFG=$c $J/02_integrate.sbatch)
      ids="$ids:$i"
    done
    e=$($SB -d afterok$ids -o $OUT/logs/eval-%j.out $J/03_evaluate.sbatch)
    echo "jobs: genome-array=$a integrate=${ids#:} eval=$e"
    ;;
  integrate)   # re-run only integration + evaluation (e.g. after fixing a failed genome)
    ids=""
    for c in $cfgs; do
      i=$($SB -o $OUT/logs/integ-$c-%j.out --export=ALL,CFG=$c $J/02_integrate.sbatch); ids="$ids:$i"
    done
    e=$($SB -d afterok$ids -o $OUT/logs/eval-%j.out $J/03_evaluate.sbatch)
    echo "jobs: integrate=${ids#:} eval=$e"
    ;;
  fetch)
    REPO=$(cd "$HERE/../../.." && pwd)
    rsync -a ibex:$OUT/results/pillar2_3_antismash_bgc_specificity_v2.tsv "$REPO/results/"
    mkdir -p "$HERE/outputs"
    rsync -a ibex:$OUT/results/{bgc_per_genome_counts_v2.tsv,bgc_regions_all_v2.tsv,bgc_timing_v2.tsv} ibex:$OUT/ref/{sha256.txt,build_ref.log,class_terms_verified.tsv} ibex:$OUT/genomes.tsv "$HERE/outputs/"
    ;;
  *) echo "usage: $0 {sync|ref|test|all|integrate|fetch}"; exit 1;;
esac
