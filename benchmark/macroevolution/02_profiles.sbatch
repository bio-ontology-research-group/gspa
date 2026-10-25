#!/bin/bash
# Species profiles (02), tree and parsimony (03), result tables (04) in one job after the extraction.
#SBATCH --job-name=gspa_macroevo_profiles
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=/encrypted/e3008/gspa/macroevo/logs/profiles_%j.out
set -euo pipefail
S=/ibex/scratch/projects/c2014/hohndor/allthebacteria
R=/encrypted/e3008/allthebacteria/repro_2026-10-19
M=/encrypted/e3008/gspa/macroevo
PY=/ibex/scratch/projects/c2014/rob/envs/gspa-neural/bin/python
GTDB=/ibex/reference/KSL/gtdb-tk/release220
: "${STEPS:=02 03 04}"; : "${SUF:=}"; : "${THR:=0.5}"   # SUF=_test writes to profiles_test/ tree_test/ results_test/; THR=0.3 SUF=_t30 for the sensitivity run
for step in $STEPS; do
  echo "[$(date '+%F %T')] step $step"
  case $step in
    02) $PY $M/code/02_species_profiles.py --presence $M/presence --manifests "$S/wave1/batches/batch_*.tsv" "$S/wave2/manifests/w2_*.tsv" \
          --insight3 $R/results_stageD/insight3_per_genome_java_stageD.tsv --out $M/profiles$SUF --threshold $THR --gtdb-metadata $M/ref/bac120_metadata_r220.tsv.gz --never-in-taxon $S/ref/never_in_taxon-propagated.tsv ;;
    03) $PY $M/code/03_tree.py --profiles $M/profiles$SUF --gtdb-metadata $M/ref/bac120_metadata_r220.tsv.gz \
          --tree $GTDB/pplacer/gtdb_r220_bac120.refpkg/gtdb_r220_bac120_decorated_unrooted.tree \
          --train-taxonomy $S/ref/train_taxonomy.tsv --taxdump-nodes $M/ref/taxdump/nodes.dmp --min-ref-proteins 1 --go-vocab $M/presence/go_vocab.txt --obo $S/ref/go.obo --out $M/tree$SUF ;;
    04) $PY $M/code/04_analyses.py --profiles $M/profiles$SUF --tree $M/tree$SUF --presence $M/presence \
          --never-in-taxon $S/ref/never_in_taxon-propagated.tsv --essential $M/ref/essential-functions/bacteria.tsv \
          --obo $S/ref/go.obo --out $M/results$SUF ;;
    05) $PY $M/code/05_rescue_floor.py --presence $M/presence --out $M/results$SUF ;;
  esac
done
echo "[$(date '+%F %T')] done"
