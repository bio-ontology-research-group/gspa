#!/bin/bash
# P12 -- 1,000-genome RefSeq bacterial structural audit (Suppl. Table S8).
# Run from the laptop (repo root); every compute step is an sbatch job on
# IBEX (or unimatrix01 for the two GOA scans).  Never computes on a login node.
#
#   bash run.sh sync          copy this directory to IBEX ($R/scripts)
#   bash run.sh goa_subset    unimatrix01: GOA rows for reference_loo9 accessions
#   bash run.sh stage         relay reference_loo9.fasta, Pfam-A.hmm, GOA subset to IBEX
#   bash run.sh prep          IBEX: diamond makedb, hmmpress, RefSeq GFF3 -> operons, taxon map
#   bash run.sh matrix        IBEX: matrix reproduction of the six tool rows
#   bash run.sh test          IBEX debug: GSPA on the first 5 genomes
#   bash run.sh gspa [CONC]   IBEX: GSPA array over 1,000 genomes (default %60)
#   bash run.sh eval tools|gspa   IBEX: structural suite + gspa evaluate, 40 chunks x 25
#   bash run.sh truth         unimatrix01: GOA experimental annotations for the RefSeq proteins
#   bash run.sh aggregate     IBEX: per-genome TSVs + 8-row table; then `fetch`
#   bash run.sh fetch         copy aggregated outputs into the repo
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$HERE/../../.." && pwd)
R=/ibex/scratch/projects/c2014/hohndor/refseq1000_repro_2026-10-19
U=/data/hohndor/gspa/refseq1000_repro_2026-10-19
C=/ibex/scratch/projects/c2014/annot_evaluation
PY_U=/storage/miniforge3/envs/metagenomics/bin/python3
cmd=${1:-help}; shift || true
case $cmd in
  sync)
    cp "$REPO"/benchmark/unified/pillar1_5_6/p12_{structural_suite,matrix_rows}.py "$HERE"/
    ssh ibex "mkdir -p $R/{bin,ref,db,logs,scripts,work,results,eval}"
    rsync -a "$HERE"/ ibex:$R/scripts/ ;;
  goa_subset)
    ssh unimatrix01 "mkdir -p $U/{scripts,logs}"; scp -q "$HERE"/00_build_goa_subset.py unimatrix01:$U/scripts/
    ssh unimatrix01 "cd $U && sbatch --parsable -p debug -x node007,unimatrix -c 2 --mem=8G -J p12_goa -o $U/logs/goa_subset_%j.log \
      --wrap 'cd $U; $PY_U scripts/00_build_goa_subset.py /data/hohndor/gspa/proteomes/reference_loo9.fasta \
      /data/hohndor/gspa/benchmark/goa_uniprot_all.gaf.gz goa_reference_loo9_subset.gaf && sha256sum \
      /data/hohndor/gspa/proteomes/reference_loo9.fasta /data/hohndor/gspa/benchmark/goa_uniprot_all.gaf.gz \
      /storage/software/databases/hmmer/Pfam-A.hmm goa_reference_loo9_subset.gaf > input_sha256.txt'" ;;
  stage)
    ssh ibex "mkdir -p $R/{bin,ref,db}; cp -n /ibex/scratch/projects/c2014/hohndor/allthebacteria/bin/gspa-unified.jar $R/bin/; \
      cd /ibex/scratch/projects/c2014/hohndor/allthebacteria/ref && cp -rn go.owl go.obo go_aspect_map.tsv kegg_modules.tsv \
      kegg_pathways.tsv ec2go.txt pfam2go.txt $R/ref/; cp -n $C/1000files $R/genomes.txt"
    scp -3 -q unimatrix01:/data/hohndor/gspa/proteomes/reference_loo9.fasta ibex:$R/db/
    scp -3 -q unimatrix01:/storage/software/databases/hmmer/Pfam-A.hmm ibex:$R/db/
    scp -3 -q unimatrix01:$U/goa_reference_loo9_subset.gaf ibex:$R/db/ ;;
  prep)    ssh ibex "cd $R && sbatch --parsable -o $R/logs/%x_%j.log scripts/01_prep_ibex.sbatch" ;;
  matrix)  ssh ibex "cd $R && sbatch --parsable -J p12_matrix -p batch -c 1 --mem=8G -t 04:00:00 -o $R/logs/%x_%j.log \
             --wrap 'python3 $R/scripts/p12_matrix_rows.py --annot-eval $C --out-json $R/results/matrix_rows.json'" ;;
  test)    ssh ibex "cd $R && sbatch --parsable -J p12_test -p debug -c 4 --mem=20G -t 02:00:00 --array=0-4 \
             -o $R/logs/%x_%A_%a.log scripts/02_gspa_genome.sh" ;;
  gspa)    CONC=${1:-60}   # 100 tasks x 10 genomes
           ssh ibex "cd $R && sbatch --parsable -J p12_gspa -p batch -c 4 --mem=20G -t 2-00:00:00 --array=0-99%$CONC --export=ALL,CHUNK=10 \
             -o $R/logs/gspa/%x_%A_%a.log scripts/02_gspa_chunk.sh" ;;
  eval)    M=${1:-tools}
           ssh ibex "cd $R && sbatch --parsable -J p12_eval_$M -p batch -c 1 --mem=12G -t 02:00:00 --array=0-39 \
             --export=ALL,METHODS=$M,CHUNK=25,GAEF_METHODS=${GAEF_METHODS:-.} -o $R/logs/%x_%A_%a.log scripts/03_eval_chunk.sh" ;;
  truth)   scp -q "$HERE"/04_goa_exp_truth.py unimatrix01:$U/
           echo "needs $U/protein_ids.tsv (genome<TAB>WP id, from done_proteins/*_protein.faa headers) and propagated_taxon_ids.tsv"
           ssh unimatrix01 "cd $U && sbatch --parsable -p debug -x node007,unimatrix -c 2 --mem=24G -J p12_truth -o $U/logs/truth_%j.log \
             --wrap 'cd $U; $PY_U 04_goa_exp_truth.py protein_ids.tsv /data/hohndor/gspa/reference/gene_refseq_uniprotkb_collab.gz \
             /data/hohndor/gspa/benchmark/goa_uniprot_all.gaf.gz propagated_taxon_ids.tsv truth'" ;;
  aggregate)
           scp -3 -q unimatrix01:$U/truth/exp_truth_summary.tsv ibex:$R/results/ || true
           ssh ibex "cd $R && sbatch --parsable -J p12_agg -p batch -c 1 --mem=16G -t 02:00:00 -o $R/logs/%x_%j.log \
             --wrap 'python3 $R/scripts/05_aggregate.py --root $R --out-dir $R/results/aggregate --truth-summary $R/results/exp_truth_summary.tsv'" ;;
  fetch)   mkdir -p "$HERE"/results
           rsync -a ibex:$R/results/aggregate/ "$HERE"/results/
           rsync -a ibex:$R/input_sha256_ibex.txt ibex:$R/tool_versions.txt ibex:$R/taxon_map.tsv "$HERE"/results/ || true
           cp "$HERE"/results/table_s8_regenerated.tsv "$REPO"/results/pillar1_suite2_1000bacteria_sota.tsv ;;
  *) sed -n 2,20p "$0" ;;
esac
