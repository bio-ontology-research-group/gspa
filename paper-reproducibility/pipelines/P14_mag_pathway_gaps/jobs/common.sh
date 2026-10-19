# Shared paths for P14 (sourced by every job script). Edit OUT to relocate.
OUT=${OUT:-/ibex/scratch/projects/c2014/hohndor/mags_repro_2026-10-19}
CODE=$OUT/code
REF=$OUT/ref
WORK=$OUT/work
RES=$OUT/results
LOGS=$OUT/logs
# Read-only inputs: Empty Quarter MAG tree (rob/metagenomics) and shared references
M=/ibex/scratch/projects/c2014/rob/metagenomics
HQ_LIST=$M/05_metabolism/hq_mag_list.tsv
GAPSEQ_OUT=$M/05_metabolism/gapseq_out
CHECKM2=$M/results/checkm2_single_sample_full.tsv
GTDB_BAC=$M/02_taxonomy/gtdbtk_out/gtdbtk.bac120.summary.tsv
GTDB_AR=$M/02_taxonomy/gtdbtk_out/gtdbtk.ar53.summary.tsv
GAPSEQ_DAT=/ibex/scratch/projects/c2014/rob/miniconda3/envs/gapseq-env/share/gapseq/dat   # gapseq 1.4.0 (same env as the GEMs)
P13_REF=/ibex/scratch/projects/c2014/hohndor/bgc_repro_2026-10-19/ref   # experimental-GO Swiss-Prot DIAMOND DB (P13 00_ref)
ATB_REF=/ibex/scratch/projects/c2014/hohndor/allthebacteria/ref          # go-dag, pfam2go, KEGG, ec2go
PFAM_HMM=/ibex/scratch/projects/c2014/rob/gspa/databases/pfam/Pfam-A.hmm  # same file as P13
JAR=/ibex/scratch/projects/c2014/hohndor/allthebacteria/bin/gspa-unified.jar   # v2.1.0
JAVA=/ibex/scratch/projects/c2014/rob/jdk-21/bin/java
PY=$CODE/mag_pathway_gaps.py
BGC=$CODE/bgc_benchmark.py   # build-claims only (copy of benchmark/unified/pillar2_fba_lro_price149/bgc_benchmark.py)
