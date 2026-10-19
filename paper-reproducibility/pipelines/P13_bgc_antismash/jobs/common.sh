# Shared paths for P13 (sourced by every job script). Edit OUT to relocate.
OUT=${OUT:-/ibex/scratch/projects/c2014/hohndor/bgc_repro_2026-10-19}
CODE=$OUT/code
REF=$OUT/ref
WORK=$OUT/work
INTEG=$OUT/integrated
RES=$OUT/results
LOGS=$OUT/logs
# Read-only inputs (colleague's METEOR tree and shared references)
ANTISMASH=/ibex/scratch/projects/c2014/kexin/funcarve/paperA_2026/results/antismash_gc
TAXTAB=/ibex/scratch/projects/c2014/kexin/funcarve/paperA_2026/results/gcf_taxonomy_gc.tsv
SPROT_DAT=/ibex/scratch/projects/c2014/rob/gspa-neural-deploy/reference/uniprot_sprot.dat.gz
GO_OBO_SRC=/ibex/scratch/projects/c2014/rob/gspa-neural-deploy/reference/go.obo   # releases/2026-03-25
PFAM_HMM=/ibex/scratch/projects/c2014/rob/gspa/databases/pfam/Pfam-A.hmm
ATB_REF=/ibex/scratch/projects/c2014/hohndor/allthebacteria/ref                   # go.owl (2026-03-25), pfam2go, KEGG, ec2go
JAR=/ibex/scratch/projects/c2014/hohndor/allthebacteria/bin/gspa-unified.jar      # v2.1.0
JAVA=/ibex/scratch/projects/c2014/rob/jdk-21/bin/java
IPS_APPL=${IPS_APPL:-Pfam,TIGRFAM,Hamap,PIRSF,SMART,CDD,ProSiteProfiles,ProSitePatterns,PRINTS,SFLD}   # Gene3D, SUPERFAMILY, PANTHER, FunFam omitted (cost; few InterPro2GO-bearing matches)
BGC=$CODE/bgc_benchmark.py
