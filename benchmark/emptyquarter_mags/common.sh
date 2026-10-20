# Shared paths for the Empty Quarter HQ/MQ bin annotation (sourced by every job script).
# Work area on the e3008 project space; edit OUT to relocate.
OUT=${OUT:-/encrypted/e3008/gspa/eq_mags}
CODE=$OUT/code
REF=$OUT/ref
WORK=${WORK:-$OUT/work}
LISTS=${LISTS:-$OUT/lists}
LOGS=${LOGS:-$OUT/logs}
RES=${RES:-$OUT/results}
BINS=${BINS:-$LISTS/bins.tsv}            # written by 00_select.py (one HQ/MQ bin per row)
CL=${CL:-$OUT/cluster}                   # pooled proteomes, linclust output, per-bin member maps
CHUNKS=${CHUNKS:-$OUT/chunks}            # representative chunks and their channel results
CBS=${CBS:-$OUT/claims_by_sample}        # expanded claims, one file per sample
CHUNK_SIZE=${CHUNK_SIZE:-50000}          # representatives per channel chunk (about 6 h on 8 CPUs; keeps the array under 500 tasks)
GENOMES_PER_TASK=${GENOMES_PER_TASK:-50} # rows of bins.tsv per array task of 01/05/06 (ibex MaxSubmitJobs 2,000 per user)
MMSEQS=/ibex/scratch/projects/c2014/rob/miniconda3/envs/mmseqs2/bin/mmseqs   # 18.8cc5c

# Read-only inputs: Empty Quarter metagenomics tree (rob/metagenomics) and shared references
M=/ibex/scratch/projects/c2014/rob/metagenomics
CHECKM2=$M/results/checkm2_single_sample_full.tsv
GTDB_BAC=$M/02_taxonomy/gtdbtk_out/gtdbtk.bac120.summary.tsv
GTDB_AR=$M/02_taxonomy/gtdbtk_out/gtdbtk.ar53.summary.tsv
GAPSEQ_OUT_HQ=$M/05_metabolism/gapseq_out          # gapseq 1.4.0 models of the 3,835 HQ MAGs (reused)
# Culture genomes: the RAK-PR catalogue (Alejandra's pangenome reference; enrichment bins, isolates and the
# thermophilic cultures, which the catalogue files as batch "soil"), their assembly registry and Bakta annotations
RAKPR=/ibex/scratch/projects/c2014/alelopezv/eqmgc
RAKPR_CATALOGUE=$RAKPR/catalogue.tsv               # 722 HQ/MQ genomes (genome_id, source_sample, batch, CheckM2, BUSCO, GTDB, tier)
RAKPR_INPUTS=$RAKPR/inputs.tsv                     # genome_id -> assembly path (uncompressed .fa)
RAKPR_BAKTA=$RAKPR/results/annotation/bakta        # <genome_id>/<genome_id>.{faa,gff3,fna} for 646 of the 722
P14_REF=/ibex/scratch/projects/c2014/hohndor/mags_repro_2026-10-19/ref   # go.owl, ec2go, KEGG, pfam2go, sprot_exp.dmnd, SEED tables, essential profiles
FBA_DATA_SRC=/encrypted/e3008/gspa/meteor108_repro_2026-10-20/data        # universal_metadata.json + fba_matrices + bounds + medium + seed maps
PFAM_HMM=/ibex/scratch/projects/c2014/rob/gspa/databases/pfam/Pfam-A.hmm
SPROT_DAT=/ibex/scratch/projects/c2014/rob/gspa-neural-deploy/reference/uniprot_sprot.dat.gz   # source of ref/sprot_ec.tsv (build_sprot_ec.py)

# Tools
JAVA=/ibex/scratch/projects/c2014/rob/jdk-21/bin/java
JAR=$OUT/bin/gspa-unified.jar                       # v2.1.0, copied at `run.sh stage`
PYRODIGAL=/ibex/scratch/projects/c2014/rob/miniconda3/envs/pyrodigal/bin/pyrodigal
GAPSEQ=/ibex/scratch/projects/c2014/rob/miniconda3/envs/gapseq-env/bin/gapseq      # fallback only (02_gapseq.sbatch)
GAPSMITH=$M/tools/gapsmith/bin/gapsmith                                               # gapsmith 0.1.0 (Rust gapseq)
GAPSMITH_DATA=$M/tools/gapsmith/data_merged
GAPSMITH_MODELS=$M/06_darkmatter/gapsmith_models                                      # 5,438 dRep-representative models (<sample>__{single,self}_SemiBin_N)
DIAMOND_FOR_GAPSMITH=/ibex/scratch/projects/c2014/rob/miniconda3/envs/bakta/bin/diamond   # 2.1.24, the aligner the earlier gapsmith run used (no module, no gcc libs)
FBA_PY=/ibex/scratch/projects/c2014/hohndor/meteor108_repro_2026-10-19/venv/bin/python   # numpy/scipy for the Stage 2B sidecar
IPS_MODULE=${IPS_MODULE:-interproscan/5.61-93.0}
IPS_APPL=${IPS_APPL:-Pfam,TIGRFAM,Hamap,PIRSF,SMART,CDD,ProSiteProfiles,ProSitePatterns,PRINTS,SFLD}   # as P13

# GSPA code (copied at stage)
PREP=$CODE/eq_prep.py
P14PY=$CODE/mag_pathway_gaps.py
BGC=$CODE/bgc_benchmark.py        # build-claims
FBA_SIDECAR=$CODE/stoichiometric_locus_refiner.py

# One row of $BINS -> shell variables (G mag_id, SAMPLE, BIN, CLASS, DOMAIN, TAXON, GRAM, PROF, FA, BD, GD, DATASET)
bin_row() {   # usage: bin_row <0-based index>
  local i=$1 row hdr
  row=$(awk -F'\t' -v i="$i" 'NR>1 && NR-2==i' "$BINS"); [[ -n "$row" ]] || return 1
  hdr=$(head -1 "$BINS")
  col() { paste <(tr '\t' '\n' <<< "$hdr") <(tr '\t' '\n' <<< "$row") | awk -F'\t' -v c="$1" '$1==c{print $2}'; }
  G=$(col mag_id); SAMPLE=$(col sample); BIN=$(col bin); CLASS=$(col class); DOMAIN=$(col domain)
  TAXON=$(col taxon); GRAM=$(col gram); PROF=$(col essential_profile)
  FA=$(col mag_fasta); BD=$(col bakta_dir); GD=$(col gapseq_dir); DATASET=$(col dataset)
}
n_tasks() {   # number of array tasks for the rows of $BINS at GENOMES_PER_TASK rows per task
  local nb=$(( $(wc -l < "$BINS") - 1 )); echo $(( (nb + GENOMES_PER_TASK - 1) / GENOMES_PER_TASK ))
}
