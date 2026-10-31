# AllTheBacteria pilot on Shaheen III (gate G0: core-hours per genome per step).
# Sourced by every job script and by driver.sh. Every path is a variable so the release candidate,
# the reference tree and the work area can move without editing the scripts.
#
# Port of the validated Empty Quarter pipeline (benchmark/emptyquarter_mags/, IBEX): same tools, versions,
# parameters and claim builders; reshaped for exclusive 192-core workq nodes and the Lustre rules
# (one directory per shard, no per-genome directories on /scratch, per-genome work in node-local tmpfs,
# outputs packed per shard).

set -a   # every variable below is exported (xargs/bash -c workers and the container need them)
# --- where things are ------------------------------------------------------------------------------
ATB=${ATB:-/project/k10018/allthebacteria}                         # archives/, pipeline/catalog, pipeline/wave2/manifests
ATB_CATALOG=${ATB_CATALOG:-$ATB/pipeline/catalog/wave2_full_2.4M_batches.tsv}
ATB_MANIFESTS=${ATB_MANIFESTS:-$ATB/pipeline/wave2/manifests}
ATB_ARCHIVES=${ATB_ARCHIVES:-$ATB/archives}
# GSPA code (python sidecars, claim builders, these scripts). Must be the tree the JAR was built from.
CODE=${CODE:-/project/k10018/code/gspa-v3}
HERE_CODE=${HERE_CODE:-$CODE/benchmark/allthebacteria/shaheen}
JAR=${JAR:-/project/k10018/gspa-v3/release-candidate-1/gspa-3.0.0.jar}
JAR_SHA256=${JAR_SHA256:-3fe6ba0a9de2f2cd09d79de05239b070332b16d7759455508891dab1fc1be064}
# Reference tree (read-only). Files the EQ run kept in $OUT/ref on IBEX; see the inventory in the report.
REF=${REF:-/project/k10018/gspa-v3/eq-ref-v1/ref}   # EQ ref dir from IBEX, hash-verified 2026-10-30
GO_OWL=${GO_OWL:-$ATB/pipeline/ref/go.owl}                         # sha256 = EQ/P14 go.owl
GO_OBO=${GO_OBO:-$REF/go.obo}                                       # EQ/P14 go.obo (ATB's go.obo differs)
GO_DAG=${GO_DAG:-$ATB/pipeline/ref/go-dag.tsv}
EC2GO=${EC2GO:-$ATB/pipeline/ref/ec2go.txt}
PFAM2GO=${PFAM2GO:-$ATB/pipeline/ref/pfam2go.txt}
KEGG_PATHWAYS=${KEGG_PATHWAYS:-$ATB/pipeline/ref/kegg_pathways.tsv}
KEGG_MODULES=${KEGG_MODULES:-$ATB/pipeline/ref/kegg_modules.tsv}
ESSENTIAL_DIR=${ESSENTIAL_DIR:-$ATB/pipeline/ref/essential-functions}
TAXON_LABELS=${TAXON_LABELS:-$ATB/pipeline/ref/taxon-constraints/ncbi-taxon-labels.tsv}
SEED_REF=${SEED_REF:-/project/k10018/gspa-v3/gapsmith-reference-bacteria-v1/reference}
SEED_REACTIONS=${SEED_REACTIONS:-$SEED_REF/seed_reactions_corrected.tsv}
DIFFUSION_METS=${DIFFUSION_METS:-$SEED_REF/diffusion_mets.tsv}
EC_ALIASES=${EC_ALIASES:-$SEED_REF/seed_Enzyme_Class_Reactions_Aliases_unique.tsv}
FBA_DATA=${FBA_DATA:-$CODE/benchmark/manifests/reference_data/psb_revision/code_snapshot/data}
SPROT_DMND=${SPROT_DMND:-$REF/sprot_exp.dmnd}   # staged from IBEX
SPROT_GO=${SPROT_GO:-$REF/sprot_exp_go.tsv}   # staged from IBEX
SPROT_EC=${SPROT_EC:-$REF/sprot_ec.tsv}   # staged from IBEX
PFAM_HMM=${PFAM_HMM:-/project/k10018/sw/db/pfam-ibex/pfam/Pfam-A.hmm}   # staged from IBEX
GAPSMITH_DATA=${GAPSMITH_DATA:-/project/k10018/gspa-v3/gapsmith-reference-bacteria-v1/reference}
# --- tools -------------------------------------------------------------------------------------------
TOOLS=${TOOLS:-/project/k10018/sw/envs/gspa-tools}                 # pyrodigal 3.7.1, mmseqs2 18.8cc5c, DIAMOND 2.1.16, HMMER 3.3.2
GAPSMITH=${GAPSMITH:-/project/k10018/gspa-v3/tools/gapsmith-linux} # sha256 8bb3cddd… (v3 cohort binary)
IPS_HOME=${IPS_HOME:-/project/k10018/sw/db/interproscan-5.61-93.0-ibex/interproscan-5.61-93.0}   # IBEX module copy without funfam/gene3d/panther/superfamily
IPS_APPL=${IPS_APPL:-Pfam,TIGRFAM,Hamap,PIRSF,SMART,CDD,ProSiteProfiles,ProSitePatterns,PRINTS,SFLD}   # EQ / P13 set
PY_HOST=${PY_HOST:-/project/k10018/sw/envs/gspa-v3/bin/python}     # stdlib-only helpers, run natively
# Container: Debian 13 sandbox + the shared package tree (JDK 21, Python 3.10 testenv with scipy/cobra/libsbml),
# the same recipe as benchmark/v3/shaheen/arm-pack.sbatch. The tree is reused, never extracted again here.
IMAGE=${IMAGE:-/project/k10018/sw/images/debian-trixie}
PKG_ROOT=${PKG_ROOT:-/scratch/project/k10018/gspa-v3/shared/arm-fs-d04c9e11ee64db47}
CTR_PY=/home/leechuck/tmp/gspa-v3-testenv/bin/python               # path inside the container
CTR_JAVA=/usr/lib/jvm/java-21-openjdk-amd64/bin/java
# --- work area (one directory per shard; see the file budget in the report) -------------------------
OUT=${OUT:-/scratch/project/k10018/gspa-v3/atb-pilot}
LISTS=$OUT/lists; SHARDS_DIR=$OUT/shards; RECEIPTS=$OUT/receipts; LOGS=$OUT/logs
LOCAL_ROOT=${LOCAL_ROOT:-/dev/shm}                                  # node-local tmpfs for per-genome work
CORES_PER_NODE=${CORES_PER_NODE:-192}                               # physical cores charged per workq node-hour
# --- pilot shape ---------------------------------------------------------------------------------------
SEED=${SEED:-20261030}
N_SHARDS=${N_SHARDS:-20}
GENOMES_PER_SHARD=${GENOMES_PER_SHARD:-500}
CHUNK_SIZE=${CHUNK_SIZE:-50000}                                     # representatives per channel chunk (plan 3.2)
INODE_LIMIT=${INODE_LIMIT:-990000}                                  # refuse to start above this many files on /scratch

set +a

nlines() { if [[ -s $1 ]]; then wc -l < "$1"; else echo 0; fi; }   # lines of a file that may not exist
shard_id() { printf 's%03d' "$1"; }

# Node-local work dir, removed on exit. Per-genome files live only here.
make_local() {
  LOCAL=$LOCAL_ROOT/atbp-${SLURM_JOB_ID:-$$}
  mkdir -p "$LOCAL/tmp" "$LOCAL/times"
  PHASES=$LOCAL/phases.tsv; VERSIONS=$LOCAL/versions.tsv; : > "$PHASES"; : > "$VERSIONS"
  export LOCAL PHASES VERSIONS
}
# Job wall clock starts here (receipts charge nodes x cores x wall).
JOB_START=${JOB_START:-$(date +%s.%N)}
# version <name> <command...>: first line of the command's output into the versions file
version() { local n=$1; shift; printf '%s\t%s\n' "$n" "$("$@" 2>&1 | head -1 | tr '\t' ' ')" >> "$VERSIONS"; }
# Shards of this job: SHARDS=s000:s001:...
shard_list() { tr ':' ' ' <<< "${SHARDS:?SHARDS=s000:s001:... required}"; }
# shard_cols <shard> <column>...: the named columns of the shard's rows of the pilot list, in the order given
shard_cols() {
  local s=$1; shift
  awk -F'\t' -v s="$s" -v want="$*" 'BEGIN{OFS="\t"; n=split(want, w, " ")}
    NR==1{for(i=1;i<=NF;i++) c[$i]=i; for(j=1;j<=n;j++) if (!(w[j] in c)) {print "no column " w[j] > "/dev/stderr"; exit 1}; next}
    $c["shard"]==s{line=$c[w[1]]; for(j=2;j<=n;j++) line=line OFS $c[w[j]]; print line}' "$LISTS/pilot_genomes.tsv"
}

# Refuse to start when the scratch file quota is nearly full (incident 2026-10-29). Unparseable output only warns.
inode_guard() {
  local used
  used=$( { lfs quota -q -u "$USER" /scratch 2>/dev/null || true; } | awk 'NF>=6 && $1 ~ /^\// {gsub(/\*/,"",$6); print $6; exit} NF>=5 && $1 !~ /^\// {gsub(/\*/,"",$5); print $5; exit}' || true)
  if [[ "$used" =~ ^[0-9]+$ ]]; then
    echo "scratch files in use: $used (limit for this pilot: $INODE_LIMIT)"
    (( used < INODE_LIMIT )) || { echo "scratch file count $used >= $INODE_LIMIT: not starting" >&2; exit 4; }
  else
    echo "WARNING: could not read the scratch file quota (lfs quota); check kuq by hand" >&2
  fi
}

# Run a command inside the container. Binds: /project/k10018 (tools, references, code, JAR) read-only, the work area and the node-local
# dir read-write, the package tree's home (testenv) and JDK read-only.
ctr() {
  singularity exec --no-home --no-mount tmp --pwd "${CTR_PWD:-$LOCAL}" \
    --bind /project/k10018:/project/k10018:ro \
    --bind "$OUT:$OUT" --bind "$LOCAL:$LOCAL" \
    --bind "$PKG_ROOT/home/leechuck:/home/leechuck:ro" \
    --bind "$PKG_ROOT/usr/lib/jvm/java-21-openjdk-amd64:/usr/lib/jvm/java-21-openjdk-amd64:ro" \
    --env TMPDIR="$LOCAL/tmp" --env JAVA_TOOL_OPTIONS="-Djava.io.tmpdir=$LOCAL/tmp" --env PYTHONDONTWRITEBYTECODE=1 \
    --env PATH="/usr/lib/jvm/java-21-openjdk-amd64/bin:/home/leechuck/tmp/gspa-v3-testenv/bin:/usr/bin:/bin:$TOOLS/bin" \
    --env LC_ALL=C.UTF-8 --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
    "$IMAGE" "$@"
}

# Step receipt (inputs, versions, hashes, /usr/bin/time -v resources), see receipt.py.
receipt() { "$PY_HOST" "$HERE_CODE/receipt.py" "$@"; }
phase() { receipt phase --file "$PHASES" --name "$1" --event "$2"; }   # phase <name> start|end

# Job receipt: receipts/<step>/<jobid>.json (+ one resources TSV with every timed process). Called from the
# EXIT trap of every job, so a failed or time-limited job still records what it was charged.
RARGS=()
job_receipt() {
  local step=$1 code=$2; shift 2
  mkdir -p "$RECEIPTS/$step"
  receipt write --step "$step" --unit "job-${SLURM_JOB_ID:-local}" --kind job \
    --out "$RECEIPTS/$step/${SLURM_JOB_ID:-local}.json" --phases "$PHASES" --times-root "$LOCAL/times" \
    --job-start "$JOB_START" --nodes "${SLURM_JOB_NUM_NODES:-1}" --cores-per-node "$CORES_PER_NODE" \
    --versions-file "$VERSIONS" --resources-tsv "$RECEIPTS/$step/${SLURM_JOB_ID:-local}.resources.tsv.gz" \
    --param "shards=${SHARDS:-}" --param "exit_code=$code" --param "local_root=$LOCAL_ROOT" "${RARGS[@]}" "$@"
}
# install_finish <step>: EXIT trap that writes the receipt and removes the node-local dir; TERM (time limit,
# scancel) exits through it as well.
install_finish() {
  STEP=$1
  finish() { local code=$?; trap - EXIT; set +e; job_receipt "$STEP" "$code"; rm -rf "$LOCAL"; exit "$code"; }
  trap finish EXIT; trap 'exit 143' TERM USR1
}
# timed <phase> <tool> <item> <command...>: run under /usr/bin/time -v, record in ${TIMES_ROOT:-$LOCAL/times}
timed() {
  local p=$1 t=$2 i=$3 r=${TIMES_ROOT:-$LOCAL/times}; shift 3
  mkdir -p "$r/$p/$t"
  if declare -F "$1" > /dev/null; then   # a shell function (ctr): GNU time can only exec programs
    /usr/bin/time -v -o "$r/$p/$t/$i.time" bash -c '"$@"' _ "$@"
  else
    /usr/bin/time -v -o "$r/$p/$t/$i.time" "$@"
  fi
}
export -f timed ctr
