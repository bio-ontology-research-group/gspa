#!/bin/bash
# Phase 10 benchmark — SLURM array job runner.
#
# Submitted as an array: SLURM_ARRAY_TASK_ID in [0..N-1] where N = |CONFIGS| * |TAGS|.
# Task id decomposes as (cfg_idx, tag_idx) = (id / |TAGS|, id % |TAGS|).
#
# One config (C5_iter_cluster_reps) runs the JAR inside Singularity to
# satisfy the "at least one config through Singularity" requirement.
# All other configs run native Java (no container overhead per job).
#SBATCH --job-name=gspa-phase10
#SBATCH --partition=debug
#SBATCH --exclude=node007
#SBATCH --output=/data/hohndor/gspa/proteomes/bench10/integrated_phase10/slurm-%A_%a.out
#SBATCH --error=/data/hohndor/gspa/proteomes/bench10/integrated_phase10/slurm-%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=4

set -euo pipefail

ROOT=/data/hohndor/gspa/proteomes/bench10
REF=/data/hohndor/gspa/reference
JAR=/data/hohndor/gspa/bin/gspa-phase10.jar
JAVA=/data/hohndor/envs/jdk-21.0.10+7/bin/java
JAVA_IMG=/data/hohndor/gspa/nf-test/singularity_cache/eclipse-temurin-21-jre.img

cd ${ROOT}

# ---- Config / tag arrays ----
CONFIGS=(C1_baseline C2_iterate C3_iter_cluster C4_iter_cluster_blastp C5_iter_cluster_reps C2_iterate_nopin)
TAGS=(vcholerae saureus spneumoniae ccrescentus rprowazekii tpallidum tthermophilus dradiodurans scoelicolor pfuriosus)

declare -A CONFIG_FLAGS=(
  [C1_baseline]=''
  [C2_iterate]='--dark-matter --iterate-gapseq'
  [C2_iterate_nopin]='--dark-matter --iterate-gapseq --gapseq-pin-promotions false'
  [C3_iter_cluster]='--dark-matter --iterate-gapseq --intragenome-cluster 0.9'
  [C4_iter_cluster_blastp]='--dark-matter --iterate-gapseq --intragenome-cluster 0.9 --gapseq-target proteome'
  [C5_iter_cluster_reps]='--dark-matter --iterate-gapseq --intragenome-cluster 0.9 --gapseq-target reps'
)

declare -A KINGDOM=(
  [vcholerae]=bacteria [saureus]=bacteria [spneumoniae]=bacteria [ccrescentus]=bacteria
  [rprowazekii]=bacteria [tpallidum]=bacteria [tthermophilus]=bacteria [dradiodurans]=bacteria
  [scoelicolor]=bacteria [pfuriosus]=archaea
)

# ---- Decompose SLURM_ARRAY_TASK_ID ----
N_TAGS=${#TAGS[@]}
TID=${SLURM_ARRAY_TASK_ID:-0}
CFG_IDX=$(( TID / N_TAGS ))
TAG_IDX=$(( TID % N_TAGS ))
CFG=${CONFIGS[$CFG_IDX]}
TAG=${TAGS[$TAG_IDX]}
K=${KINGDOM[$TAG]}
FLAGS=${CONFIG_FLAGS[$CFG]}

echo "=== task ${TID}: ${CFG} / ${TAG} on $(hostname) ==="
echo "    flags: ${FLAGS:-(none)}"

# Pick claims file (dpi_merged preferred, falls back to dp_claims)
CLAIMS=""
for cand in claims/${TAG}_dpi_merged.jsonl claims/${TAG}_merged.jsonl claims/${TAG}_dp_claims.jsonl; do
  [[ -s ${cand} ]] && { CLAIMS=${cand}; break; }
done
[[ -n ${CLAIMS} ]] || { echo "no claims for ${TAG}"; exit 2; }

OUT_DIR=integrated_phase10/${CFG}
mkdir -p ${OUT_DIR}
OUT_TSV=${OUT_DIR}/${TAG}_integrated.tsv
SUG_TSV=${OUT_DIR}/${TAG}_suggestions.tsv
LOG=${OUT_DIR}/${TAG}.log

[[ -s ${OUT_TSV} ]] && { echo "already done: ${OUT_TSV}"; exit 0; }

OP_ARG=""
GP_ARG=""
[[ -s ops/${TAG}_operons_up.tsv ]] && OP_ARG="--operons ops/${TAG}_operons_up.tsv"
[[ -s gaps/${TAG}_gaps.jsonl ]] && GP_ARG="--gaps gaps/${TAG}_gaps.jsonl"

# Build the java invocation; wrap in singularity for C5.
INVOKE_CMD=("${JAVA}" -jar "${JAR}")
if [[ ${CFG} == C5_iter_cluster_reps ]]; then
  # Stage image to /tmp with verification + retry to handle GlusterFS
  # read-consistency issues where a worker occasionally sees a corrupt
  # version of the cached file.
  STAGED_IMG=/tmp/gspa_${SLURM_JOB_ID:-$$}_jre.img
  trap "rm -f ${STAGED_IMG}" EXIT
  STAGED=0
  for attempt in 1 2 3 4 5; do
    # Pre-warm: full cat read forces FUSE to pull fresh bytes from the
    # authoritative brick.
    cat ${JAVA_IMG} > /dev/null 2>&1 || true
    rm -f ${STAGED_IMG}
    cp ${JAVA_IMG} ${STAGED_IMG}
    if singularity inspect ${STAGED_IMG} >/dev/null 2>&1; then
      echo "    [singularity mode via ${STAGED_IMG} (stage attempt ${attempt})]"
      STAGED=1
      break
    fi
    echo "    stage attempt ${attempt} produced corrupt image; retrying"
    sleep 2
  done
  if [[ ${STAGED} -eq 0 ]]; then
    echo "    [singularity staging failed after 5 attempts; falling back to native Java]"
    INVOKE_CMD=("${JAVA}" -jar "${JAR}")
  else
    INVOKE_CMD=(singularity exec --bind /data:/data ${STAGED_IMG} java -jar "${JAR}")
  fi
fi

time "${INVOKE_CMD[@]}" integrate \
  --claims ${CLAIMS} \
  --out ${OUT_TSV} \
  --suggestions-out ${SUG_TSV} \
  --go-owl ${REF}/go.owl \
  --lite \
  --essential-profile ${K} \
  --pathways ${REF}/kegg_pathways.tsv \
  --ec2go ${REF}/ec2go.txt \
  --enable-priors essentiality,coherence,gap_filling,genomic_context \
  ${OP_ARG} ${GP_ARG} \
  ${FLAGS} \
  > ${LOG} 2>&1

RC=$?
if [[ -s ${OUT_TSV} ]]; then
  echo "  wrote $(wc -l < ${OUT_TSV}) rows"
  exit 0
else
  echo "  no output (rc=${RC}); see ${LOG}"
  exit 1
fi
