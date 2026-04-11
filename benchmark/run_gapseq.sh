#!/bin/bash
# Run gapseq find on all 9 genomes, using the writable gapseq install
# at /data/hohndor/envs/gapseq-rw. gapseq's uniprot.sh needs to write into
# its own share/dat tree, so the conda-installed read-only copy fails.
set -u
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate gapseq-env

GAPSEQ=/data/hohndor/envs/gapseq-rw/gapseq
ROOT=/data/hohndor/gspa/proteomes
OUT=${ROOT}/gapseq
mkdir -p ${OUT}

GENOMES=(ecoli hpylori mgenitalium mjannaschii ecolo157 bsubtilis mtb synechocystis paeruginosa)

MAX_JOBS=3
running=0

run_one() {
  local tag=$1
  local dir=${OUT}/${tag}
  mkdir -p ${dir}
  if [[ -s ${dir}/${tag}-all-Pathways.tbl && -s ${dir}/${tag}-all-Reactions.tbl ]]; then
    echo "  ${tag}: already done, skipping"
    return 0
  fi
  # Clean any stale partial outputs
  rm -f ${dir}/*.fasta ${dir}/*.log 2>/dev/null
  cp ${ROOT}/${tag}_genomic.fna ${dir}/${tag}.fna
  pushd ${dir} >/dev/null
  echo "  [${tag}] starting gapseq find..."
  # -p all: all MetaCyc pathways. gapseq parallelizes internally; we
  # limit concurrency via MAX_JOBS above.
  time ${GAPSEQ} find -p all -v 0 ${tag}.fna > gapseq_find.log 2>&1
  rc=$?
  echo "  [${tag}] done rc=${rc}"
  popd >/dev/null
  return ${rc}
}

for tag in "${GENOMES[@]}"; do
  (run_one ${tag}) &
  running=$((running + 1))
  if (( running >= MAX_JOBS )); then
    wait -n
    running=$((running - 1))
  fi
done
wait
echo "=== ALL GAPSEQ DONE ==="
for tag in "${GENOMES[@]}"; do
  pt=${OUT}/${tag}/${tag}-all-Pathways.tbl
  rt=${OUT}/${tag}/${tag}-all-Reactions.tbl
  if [[ -s ${pt} && -s ${rt} ]]; then
    echo "  ${tag}: $(grep -vc '^#' ${pt}) pathways, $(grep -vc '^#' ${rt}) reactions"
  else
    echo "  ${tag}: FAILED or missing"
  fi
done
