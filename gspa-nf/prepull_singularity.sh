#!/bin/bash
# Pre-pull singularity images for the GSPA Nextflow pipeline
set -euo pipefail
export SINGULARITY_CACHEDIR=/data/hohndor/gspa/nf-test/singularity_cache
mkdir -p ${SINGULARITY_CACHEDIR}
cd ${SINGULARITY_CACHEDIR}

# Images we need for the M. genitalium D+P+barrnap+minced run
IMAGES=(
  'docker://quay.io/biocontainers/pyrodigal:3.7.1--py312h247cb63_1'
  'docker://quay.io/biocontainers/diamond:2.1.9--h43eeafb_0'
  'docker://quay.io/biocontainers/hmmer:3.4--hb6cb901_4'
  'docker://quay.io/biocontainers/barrnap:0.8--0'
  'docker://quay.io/biocontainers/minced:0.3.0--0'
  'docker://python:3.12-slim'
)

for img in "${IMAGES[@]}"; do
  # Nextflow uses a specific name format for the cache file
  name=$(echo "${img}" | sed 's|docker://||; s|/|-|g; s|:|-|g').img
  if [[ -s "${name}" ]]; then
    echo "${name}: already cached ($(ls -lh ${name} | awk '{print $5}'))"
    continue
  fi
  echo "Pulling ${img} -> ${name}..."
  rm -f "${name}"
  singularity pull --name "${name}" "${img}" 2>&1 | tail -3
  echo "  done: $(ls -lh ${name} 2>/dev/null | awk '{print $5}')"
done

echo "All images pulled."
ls -lh ${SINGULARITY_CACHEDIR}/
