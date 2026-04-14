#!/bin/bash
# Stage all singularity images cleanly via /tmp (avoids GlusterFS cache corruption)
set -euo pipefail

CACHE=/data/hohndor/gspa/nf-test/singularity_cache
mkdir -p ${CACHE}

declare -A IMAGES=(
  [quay.io-biocontainers-pyrodigal-3.7.1--py312h247cb63_1.img]='docker://quay.io/biocontainers/pyrodigal:3.7.1--py312h247cb63_1'
  [quay.io-biocontainers-diamond-2.1.9--h43eeafb_0.img]='docker://quay.io/biocontainers/diamond:2.1.9--h43eeafb_0'
  [quay.io-biocontainers-hmmer-3.4--hb6cb901_4.img]='docker://quay.io/biocontainers/hmmer:3.4--hb6cb901_4'
  [quay.io-biocontainers-barrnap-0.8--0.img]='docker://quay.io/biocontainers/barrnap:0.8--0'
  [quay.io-biocontainers-minced-0.3.0--0.img]='docker://quay.io/biocontainers/minced:0.3.0--0'
  [python-3.12.img]='docker://python:3.12'
)

for name in "${!IMAGES[@]}"; do
  uri="${IMAGES[$name]}"
  target="${CACHE}/${name}"
  tmp_file="/tmp/$$_${name}"

  # Check current state
  size=$(stat -c %s "${target}" 2>/dev/null || echo 0)
  # Test if it is readable by singularity
  if [[ ${size} -gt 1000000 ]] && singularity inspect "${target}" >/dev/null 2>&1; then
    echo "OK: ${name} (${size} bytes, inspectable)"
    continue
  fi

  echo "Re-staging: ${name}"
  rm -f "${target}" "${tmp_file}"
  SINGULARITY_CACHEDIR=/tmp/sc_$$ singularity pull --force --name "${tmp_file}" "${uri}" 2>&1 | tail -2
  cp "${tmp_file}" "${target}"
  rm -f "${tmp_file}"
  new_size=$(stat -c %s "${target}" 2>/dev/null || echo 0)
  echo "  -> ${new_size} bytes"
done

# Clean /tmp scratch caches
rm -rf /tmp/sc_$$

echo ""
echo "=== Final state ==="
for name in "${!IMAGES[@]}"; do
  target="${CACHE}/${name}"
  singularity inspect "${target}" >/dev/null 2>&1 && status="OK" || status="BROKEN"
  size=$(stat -c %s "${target}" 2>/dev/null || echo 0)
  printf "  %-7s %12d  %s\n" "${status}" "${size}" "${name}"
done
