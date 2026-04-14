#!/bin/bash
# Run GSPA Nextflow on M. genitalium via SLURM + singularity
set -euo pipefail
cd /data/hohndor/gspa/nf-test

# Java 21 needed for Nextflow
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate java21

NEXTFLOW=/home/hohndor/bin/nextflow
INPUT=/data/hohndor/gspa/proteomes/mgenitalium_genomic.fna
DIAMOND_DB=/data/hohndor/gspa/proteomes/reference_loo9.dmnd
PFAM_DB=/storage/software/databases/hmmer/Pfam-A.hmm
OUTDIR=/data/hohndor/gspa/nf-test/results

# Clean previous run
rm -rf work .nextflow* "${OUTDIR}/mgenitalium" 2>/dev/null || true

${NEXTFLOW} run gspa-nf/main.nf \
  -c gspa-nf/nextflow.config \
  -c slurm_singularity.config \
  --input "${INPUT}" \
  --diamond_db "${DIAMOND_DB}" \
  --pfam_db "${PFAM_DB}" \
  --outdir "${OUTDIR}" \
  --kingdom bacteria \
  -profile singularity \
  -with-report results/report.html \
  -with-trace results/trace.txt
