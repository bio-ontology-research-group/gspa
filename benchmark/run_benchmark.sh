#!/bin/bash
# GSPA Weight Optimization Benchmark
# Run on unimatrix01: bash run_benchmark.sh
#
# This script:
# 1. Prepares a benchmark dataset from Swiss-Prot + GOA
# 2. Submits predictor jobs via SLURM
# 3. Computes optimal combination weights

set -euo pipefail

WORKDIR="/data/hohndor/gspa/benchmark"
SCRIPTDIR="$(cd "$(dirname "$0")" && pwd)"

echo "=== GSPA Weight Optimization Benchmark ==="
echo "Work dir: ${WORKDIR}"
echo "Scripts: ${SCRIPTDIR}"
echo ""

# Step 0: Check prerequisites
echo "--- Checking prerequisites ---"
for tool in diamond hmmsearch mmseqs python3; do
    which "$tool" > /dev/null 2>&1 || { echo "ERROR: $tool not found"; exit 1; }
done

# Check databases
[ -f "/storage/software/databases/hmmer/Pfam-A.hmm" ] || { echo "ERROR: Pfam DB not found"; exit 1; }
[ -f "/storage/software/databases/eggnog/eggnog.db" ] || echo "WARNING: eggNOG DB not found, will skip"
echo "Prerequisites OK"

# Step 1: Download data if needed
echo ""
echo "--- Step 1: Download Swiss-Prot + GOA ---"
mkdir -p "${WORKDIR}"
cd "${WORKDIR}"

if [ ! -f "uniprot_sprot.fasta.gz" ]; then
    echo "Downloading Swiss-Prot..."
    wget -q "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
fi

if [ ! -f "goa_uniprot_all.gaf.gz" ]; then
    echo "Downloading GOA..."
    wget -q "ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz"
fi

echo "Data files:"
ls -lh uniprot_sprot.fasta.gz goa_uniprot_all.gaf.gz

# Step 2: Prepare benchmark dataset
echo ""
echo "--- Step 2: Prepare benchmark dataset ---"
if [ ! -f "benchmark_data/test_proteins.fasta" ]; then
    python3 "${SCRIPTDIR}/01_prepare_benchmark.py" \
        --sprot uniprot_sprot.fasta.gz \
        --goa goa_uniprot_all.gaf.gz \
        --outdir benchmark_data \
        --n-test 5000
else
    echo "Benchmark data already prepared"
    echo "Test proteins: $(grep -c '^>' benchmark_data/test_proteins.fasta)"
fi

# Step 3: Submit predictor jobs via SLURM
echo ""
echo "--- Step 3: Run predictors ---"
JOBID=$(sbatch --parsable "${SCRIPTDIR}/02_run_predictors.sh" "${WORKDIR}")
echo "Submitted SLURM job: ${JOBID}"
echo "Monitor with: squeue -j ${JOBID}"
echo "Log: benchmark_${JOBID}.log"

# Wait for job to complete
echo "Waiting for job to complete..."
while squeue -j "${JOBID}" -h 2>/dev/null | grep -q "${JOBID}"; do
    sleep 30
    echo -n "."
done
echo ""
echo "Job ${JOBID} complete"

# Step 4: Compute weights
echo ""
echo "--- Step 4: Compute optimal weights ---"
python3 "${SCRIPTDIR}/03_compute_weights.py" \
    --results-dir "${WORKDIR}/predictor_results" \
    --data-dir "${WORKDIR}/benchmark_data" \
    --goa "${WORKDIR}/goa_uniprot_all.gaf.gz" \
    --output "${WORKDIR}/optimal_weights.json" \
    --weight-step 0.1

echo ""
echo "=== Benchmark complete ==="
echo "Results: ${WORKDIR}/optimal_weights.json"
cat "${WORKDIR}/optimal_weights.json"
