#!/bin/bash
#SBATCH --job-name=gspa_benchmark
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=12:00:00
#SBATCH --output=benchmark_%j.log

# GSPA Benchmark Step 2: Run all available predictors on test proteins
#
# Usage: sbatch 02_run_predictors.sh /data/hohndor/gspa/benchmark

set -euo pipefail

WORKDIR="${1:-/data/hohndor/gspa/benchmark}"
DATADIR="${WORKDIR}/benchmark_data"
RESULTSDIR="${WORKDIR}/predictor_results"
THREADS=${SLURM_CPUS_PER_TASK:-16}

mkdir -p "${RESULTSDIR}"

echo "=== GSPA Benchmark: Running predictors ==="
echo "Work dir: ${WORKDIR}"
echo "Threads: ${THREADS}"
echo "Date: $(date)"

# Paths to databases
DIAMOND_DB="${WORKDIR}/benchmark_data/reference_db.dmnd"
PFAM_DB="/storage/software/databases/hmmer/Pfam-A.hmm"
FOLDSEEK_DB="/data/databases/foldseek/pdb"  # Use PDB structures DB
# eggNOG: system has eggnog.db but needs eggnog_proteins.dmnd too
# If the proteins DMND is at the custom path, symlink it
EGGNOG_DB="/storage/software/databases/eggnog"
EGGNOG_CUSTOM="/data/hohndor/gspa/databases/eggnog/eggnog_proteins.dmnd"
if [ -f "${EGGNOG_CUSTOM}" ] && [ ! -f "${EGGNOG_DB}/eggnog_proteins.dmnd" ]; then
    # Create a merged dir with symlinks
    EGGNOG_MERGED="${WORKDIR}/eggnog_merged"
    mkdir -p "${EGGNOG_MERGED}"
    ln -sf "${EGGNOG_DB}/eggnog.db" "${EGGNOG_MERGED}/"
    ln -sf "${EGGNOG_DB}/eggnog.taxa.db" "${EGGNOG_MERGED}/"
    ln -sf "${EGGNOG_DB}/eggnog.taxa.db.traverse.pkl" "${EGGNOG_MERGED}/"
    ln -sf "${EGGNOG_CUSTOM}" "${EGGNOG_MERGED}/"
    EGGNOG_DB="${EGGNOG_MERGED}"
fi

TEST_FASTA="${DATADIR}/test_proteins.fasta"
NPROTEINS=$(grep -c "^>" "${TEST_FASTA}")
echo "Test proteins: ${NPROTEINS}"

# --- 1. Build reference DIAMOND DB (Swiss-Prot minus test proteins) ---
echo ""
echo "=== Building reference DIAMOND database ==="
if [ ! -f "${DIAMOND_DB}" ]; then
    diamond makedb --in "${DATADIR}/reference_db.fasta" --db "${DIAMOND_DB%.dmnd}" --threads "${THREADS}"
fi
echo "DIAMOND DB ready: $(ls -lh ${DIAMOND_DB})"

# --- 2. DIAMOND blastp ---
echo ""
echo "=== Running DIAMOND ==="
diamond blastp \
    --db "${DIAMOND_DB}" \
    --query "${TEST_FASTA}" \
    --out "${RESULTSDIR}/diamond_results.tsv" \
    --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore stitle \
    --evalue 1e-5 --max-target-seqs 50 \
    --query-cover 50 --subject-cover 50 --id 30 \
    --threads "${THREADS}"
echo "DIAMOND hits: $(wc -l < ${RESULTSDIR}/diamond_results.tsv)"

# --- 3. HMMER/Pfam ---
echo ""
echo "=== Running HMMER/Pfam ==="
hmmsearch \
    --domtblout "${RESULTSDIR}/pfam_results.domtbl" \
    --noali -E 1e-5 \
    --cpu "${THREADS}" \
    "${PFAM_DB}" "${TEST_FASTA}" > /dev/null
echo "Pfam hits: $(grep -vc '^#' ${RESULTSDIR}/pfam_results.domtbl)"

# --- 4. MMseqs2 ---
echo ""
echo "=== Running MMseqs2 ==="
MMSEQS_TMP="${RESULTSDIR}/mmseqs_tmp"
mkdir -p "${MMSEQS_TMP}"
# Build MMseqs2 DB from reference if not exists
if [ ! -f "${WORKDIR}/benchmark_data/reference_db_mmseqs" ]; then
    mmseqs createdb "${DATADIR}/reference_db.fasta" "${WORKDIR}/benchmark_data/reference_db_mmseqs"
fi
mmseqs createdb "${TEST_FASTA}" "${MMSEQS_TMP}/queryDB"
mmseqs search "${MMSEQS_TMP}/queryDB" "${WORKDIR}/benchmark_data/reference_db_mmseqs" \
    "${MMSEQS_TMP}/resultDB" "${MMSEQS_TMP}/tmp" \
    -e 1e-5 --max-seqs 50 --threads "${THREADS}"
mmseqs convertalis "${MMSEQS_TMP}/queryDB" "${WORKDIR}/benchmark_data/reference_db_mmseqs" \
    "${MMSEQS_TMP}/resultDB" "${RESULTSDIR}/mmseqs2_results.tsv" \
    --format-output "query,target,pident,alnlen,qlen,tlen,evalue,bits"
echo "MMseqs2 hits: $(wc -l < ${RESULTSDIR}/mmseqs2_results.tsv)"

# --- 5. FoldSeek (ProstT5 if available, otherwise skip) ---
echo ""
echo "=== Running FoldSeek ==="
if command -v foldseek &>/dev/null && [ -d "/data/databases/foldseek" ]; then
    # Use PDB database with sequence search
    foldseek easy-search \
        "${TEST_FASTA}" \
        "/data/databases/foldseek/pdb" \
        "${RESULTSDIR}/foldseek_results.tsv" \
        "${RESULTSDIR}/foldseek_tmp" \
        --format-output "query,target,pident,evalue,bits,theader" \
        --threads "${THREADS}" \
        -e 1e-3 2>&1 | tail -5 || echo "FoldSeek failed (may need ProstT5 model)"
    echo "FoldSeek hits: $(wc -l < ${RESULTSDIR}/foldseek_results.tsv 2>/dev/null || echo 0)"
else
    echo "FoldSeek not available, skipping"
fi

# --- 6. eggNOG-mapper ---
echo ""
echo "=== Running eggNOG-mapper ==="
if [ -f "${EGGNOG_DB}/eggnog.db" ]; then
    # Need emapper.py — check if it's in the metagenomics env or a separate one
    if command -v emapper.py &>/dev/null; then
        emapper.py \
            -i "${TEST_FASTA}" \
            -o "${RESULTSDIR}/eggnog_results" \
            --data_dir "${EGGNOG_DB}" \
            -m diamond \
            --cpu "${THREADS}" \
            --itype proteins \
            --override 2>&1 | tail -5
        echo "eggNOG hits: $(grep -vc '^#' ${RESULTSDIR}/eggnog_results.emapper.annotations 2>/dev/null || echo 0)"
    else
        echo "emapper.py not found, skipping"
    fi
else
    echo "eggNOG DB not found at ${EGGNOG_DB}, skipping"
fi

echo ""
echo "=== All predictors complete ==="
echo "Results in: ${RESULTSDIR}"
ls -lh "${RESULTSDIR}/"
echo "Date: $(date)"
