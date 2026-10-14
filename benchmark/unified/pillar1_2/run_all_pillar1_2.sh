#!/bin/bash
set -euo pipefail

PYTHON="/storage/miniforge3/envs/metagenomics/bin/python3"
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=========================================================="
echo "Phase 3 Unified Benchmark: Pillar 1 & Pillar 2"
echo "Host: $(hostname) | Date: $(date)"
echo "Directory: ${DIR}"
echo "=========================================================="

echo ""
echo ">>> STEP 1: Running 9-Genome x 9-Ablation Integration Grid..."
${PYTHON} "${DIR}/01_run_ablation_grid.py" -j 8

echo ""
echo ">>> STEP 2: Running CAFA F-max Evaluations (Micro & Protein-Centric, Bootstrap=50)..."
${PYTHON} "${DIR}/02_run_cafa_eval.py" -j 8 --n-bootstrap 50

echo ""
echo ">>> STEP 3: Running GAEF Genome-Scale Quality Evaluations..."
${PYTHON} "${DIR}/03_run_gaef_eval.py" -j 8

echo ""
echo ">>> STEP 4: Running Cross-Clade Taxon & Anatomical Perturbation Experiment..."
${PYTHON} "${DIR}/04_run_perturbation.py" -j 8

echo ""
echo ">>> STEP 5: Compiling Results Tables & Markdown Report..."
${PYTHON} "${DIR}/05_compile_results.py"

echo ""
echo "=========================================================="
echo "Pillar 1 & Pillar 2 Benchmark Run Finished!"
echo "Date: $(date)"
echo "=========================================================="
