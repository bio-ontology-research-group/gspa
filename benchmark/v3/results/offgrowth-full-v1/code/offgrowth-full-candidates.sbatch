#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=offgrowth-full-candidates
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --exclusive
#SBATCH --time=03:00:00
#SBATCH --no-requeue
#SBATCH --exclude=nid00022,nid00497,nid01033
#
# Full-panel off-growth run (tasks/v3-offgrowth-full-panel-2026-11-02.md), candidate step on one node:
#   1. leakage search: every hidden protein against the repaired genus-fenced Swiss-Prot reference
#      (DIAMOND --sensitive, e <= 1e-5, unlimited targets); exact / >= 95 % identical entries are removed;
#   2. verification search against the final reference (must find none);
#   3. the homology-baseline adapter search (training parameters) of all 145 proteomes against the final
#      reference, 48 genomes in parallel x 4 threads, and one compact EC candidate table per genome.
# Submit (login node):  sbatch --output=$OUT_ROOT/slurm/candidates-%j.log --export=ALL,REPO_ROOT=...,OUT_ROOT=... <this file>
set -euo pipefail
: "${REPO_ROOT:?REPO_ROOT must be set explicitly}"
: "${OUT_ROOT:?OUT_ROOT must be set explicitly}"
PY=${PYTHON_BIN:-/project/k10018/sw/envs/gspa-v3/bin/python3}
DIAMOND_BIN=${DIAMOND_BIN:-/project/k10018/sw/envs/diamond-2.1.11/bin}
export PATH="$DIAMOND_BIN:$PATH" PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1
export PYTHONPATH="$REPO_ROOT/benchmark/v3"
IN=$OUT_ROOT/inputs
W=$OUT_ROOT/candidate-search
mkdir -p "$W/leakage" "$W/final-reference" "$OUT_ROOT/candidates"
FP=(-I -c "import sys; sys.path.insert(0, '$REPO_ROOT/benchmark/v3'); from offgrowth import full_panel as F; sys.exit(F.main())")
{ date --iso-8601=seconds; hostname; diamond version; sha256sum "$DIAMOND_BIN/diamond"; } > "$W/environment.txt"
sha256sum "$IN"/* > "$W/inputs.sha256"
OUTFMT=(--outfmt 6 qseqid sseqid pident qcovhsp scovhsp evalue bitscore)

# 1. leakage search on the fenced reference
diamond makedb --in "$IN/reference.faa" -d "$W/leakage/fenced" --threads 64 > "$W/leakage/makedb.log" 2>&1
diamond blastp -q "$IN/hidden_proteins.faa" -d "$W/leakage/fenced.dmnd" -o "$W/leakage/hits_fenced.tsv" --sensitive \
  --evalue 1e-5 --max-target-seqs 0 --threads 64 "${OUTFMT[@]}" > "$W/leakage/blastp_fenced.log" 2>&1
"$PY" "${FP[@]}" leakage --hidden-faa "$IN/hidden_proteins.faa" --reference-faa "$IN/reference.faa" \
  --reference-labels "$IN/reference_labels.tsv" --hits "$W/leakage/hits_fenced.tsv" --out-dir "$W/final-reference" \
  --audit "$W/leakage/leakage_audit.json"

# 2. verification on the final reference
diamond makedb --in "$W/final-reference/reference.faa" -d "$W/leakage/final" --threads 64 > "$W/leakage/makedb_final.log" 2>&1
diamond blastp -q "$IN/hidden_proteins.faa" -d "$W/leakage/final.dmnd" -o "$W/leakage/hits_final.tsv" --sensitive \
  --evalue 1e-5 --max-target-seqs 0 --threads 64 "${OUTFMT[@]}" > "$W/leakage/blastp_final.log" 2>&1
"$PY" "${FP[@]}" leakage --verify --hidden-faa "$IN/hidden_proteins.faa" --reference-faa "$W/final-reference/reference.faa" \
  --hits "$W/leakage/hits_final.tsv" --audit "$W/leakage/leakage_verify.json"
rm -f "$W/leakage/"*.dmnd

# 3. candidate search, 145 genomes
/usr/bin/time -v -o "$W/search_time_v.txt" "$PY" "${FP[@]}" search --genomes "$IN/genomes_shaheen.tsv" \
  --reference-faa "$W/final-reference/reference.faa" --reference-labels "$W/final-reference/reference_labels.tsv" \
  --out-dir "$W/genomes" --cand-dir "$OUT_ROOT/candidates" --receipt "$W/search_receipt.json" --workers 48 \
  > "$W/search.log" 2>&1
( cd "$W" && tar -cf genomes.tar genomes && rm -rf genomes )
sha256sum "$OUT_ROOT/candidates/"*.candidates.tsv > "$OUT_ROOT/candidates/candidates.sha256"
date --iso-8601=seconds >> "$W/environment.txt"
