#!/usr/bin/env bash
# Shared configuration for rebuilding the 17,055 Empty Quarter metabolic models
# and the 145 benchmark metabolic models (22 Price + 18 cohort + 105 peers)
# with gapsmith 0.1.0 on the Shaheen III supercomputer.

set -euo pipefail

export SHAHEEN_SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export EQ_MAGS_DIR="$(cd "$SHAHEEN_SCRIPTS_DIR/.." && pwd)"
export REPO_ROOT="$(cd "$EQ_MAGS_DIR/../.." && pwd)"

# Read-only Empty Quarter inputs on Shaheen (/project/k10018/empty-quarter)
export EQ_ROOT="${EQ_ROOT:-/project/k10018/empty-quarter}"
export EQ_SITES_META="${EQ_SITES_META:-$EQ_ROOT/meta/sites.tsv}"
export EQ_SAMPLES_META="${EQ_SAMPLES_META:-$EQ_ROOT/meta/samples.tsv}"
export EQ_SITES_TAR_DIR="${EQ_SITES_TAR_DIR:-$EQ_ROOT/data/sites-tar}"
# RAK-PR culture genome catalogue (722 genomes: catalogue.tsv, inputs.tsv,
# results/annotation/bakta/<gid>/<gid>.{faa,gff3,fna}, and raw assemblies)
export EQ_RAKPR_DIR="${EQ_RAKPR_DIR:-$EQ_ROOT/data/rakpr}"

# Reference per-bin table from the October IBEX build (17,055 rows)
export EQ_PER_BIN_TSV="${EQ_PER_BIN_TSV:-$EQ_MAGS_DIR/results/eq_mags_per_bin.tsv}"

# Benchmark genome inputs on Shaheen
export PRICE_DIR="${PRICE_DIR:-/project/k10018/gspa-v3/price-full-matrix-v2}"
export PRICE_TARGETS_JSON="${PRICE_TARGETS_JSON:-$PRICE_DIR/targets.json}"
export PRICE_LAYOUT_REPORT="${PRICE_LAYOUT_REPORT:-$PRICE_DIR/layout/report.json}"
export PRICE_GENOMES_DIR="${PRICE_GENOMES_DIR:-/project/k10018/gspa-v3/price-raw-sources/genomes}"
export WP4_SEARCH_DIR="${WP4_SEARCH_DIR:-/scratch/project/k10018/gspa-v3/out/wp4-search-v1}"
export PEERS_TSV="${PEERS_TSV:-$REPO_ROOT/benchmark/v3/results/peer-selection-v1/peers.tsv}"

# Tool binaries on Shaheen (/project/k10018)
# October IBEX tools (sha256 in /project/k10018/sw/gapsmith-eq-ibex/bins.sha256): gapsmith 7b7ff538…, DIAMOND 2.1.24 (bakta env),
# each behind a wrapper that sets its own library path; reference = October data_merged, one tar unpacked per node.
export EQ_IBEX_TOOLS="${EQ_IBEX_TOOLS:-/project/k10018/sw/gapsmith-eq-ibex}"
export GAPSMITH_BIN="${GAPSMITH_BIN:-$EQ_IBEX_TOOLS/wrap/gapsmith}"
export GAPSMITH_DATA_TAR="${GAPSMITH_DATA_TAR:-$EQ_IBEX_TOOLS/data_merged.tar}"
export GAPSMITH_BACTERIA_REF="${GAPSMITH_BACTERIA_REF:-/project/k10018/gspa-v3/gapsmith-reference-bacteria-v1/reference}"
export GAPSMITH_ARCHAEA_REF="${GAPSMITH_ARCHAEA_REF:-/project/k10018/gspa-v3/archaea-reference-v1/reference}"
export PYRODIGAL_BIN="${PYRODIGAL_BIN:-/project/k10018/sw/envs/gspa-tools/bin/pyrodigal}"
export EXPECTED_PYRODIGAL_VERSION="${EXPECTED_PYRODIGAL_VERSION:-3.7.1}"
export EXPECTED_GAPSMITH_VERSION="${EXPECTED_GAPSMITH_VERSION:-0.1.0}"
export DIAMOND_BIN="${DIAMOND_BIN:-$EQ_IBEX_TOOLS/wrap/diamond}"
export PYTHON_BIN="${PYTHON_BIN:-/project/k10018/sw/envs/gspa-tools/bin/python3}"

# Writable scratch workspace on Shaheen (/scratch/project/k10018)
export WORK_ROOT="${WORK_ROOT:-/scratch/project/k10018/gspa-v3/eq-models-v1}"
# Always derived from WORK_ROOT: an inherited value (sbatch --export=ALL from a shell that sourced this file
# for another mode) once sent the benchmark job to the Empty Quarter batches (2026-10-30).
export GAPSMITH_DATA_MERGED="$WORK_ROOT/gapsmith_data_merged"
export BATCHES_DIR="$WORK_ROOT/batches"
export LOG_DIR="$WORK_ROOT/logs"

# Slurm and node packing parameters (Shaheen III workq: 192 physical cores, ~376 GB RAM)
export SLURM_ACCOUNT="${SLURM_ACCOUNT:-k10018}"
export SLURM_PARTITION="${SLURM_PARTITION:-workq}"
export SLURM_EXCLUDE="${SLURM_EXCLUDE:-nid00022,nid00497,nid01033}"
export SLURM_TIME="${SLURM_TIME:-24:00:00}"
export CORES_PER_NODE="${CORES_PER_NODE:-192}"
export THREADS_PER_GENOME="${THREADS_PER_GENOME:-4}"
export JOBS_PER_NODE="${JOBS_PER_NODE:-48}"
export BATCH_SIZE="${BATCH_SIZE:-500}"
