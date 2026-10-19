#!/bin/bash
# P07 -- 108-genome stoichiometric invariant audit (Fig. 6a, Suppl. S12).
# Hosts: ibex (extraction of METEOR v8 drafts / solutions / Grow-Support repair; SLURM only, never the
#        login node) and local (LP audit + GSPA refiner; ~16 cores, ~30-60 min).
#
#   bash run.sh stage     # local: copy ibex/* to $D/scripts on ibex
#   bash run.sh env       # local: sbatch ibex/setup_env.sbatch (venv with the archived package versions + schema probe)
#   bash run.sh extract   # local: sbatch the 108-task array (ibex/extract_array.sbatch)
#   bash run.sh merge     # local: sbatch ibex/merge.sbatch -> meteor108_drafts_compact.pkl.gz + manifest_sha256.tsv
#   bash run.sh fetch     # local: copy the compact archive + sha256 manifest into the repo
#   bash run.sh audit     # local: run Experiment 2A -> results/pillar2_invariant_audit_{108,summary}.tsv
#   bash run.sh all-local # fetch + audit (when the ibex side is already done)
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd); root=$(cd "$here/../../.." && pwd)
D=/ibex/scratch/projects/c2014/hohndor/meteor108_repro_2026-10-19   # our write area on ibex
REFD=$root/benchmark/manifests/reference_data
mode=${1:-help}
case "$mode" in
  stage)   ssh ibex "mkdir -p $D/scripts $D/logs $D/per_genome"; scp -q "$here"/ibex/* ibex:$D/scripts/ ;;
  env)     ssh ibex "sbatch $D/scripts/setup_env.sbatch" ;;
  extract) ssh ibex "sbatch $D/scripts/extract_array.sbatch" ;;
  merge)   ssh ibex "sbatch $D/scripts/merge.sbatch" ;;
  fetch)
    scp -q ibex:$D/meteor108_drafts_compact.pkl.gz "$REFD/meteor108_drafts_compact.pkl.gz"
    scp -q ibex:$D/manifest_sha256.tsv "$here/ibex_sha256.tsv"   # raw; converted into manifest_rows.tsv
    sha256sum "$REFD/meteor108_drafts_compact.pkl.gz" ;;
  audit)
    cd "$root"
    python3 - <<'PY'
import sys
sys.path.insert(0, "benchmark/unified/pillar2_fba_lro_price149")
import run_pillar2_all as m
ref = m.StoichiometricLocusRefiner(m.REF_DIR / "psb_revision" / "code_snapshot" / "data")
m.run_108_genome_invariant_audit(ref)
PY
    ;;
  all-local) "$0" fetch; "$0" audit ;;
  *) sed -n '2,15p' "$0" ;;
esac
