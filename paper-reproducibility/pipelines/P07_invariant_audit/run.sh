#!/bin/bash
# P07 -- 108-genome stoichiometric invariant audit (Fig. 6a, Suppl. S12).
# Hosts: ibex (extraction of the threshold drafts, METEOR solutions and Grow-Support repair from the archived
#        third-party results; SLURM only, never the login node) and local (LP audit + GSPA refiner; ~16 cores).
# All computation is GSPA code: scripts/meteor_threshold_repair.py (comparator) and
# scripts/stoichiometric_locus_refiner.py (GSPA Stage 2B); the archived METEOR tree is read as data only.
#
#   bash run.sh stage     # local: copy ibex/*, the comparator module, the panel and the data files to $W on ibex
#   bash run.sh extract   # local: sbatch the 108-task array (ibex/extract_array.sbatch)
#   bash run.sh merge     # local: sbatch ibex/merge.sbatch -> meteor108_drafts_compact.pkl.gz + manifest_sha256.tsv
#   bash run.sh fetch     # local: copy the compact archive + sha256 manifest into the repo
#   bash run.sh smoke     # local: smoke_test_reimpl.py -> results/smoke_test_reimpl.{tsv,summary.json}
#   bash run.sh audit     # local: run Experiment 2A -> results/pillar2_invariant_audit_{108,summary}.tsv
#   bash run.sh all-local # fetch + audit (when the ibex side is already done)
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd); root=$(cd "$here/../../.." && pwd)
W=/encrypted/e3008/gspa/meteor108_repro_2026-10-20          # our write area on ibex (e3008 project space)
METEOR_ROOT=${METEOR_ROOT:-/ibex/scratch/projects/c2014/kexin/funcarve}   # archived third-party results, read only
REFD=$root/benchmark/manifests/reference_data
DATA=$REFD/psb_revision/code_snapshot/data
mode=${1:-help}
case "$mode" in
  stage)
    ssh ibex "mkdir -p $W/scripts $W/data $W/logs $W/per_genome"
    scp -q "$here"/ibex/* "$root/scripts/meteor_threshold_repair.py" "$here/panel108_gram.tsv" ibex:$W/scripts/
    scp -q "$DATA"/fba_matrices_v6RE.pkl "$DATA"/tight_bounds_v6_neg.pkl "$DATA"/tight_bounds_v6_pos.pkl \
      "$DATA"/medium.pkl "$DATA"/seedr2ec.pkl "$DATA"/seedec2r.pkl "$DATA"/all_ancestors.txt "$DATA"/enzymeobsolete.pkl \
      "$REFD/universal_metadata.json" ibex:$W/data/ ;;
  extract) ssh ibex "METEOR_ROOT=$METEOR_ROOT sbatch --export=ALL,W=$W,METEOR_ROOT=$METEOR_ROOT $W/scripts/extract_array.sbatch" ;;
  merge)   ssh ibex "sbatch --export=ALL,W=$W $W/scripts/merge.sbatch" ;;
  fetch)
    scp -q ibex:$W/meteor108_drafts_compact.pkl.gz "$REFD/meteor108_drafts_compact.pkl.gz"
    scp -q ibex:$W/manifest_sha256.tsv "$here/ibex_sha256.tsv"   # raw; converted into manifest_rows.tsv
    sha256sum "$REFD/meteor108_drafts_compact.pkl.gz" ;;
  smoke)   python3 "$here/smoke_test_reimpl.py" "${@:2}" ;;
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
