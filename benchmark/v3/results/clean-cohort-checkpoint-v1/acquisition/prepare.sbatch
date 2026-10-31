#!/bin/bash
#SBATCH --job-name=gspa-clean-cohort-prep
#SBATCH --partition=debug
#SBATCH --nodelist=node004
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:40:00
#SBATCH --output=/data/hohndor/gspa-v3/clean-cohort-v2/slurm/%j.log
set -euo pipefail
root=/data/hohndor/gspa-v3/clean-cohort-v2
export TMPDIR=/scratch/$SLURM_JOB_ID
mkdir -p "$TMPDIR" "$root/setup/$SLURM_JOB_ID"
finish() { code=$?; trap - EXIT; printf '%s\n' "$code" > "$root/setup/$SLURM_JOB_ID/exit-code.txt"; exit "$code"; }
trap finish EXIT
cd "$root/source"
sha256sum -c "$root/source.sha256" > "$root/setup/$SLURM_JOB_ID/source-check.log"
python3 -m unittest discover -s benchmark/v3/tests -p test_observe_clean_cohort.py -v > "$root/setup/$SLURM_JOB_ID/observation-tests.log" 2>&1
python3 -m unittest discover -s benchmark/neural/tests -p test_clean_contract.py -v > "$root/setup/$SLURM_JOB_ID/adapter-tests.log" 2>&1
test ! -e "$root/env"
/storage/miniforge3/envs/test/bin/python -m venv --system-site-packages "$root/env"
"$root/env/bin/python" -m pip install numpy==1.22.3 pandas==1.4.2 scikit-learn==1.2.0 scipy==1.7.3 matplotlib==3.7.0 tqdm==4.64.0 fair-esm==2.0.0 > "$root/setup/$SLURM_JOB_ID/dependencies.log" 2>&1
python3 benchmark/v3/observe_clean_cohort.py prepare \
 --cohort /data/hohndor/gspa-v3/cohort-pathway-search-v1/manifest.json \
 --expected-cohort b9491693f7b34b29abe99c7f077c44f2b70b3003db258816d487d3badf2e8fa9 \
 --app /data/hohndor/gspa-neural/models/CLEAN/app \
 --torch-home /data/hohndor/gspa-v3/clean-contract-v1/torch-cache \
 --python "$root/env/bin/python" --out "$root/collection"
sha256sum "$root/collection/plan.json" > "$root/plan.sha256"
