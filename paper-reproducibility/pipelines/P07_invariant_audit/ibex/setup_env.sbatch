#!/bin/bash
#SBATCH --job-name=p07_env
#SBATCH --partition=batch
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/encrypted/e3008/gspa/meteor108_repro_2026-10-20/logs/env_%j.out
set -euo pipefail
# Python env for the P07 extraction (numpy/scipy/pandas; cobra/pulp/highspy are kept so that the archived
# METEOR MILP pickles and bounds files load, nothing of METEOR's own code is imported).
D=${W:-/encrypted/e3008/gspa/meteor108_repro_2026-10-20}
module load python/3.11.0
python3 -m venv $D/venv
source $D/venv/bin/activate
pip install -q --upgrade pip
pip install -q "numpy==1.23.4" "pandas==1.5.3" "scipy==1.13.1" "cobra==0.31.1" "optlang==1.9.1" "pulp==3.3.1" "highspy==1.15.1"
pip freeze > $D/venv_freeze.txt
python3 -u $D/scripts/inspect_meteor.py
