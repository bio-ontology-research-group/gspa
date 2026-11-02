#!/usr/bin/python3
"""Pass-through wrapper for run_genomic_predictors.py that keeps a copy of the sidecar TSV the
adapter is about to parse (the adapter deletes its scratch dir). Needed for tools that are not
deterministic between runs, so the oracle is read from the same run that produced the claims."""
import os, shutil, subprocess, sys
args = sys.argv[1:]
rows = open(args[args.index('--manifest') + 1]).read().splitlines()
row = dict(zip(rows[0].split('\t'), rows[1].split('\t')))
r = subprocess.run([sys.executable, os.environ['GSPA_CAPTURE_REAL']] + args)
src = os.path.join(row['output_dir'], row['tag'] + '.' + args[args.index('--predictor') + 1] + '.genomic.tsv')
if os.path.exists(src):
    shutil.copy(src, os.path.join(os.environ['GSPA_CAPTURE_DIR'], os.environ['GSPA_CAPTURE_NAME']))
sys.exit(r.returncode)
