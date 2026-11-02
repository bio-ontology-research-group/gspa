#!/bin/bash
#SBATCH --job-name=gspa-cert-gapseq-h
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=04:00:00
#SBATCH --exclude=nid00022,nid00497,nid01033
#SBATCH --output=/scratch/hohndor/gspa-cert/logs/cert-gapseq-harness-%j.out
# Harness-only rerun for the gapseq nucleotide route on the genuine Reactions.tbl already captured.
# (second submission 2026-11-02: harness honours the spec's cli_threads/timeout_seconds; job 15745883
# had run gapseq with -K 1 and hit the harness's 3600 s limit before any check.)
set -u
ROOT=/scratch/hohndor/gspa-cert
export TMPDIR=$ROOT/tmp; mkdir -p $TMPDIR
PY=/project/k10018/sw/miniforge3/bin/python3
GAPSEQ=/project/k10018/sw/envs/gapseq/bin/gapseq
export PATH=/project/k10018/sw/envs/gapseq/bin:/usr/bin:/bin
JAVA=/project/k10018/sw/jdk/current/bin/java
NATIVE=$(ls $ROOT/results/gapseq/observed/out/*-Reactions.tbl | head -1)
$PY $ROOT/make_genomic_specs.py $ROOT gapseq $GAPSEQ /project/k10018/sw/envs/gapseq/bin x $NATIVE $PY
[ -d $ROOT/results/gapseq/after ] && mv $ROOT/results/gapseq/after $ROOT/results/gapseq/after.$(date +%H%M%S)
$PY $ROOT/check_adapter_route.py --spec $ROOT/results/gapseq/spec.json --jar $ROOT/gspa-after.jar --java $JAVA --out $ROOT/results/gapseq/after > $ROOT/results/gapseq/after.log 2>&1
echo "gapseq harness rc=$?"; grep "^{" $ROOT/results/gapseq/after.log
echo CERT_GAPSEQ_HARNESS_DONE
