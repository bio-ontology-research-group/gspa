#!/bin/bash
#SBATCH --job-name=gspa-cert-genomic-h
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=03:00:00
#SBATCH --exclude=nid00022
#SBATCH --output=/scratch/hohndor/gspa-cert/logs/cert-genomic-harness-%j.out
# Harness-only rerun (after a harness fix) on the genuine sidecar outputs already captured by adapter-cert-genomic.sbatch.
set -u
module load singularity
ROOT=/scratch/hohndor/gspa-cert
export TMPDIR=$ROOT/tmp; mkdir -p $TMPDIR
PY=/project/k10018/sw/miniforge3/bin/python3
export PATH=$(dirname $(which singularity)):$(dirname $PY):/usr/bin:/bin
JAVA=/project/k10018/sw/jdk/current/bin/java
SING_DIR=$(dirname $(which singularity))
for P in ${PREDICTORS:-genomad checkv}; do
  OBS=$ROOT/results/$P/observed
  if [ $P = genomad ]; then SIF=/project/k10018/sw/images/eq-annot/genomad.sif; DB=$ROOT/db/genomad_db; else SIF=$ROOT/images/checkv.sif; DB=$ROOT/db/checkv-db-v1.5; fi
  [ -s $OBS/query.$P.genomic.tsv ] || { echo "$P: no native TSV"; continue; }
  $PY $ROOT/make_genomic_specs.py $ROOT $P $SIF $DB $SING_DIR $OBS/query.$P.genomic.tsv $PY
  [ -d $ROOT/results/$P/after ] && mv $ROOT/results/$P/after $ROOT/results/$P/after.$(date +%H%M%S)
  $PY $ROOT/check_adapter_route.py --spec $ROOT/results/$P/spec.json --jar $ROOT/gspa-after.jar --java $JAVA --out $ROOT/results/$P/after > $ROOT/results/$P/after.log 2>&1
  echo "$P harness rc=$?"; grep "^{" $ROOT/results/$P/after.log
done
echo CERT_GENOMIC_HARNESS_DONE
