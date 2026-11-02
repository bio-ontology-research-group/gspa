#!/bin/bash
#SBATCH --job-name=gspa-cert-gapseq
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=10:00:00
#SBATCH --exclude=nid00022
#SBATCH --output=/scratch/hohndor/gspa-cert/logs/cert-gapseq-%j.out
# Row-3 certification of the gapseq nucleotide route (GapseqPredictor tool=gapseq target=genome)
# on BW25113: 1. genuine `gapseq find -p all -t Bacteria -M nucl` with the adapter's exact
# arguments (native Reactions.tbl kept), 2. the generic route harness.
set -u
ROOT=/scratch/hohndor/gspa-cert
export TMPDIR=$ROOT/tmp; mkdir -p $TMPDIR $ROOT/results/gapseq/observed
PY=/project/k10018/sw/miniforge3/bin/python3
GAPSEQ=/project/k10018/sw/envs/gapseq/bin/gapseq
export PATH=/project/k10018/sw/envs/gapseq/bin:/usr/bin:/bin
JAVA=/project/k10018/sw/jdk/current/bin/java
OBS=$ROOT/results/gapseq/observed
$GAPSEQ -v > $OBS/tool_version.txt 2>&1
echo "$GAPSEQ find -p all -t Bacteria -M nucl -K 32 -f $OBS/out $ROOT/bw25113/GCF_000750555.1.fna" > $OBS/command.txt
mkdir -p $OBS/out; t0=$(date +%s.%N)
cd $OBS && $GAPSEQ find -p all -t Bacteria -M nucl -K 32 -f $OBS/out $ROOT/bw25113/GCF_000750555.1.fna > $OBS/stdout.txt 2> $OBS/stderr.txt
rc=$?
echo "{\"exit_code\": $rc, \"seconds\": $(echo "$(date +%s.%N) - $t0" | bc), \"date\": \"$(date -Is)\", \"node\": \"$(hostname)\"}" > $OBS/run.json
echo "gapseq direct exit=$rc"; ls -la $OBS/out
NATIVE=$(ls $OBS/out/*-Reactions.tbl 2>/dev/null | head -1)
[ -n "$NATIVE" ] || { echo "no Reactions.tbl"; exit 1; }
$PY $ROOT/make_genomic_specs.py $ROOT gapseq $GAPSEQ /project/k10018/sw/envs/gapseq/bin x $NATIVE $PY
[ -d $ROOT/results/gapseq/after ] && mv $ROOT/results/gapseq/after $ROOT/results/gapseq/after.$(date +%H%M%S)
$PY $ROOT/check_adapter_route.py --spec $ROOT/results/gapseq/spec.json --jar $ROOT/gspa-after.jar --java $JAVA --out $ROOT/results/gapseq/after > $ROOT/results/gapseq/after.log 2>&1
echo "gapseq harness rc=$?"; grep "^{" $ROOT/results/gapseq/after.log
echo CERT_GAPSEQ_DONE
