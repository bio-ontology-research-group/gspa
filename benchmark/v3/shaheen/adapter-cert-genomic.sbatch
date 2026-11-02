#!/bin/bash
#SBATCH --job-name=gspa-cert-genomic
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=04:00:00
#SBATCH --exclude=nid00022
#SBATCH --output=/scratch/hohndor/gspa-cert/logs/cert-genomic-%j.out
# Row-3 adapter certification of the geNomad and CheckV sidecar adapters on BW25113:
# 1. genuine direct sidecar run (native TSV = oracle), 2. generic route harness
# (annotate / integrated / direct-integrate, injected failures, attribution).
set -u
module load singularity

ROOT=/scratch/hohndor/gspa-cert
export TMPDIR=$ROOT/tmp; mkdir -p $TMPDIR $ROOT/results/genomad $ROOT/results/checkv
PY=/project/k10018/sw/miniforge3/bin/python3; PY_DIR=$(dirname $PY)   # explicit 3.12: ~/.pyenv shims (3.6) shadow the python module on compute nodes
export PATH=$(dirname $(which singularity)):$PY_DIR:/usr/bin:/bin
$PY --version > $ROOT/results/python_version.txt
JAVA=/project/k10018/sw/jdk/current/bin/java
SING_DIR=$(dirname $(which singularity))
singularity --version > $ROOT/results/singularity_version.txt
# CheckV's database tarball ships checkv_reps.faa; `checkv download_database` would build the DIAMOND index.
CVDB=$ROOT/db/checkv-db-v1.5
if [ ! -f $CVDB/genome_db/checkv_reps.dmnd ]; then
  singularity exec --bind $CVDB:/db $ROOT/images/checkv.sif diamond makedb --in /db/genome_db/checkv_reps.faa --db /db/genome_db/checkv_reps --threads 32 > $ROOT/logs/checkv_makedb.log 2>&1
  echo "checkv diamond makedb rc=$?"
fi
for P in genomad checkv; do
  OBS=$ROOT/results/$P/observed; mkdir -p $OBS
  if [ $P = genomad ]; then SIF=/project/k10018/sw/images/eq-annot/genomad.sif; DB=$ROOT/db/genomad_db; EXTRA="--genomad-sif $SIF"
  else SIF=$ROOT/images/checkv.sif; DB=$ROOT/db/checkv-db-v1.5; EXTRA="--checkv-sif $SIF --threads 32"; fi
  singularity exec $SIF $P --version > $OBS/tool_version.txt 2>&1 || singularity exec $SIF $P -h 2>&1 | head -3 > $OBS/tool_version.txt
  printf "tag\tgenome_fasta\tgff_path\toutput_dir\nquery\t$ROOT/bw25113/GCF_000750555.1.fna\t$ROOT/bw25113/GCF_000750555.1.gff\t$OBS\n" > $OBS/manifest.tsv
  echo "$PY $ROOT/run_genomic_predictors.py --predictor $P --manifest $OBS/manifest.tsv --min-score $([ $P = genomad ] && echo 0.7 || echo 0.5) --db-path $DB $EXTRA" > $OBS/command.txt
  t0=$(date +%s.%N)
  $PY $ROOT/run_genomic_predictors.py --predictor $P --manifest $OBS/manifest.tsv --min-score $([ $P = genomad ] && echo 0.7 || echo 0.5) --db-path $DB $EXTRA > $OBS/sidecar.log 2>&1
  rc=$?
  echo "{\"exit_code\": $rc, \"seconds\": $(echo "$(date +%s.%N) - $t0" | bc), \"date\": \"$(date -Is)\", \"node\": \"$(hostname)\"}" > $OBS/run.json
  echo "$P direct sidecar exit=$rc"; ls -la $OBS
  if [ ! -s $OBS/query.$P.genomic.tsv ]; then echo "$P: no native TSV, harness skipped"; continue; fi
  $PY $ROOT/make_genomic_specs.py $ROOT $P $SIF $DB $SING_DIR $OBS/query.$P.genomic.tsv $PY
  rm -rf $ROOT/results/$P/after
  $PY $ROOT/check_adapter_route.py --spec $ROOT/results/$P/spec.json --jar $ROOT/gspa-after.jar --java $JAVA --out $ROOT/results/$P/after > $ROOT/results/$P/after.log 2>&1
  echo "$P harness rc=$?"; tail -3 $ROOT/results/$P/after.log
done
echo CERT_GENOMIC_DONE
