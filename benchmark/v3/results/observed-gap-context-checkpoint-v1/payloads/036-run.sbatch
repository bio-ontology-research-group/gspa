#!/bin/bash
#SBATCH --job-name=gspa-gap-context
#SBATCH --partition=debug
#SBATCH --nodelist=node007
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=/data/hohndor/gspa-v3/observed-gap-context-modes-v1-relocation/slurm/%j.log
set -euo pipefail
package=/data/hohndor/gspa-v3/observed-gap-context-modes-v1-relocation
work=/scratch/$SLURM_JOB_ID/gspa-context
jobout=$package/jobs/$SLURM_JOB_ID
mkdir -p "$work" "$jobout"
export TMPDIR=/scratch/$SLURM_JOB_ID/tmp
mkdir -p "$TMPDIR"
finish() {
 code=$?
 trap - EXIT
 set +e
 mkdir -p "$jobout/results"
 for name in observed-gap-context-modes-v1 observed-gap-context-modes-v1-parity-homology observed-gap-context-modes-v1-parity-family_excluded observed-gap-context-model-replay-v1 observed-gap-context-model-analysis-v1; do
  src="$work/fs/home/leechuck/.cache/gspa-v3-data/$name"
  if [ -d "$src" ]; then rsync -a "$src" "$jobout/results/" || code=90; fi
 done
 if [ -d "$work/fs/home/leechuck/tmp/gspa-context-remote-execution" ]; then
  rsync -a "$work/fs/home/leechuck/tmp/gspa-context-remote-execution/" "$jobout/logs/" || code=91
 fi
 printf '%s\n' "$code" > "$jobout/exit-code.txt"
 date --iso-8601=seconds > "$jobout/finished-at.txt"
 exit "$code"
}
trap finish EXIT
hostname > "$jobout/hostname.txt"
date --iso-8601=seconds > "$jobout/started-at.txt"
cp "$package/run.sbatch" "$package/run-remote.py" "$jobout/"
mkdir -p "$work/fs"
zstd -dc "$package/package.tar.zst" | tar -xf - -C "$work/fs"
mkdir -p "$work/fs/home/leechuck/tmp"
cp "$package/run-remote.py" "$work/fs/home/leechuck/tmp/run-context-remote.py"
container=(singularity exec --no-home --no-mount tmp --pwd /home/leechuck/Public/software/gspa
 --bind /usr,/lib,/lib64,/bin,/etc
 --bind "$work/fs/home/leechuck:/home/leechuck"
 --bind "$work/fs/usr/lib/jvm/java-21-openjdk-amd64:/usr/lib/jvm/java-21-openjdk-amd64"
 --env TMPDIR=/home/leechuck/tmp
 --env JAVA_TOOL_OPTIONS=-Djava.io.tmpdir=/home/leechuck/tmp
 --env PATH=/usr/lib/jvm/java-21-openjdk-amd64/bin:/usr/bin:/bin
 --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1
 "$package/container-root")
"${container[@]}" /usr/lib/jvm/java-21-openjdk-amd64/bin/java -version > "$jobout/java-version.txt" 2>&1
"${container[@]}" /home/leechuck/tmp/gspa-v3-testenv/bin/python -c 'import sys,numpy;print(sys.version);print(numpy.__version__)' > "$jobout/python-version.txt" 2>&1
"${container[@]}" /home/leechuck/tmp/gspa-v3-testenv/bin/python /home/leechuck/tmp/run-context-remote.py
