#!/bin/bash
#SBATCH --job-name=gspa-context-resume
#SBATCH --partition=debug
#SBATCH --nodelist=node002
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=/data/hohndor/gspa-v3/context-resume-v1/slurm/%j.log
set -euo pipefail
package=/data/hohndor/gspa-v3/observed-gap-context-modes-v1-relocation
current=/data/hohndor/gspa-v3/context-resume-v1
work=/scratch/$SLURM_JOB_ID/gspa-context
jobout=$current/jobs/$SLURM_JOB_ID
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
  if [ -d "$src" ]; then
   base="$package/jobs/10163/results/$name"
   mkdir -p "$jobout/results/$name"
   if [ -d "$base" ]; then
    rsync -a --compare-dest="$base" "$src/" "$jobout/results/$name/" || code=90
   else
    rsync -a "$src/" "$jobout/results/$name/" || code=90
   fi
  fi
 done
 if [ -d "$work/fs/home/leechuck/tmp/gspa-context-resume-execution" ]; then
  rsync -a "$work/fs/home/leechuck/tmp/gspa-context-resume-execution/" "$jobout/logs/" || code=91
 fi
 printf '%s\n' "$code" > "$jobout/exit-code.txt"
 date --iso-8601=seconds > "$jobout/finished-at.txt"
 exit "$code"
}
trap finish EXIT
hostname > "$jobout/hostname.txt"
date --iso-8601=seconds > "$jobout/started-at.txt"
cp "$current/run.sbatch" "$current/run-resume.py" "$current/verify-supplement.py" "$current/supplement.json" "$jobout/"
printf '%s\n' "$package/jobs/10163/results" > "$jobout/unchanged-export-base.txt"
sha256sum "$package/package.tar.zst" "$current/supplement.tar" "$current/supplement.json" "$current/run-resume.py" "$current/verify-supplement.py" > "$jobout/packages.sha256"
mkdir -p "$work/fs"
zstd -dc "$package/package.tar.zst" | tar -xf - -C "$work/fs"
mkdir -p "$work/fs/home/leechuck/tmp"
rsync -a "$package/jobs/10163/results/observed-gap-context-modes-v1/" "$work/fs/home/leechuck/.cache/gspa-v3-data/observed-gap-context-modes-v1/"
tar -xf "$current/supplement.tar" -C "$work/fs"
cp "$current/supplement.json" "$work/fs/home/leechuck/tmp/context-input-supplement.json"
cp "$current/verify-supplement.py" "$work/fs/home/leechuck/tmp/verify-context-supplement.py"
cp "$current/run-resume.py" "$work/fs/home/leechuck/tmp/run-context-remote.py"
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
