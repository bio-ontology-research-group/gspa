#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=shared
#SBATCH --job-name=gspa-arm
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=15G
#SBATCH --time=06:00:00
#SBATCH --no-requeue
#SBATCH --output=/scratch/project/k10018/gspa-v3/out/arm-logs/%A_%a.log
# Shaheen III port of the unimatrix arm-array recipes (role-genome-arm-v1/arm-v2.sbatch, job 10666;
# wp4-pilot-v1/wp4-pilot.sbatch, job 10686): one array task per arm of an arms.tsv, the recorded command JSON
# executed unchanged inside the same container (package + supplements, JAR and models from the data mirror).
# Required environment: BUNDLE (e.g. role-genome-arm-v1), ARMS (arms.tsv path relative to the bundle root),
# ARM_COLS (number of leading columns that name the arm: 2 for ctx,kind; 3 for regime,ctx,kind).
# The data mirror /project/k10018/gspa-v3 is bound read-only at the unimatrix path /data/hohndor/gspa-v3 and the
# bundle's writable copy on /scratch is bound at the bundle's original path, so every recorded path resolves.
set -euo pipefail
: "${BUNDLE:?}" "${ARMS:?}" "${ARM_COLS:?}"
module load singularity
# Container root: a plain Debian 13 sandbox (the unimatrix recipe bound the Debian host's /usr,/lib,/bin,/etc into an
# empty root; Shaheen compute nodes run SLES without /usr/lib/jvm, so the JDK bind point cannot be created there).
# Java and Python still come from the package; rsync from the gspa-tools conda env.
image=/project/k10018/sw/images/debian-trixie
data_root=/project/k10018/gspa-v3
logical_data=/data/hohndor/gspa-v3
out_root=/scratch/project/k10018/gspa-v3/out/$BUNDLE
root=$logical_data/$BUNDLE
work=/scratch/project/k10018/gspa-v3/work/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}/arm
mkdir -p "$work/fs" "$out_root/runs"
export TMPDIR=/scratch/project/k10018/gspa-v3/work/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}/tmp; mkdir -p "$TMPDIR"
# Writable bundle root: inputs copied once from the mirror (never runs/, slurm/, jobs/).
rsync -a --exclude runs/ --exclude slurm/ --exclude jobs/ --exclude 'runs-*/' --exclude 'mmseqs/tmp/' ${EXTRA_EXCLUDES:-} "$data_root/$BUNDLE/" "$out_root/"
row=$(awk -v i=$((SLURM_ARRAY_TASK_ID+1)) -F'\t' 'NR>1{n++} n==i{print; exit}' "$out_root/$ARMS")
[ -n "$row" ] || { echo "no arm row for task $SLURM_ARRAY_TASK_ID in $ARMS" >&2; exit 2; }
name=$(cut -f1-"$ARM_COLS" <<<"$row" | tr '\t' '/'); manifest=$(cut -f$((ARM_COLS+1)) <<<"$row"); command=$(cut -f$((ARM_COLS+2)) <<<"$row")
arm=$out_root/runs/$name; mkdir -p "$arm"
echo "task $SLURM_ARRAY_TASK_ID arm=$name manifest=$manifest command=$command host=$(hostname)" | tee "$arm/task.txt"
finish() { code=$?; trap - EXIT; set +e; printf '%s\n' "$code" > "$arm/exit-code.txt"; date --iso-8601=seconds > "$arm/finished-at.txt"
 [ "$code" -eq 0 ] && [ "${KEEP_WORK:-0}" != 1 ] && rm -rf "/scratch/project/k10018/gspa-v3/work/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"; exit "$code"; }
trap finish EXIT
date --iso-8601=seconds > "$arm/started-at.txt"; cp "$0" "$arm/executed.sbatch"
zstd -dc "$data_root/observed-gap-context-modes-v1-relocation/package.tar.zst" | tar -xf - -C "$work/fs"
tar -xf "$data_root/context-resume-v1/supplement.tar" -C "$work/fs"
tar -xf "$data_root/source-context-full-v1/parent-supplement.tar" -C "$work/fs"
mkdir -p "$work/fs/home/leechuck/tmp"
container=(singularity exec --no-home --no-mount tmp --pwd /home/leechuck
 --bind /project/k10018/sw/envs/gspa-tools:/project/k10018/sw/envs/gspa-tools:ro
 --bind "$data_root:$logical_data:ro" --bind "$out_root:$root"
 --bind "$work/fs/home/leechuck:/home/leechuck"
 --bind "$work/fs/usr/lib/jvm/java-21-openjdk-amd64:/usr/lib/jvm/java-21-openjdk-amd64"
 --env TMPDIR=/home/leechuck/tmp --env JAVA_TOOL_OPTIONS=-Djava.io.tmpdir=/home/leechuck/tmp
 --env PATH=/usr/lib/jvm/java-21-openjdk-amd64/bin:/usr/bin:/bin:/project/k10018/sw/envs/gspa-tools/bin --env LC_ALL=C.UTF-8 --env HOME=/home/leechuck --env XDG_CONFIG_HOME=/home/leechuck/.config --env XDG_CACHE_HOME=/home/leechuck/.cache --env XDG_DATA_HOME=/home/leechuck/.local/share
 --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1
 "$image")
printf '%s\n' "${container[@]}" > "$arm/container-command.txt"
"${container[@]}" java -version > "$arm/java-version.txt" 2>&1
# The command JSON lives at a /data path; read it through the mirror.
mapfile -t cmd < <(python3 -c 'import json,sys; [print(t) for t in json.load(open(sys.argv[1]))]' "${command/#$logical_data/$data_root}")
printf '%s\n' "${cmd[@]}" > "$arm/command.txt"
/usr/bin/time -v -o "$arm/resources.txt" "${container[@]}" "${cmd[@]}" > "$arm/integrate.log" 2>&1
echo done
