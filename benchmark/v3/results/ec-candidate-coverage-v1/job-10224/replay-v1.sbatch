#!/bin/bash
#SBATCH --job-name=gspa-ec-reference-report
#SBATCH --partition=debug
#SBATCH --nodelist=node002
#SBATCH --cpus-per-task=3
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/data/hohndor/gspa-v3/candidate-coverage-v1/slurm/%j.log
set -euo pipefail
current=/data/hohndor/gspa-v3/candidate-coverage-v1
package=/data/hohndor/gspa-v3/observed-gap-context-modes-v1-relocation
work=/scratch/$SLURM_JOB_ID/ec-reference-report
jobout=$current/jobs/$SLURM_JOB_ID
mkdir -p "$work/fs" "$jobout"
export TMPDIR=/scratch/$SLURM_JOB_ID/tmp
mkdir -p "$TMPDIR"
finish() {
 code=$?; trap - EXIT; set +e
 if [ -d "$work/fs/home/leechuck/tmp/full-reference-panel" ]; then rsync -a "$work/fs/home/leechuck/tmp/full-reference-panel/" "$jobout/results/" || code=90; fi
 if [ -d "$work/fs/home/leechuck/Public/software/gspa/benchmark/v3/reference-panel-replay-inputs" ]; then rsync -a "$work/fs/home/leechuck/Public/software/gspa/benchmark/v3/reference-panel-replay-inputs/" "$jobout/inputs/" || code=91; fi
 printf '%s\n' "$code" > "$jobout/exit-code.txt"
 date --iso-8601=seconds > "$jobout/finished-at.txt"
 exit "$code"
}
trap finish EXIT
hostname > "$jobout/hostname.txt"
date --iso-8601=seconds > "$jobout/started-at.txt"
cp "$current/replay-v1.sbatch" "$current/prepare-replay.py" "$current/replay-report.py" "$jobout/"
sha256sum "$package/package.tar.zst" "$current/source-v1.tar.gz" "$current/prepare-replay.py" "$current/replay-report.py" > "$jobout/packages.sha256"
zstd -dc "$package/package.tar.zst" | tar -xf - -C "$work/fs"
tar -xf /data/hohndor/gspa-v3/context-resume-v1/supplement.tar -C "$work/fs"
tar -xf "$current/source-v1.tar.gz" -C "$work/fs"
cp "$current/replay-report.py" "$work/fs/home/leechuck/tmp/replay-report.py"
python3 "$current/prepare-replay.py" "$work/fs" > "$jobout/prepare.log" 2>&1
container=(singularity exec --no-home --no-mount tmp --pwd /home/leechuck/Public/software/gspa
 --bind /usr,/lib,/lib64,/bin,/etc --bind "$work/fs/home/leechuck:/home/leechuck"
 --env TMPDIR=/home/leechuck/tmp --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1
 "$package/container-root")
"${container[@]}" /home/leechuck/tmp/gspa-v3-testenv/bin/python -m unittest -v benchmark.v3.tests.test_observed_gap_context_analysis benchmark.v3.tests.test_ec_reference_panel > "$jobout/tests.log" 2>&1
"${container[@]}" /home/leechuck/tmp/gspa-v3-testenv/bin/python /home/leechuck/tmp/replay-report.py > "$jobout/replay.log" 2>&1
