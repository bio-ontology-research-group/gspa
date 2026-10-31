#!/bin/bash
#SBATCH --job-name=gspa-knockout-identity
#SBATCH --partition=debug
#SBATCH --cpus-per-task=2
#SBATCH --nodelist=node002
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/data/hohndor/gspa-v3/knockout-identity-v1/slurm/%j.log
set -euo pipefail
package=/data/hohndor/gspa-v3/empty-input-v1
current=/data/hohndor/gspa-v3/knockout-identity-v1
work=/scratch/$SLURM_JOB_ID/gspa-knockout-identity
jobout=$current/jobs/$SLURM_JOB_ID
mkdir -p "$work/fs" "$jobout"
export TMPDIR=/scratch/$SLURM_JOB_ID/tmp
mkdir -p "$TMPDIR"
finish() { code=$?; trap - EXIT; if [ -d "$work/fs/home/leechuck/tmp/germicidin-diagnostic" ]; then rsync -a "$work/fs/home/leechuck/tmp/germicidin-diagnostic/" "$jobout/diagnostic/"; fi; printf '%s\n' "$code" > "$jobout/exit-code.txt"; exit "$code"; }
trap finish EXIT
hostname > "$jobout/hostname.txt"
cp "$current/validate-v4.sbatch" "$jobout/"
sha256sum "$package/before.tar.zst" "$package/python-runtime.tar.zst" "$current/source-v4.tar.gz" > "$jobout/packages.sha256"
zstd -dc "$package/before.tar.zst" | tar -xf - -C "$work/fs"
zstd -dc "$package/python-runtime.tar.zst" | tar -xf - -C "$work/fs"
tar -xf "$current/source-v4.tar.gz" -C "$work/fs"
mkdir -p "$work/fs/home/leechuck/tmp"
container=(singularity exec --no-home --no-mount tmp --pwd /home/leechuck/Public/software/gspa
 --bind /usr,/lib,/lib64,/bin,/etc --bind "$work/fs/home/leechuck:/home/leechuck"
 --env PYTHONPATH=/home/leechuck/Public/software/gspa/scripts/tests --env TMPDIR=/home/leechuck/tmp --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1
 --env PATH=/home/leechuck/tmp/gspa-v3-testenv/bin:/usr/bin:/bin
 /data/hohndor/gspa-v3/observed-gap-context-modes-v1-relocation/container-root)
repo="$work/fs/home/leechuck/Public/software/gspa"
cp "$repo/scripts/metabolic_experiments.py" "$work/after.py"
cp "$repo/scripts/community_metabolism.py" "$work/after-community.py"
cp "$work/fs/home/leechuck/tmp/before-knockout/metabolic_experiments.py" "$repo/scripts/metabolic_experiments.py"
cp "$work/fs/home/leechuck/tmp/before-knockout/community_metabolism.py" "$repo/scripts/community_metabolism.py"
set +e
"${container[@]}" python -m unittest -v scripts.tests.test_knockout_identity scripts.tests.test_community_knockout_identity > "$jobout/before-tests.log" 2>&1
before_code=$?
set -e
printf '%s\n' "$before_code" > "$jobout/before-exit-code.txt"
cp "$work/after.py" "$repo/scripts/metabolic_experiments.py"
cp "$work/after-community.py" "$repo/scripts/community_metabolism.py"
"${container[@]}" python -m unittest -v scripts.tests.test_metabolic_experiments scripts.tests.test_masked_assignment_experiments scripts.tests.test_absence_experiments scripts.tests.test_metabolic_planner_feedback scripts.tests.test_knockout_identity scripts.tests.test_community_knockout_identity scripts.tests.test_community_metabolism > "$jobout/tests.log" 2>&1
"${container[@]}" python -m unittest -v benchmark.v3.tests.test_run_metabolic_benchmark benchmark.v3.tests.test_check_metabolic_cobra benchmark.v3.tests.test_metabolic_protocols benchmark.v3.tests.test_planner_benchmark benchmark.v3.tests.test_run_community_benchmark benchmark.v3.tests.test_check_community_cobra benchmark.v3.tests.test_nutrient_community > "$jobout/benchmark-tests.log" 2>&1


"${container[@]}" python benchmark/v3/run_germicidin_diagnostic.py --model /home/leechuck/tmp/germicidin-input/Sco-GEM.xml --medium /home/leechuck/tmp/germicidin-input/minimal_medium.csv --before /home/leechuck/tmp/before-knockout/metabolic_experiments.py --out /home/leechuck/tmp/germicidin-diagnostic > "$jobout/diagnostic.log" 2>&1
