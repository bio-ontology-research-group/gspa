#!/bin/bash
#SBATCH --job-name=gspa-extracellular-products
#SBATCH --partition=debug
#SBATCH --nodelist=node004
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:40:00
#SBATCH --no-requeue
#SBATCH --output=/data/hohndor/gspa-v3/extracellular-products-v2/slurm/%j.log
set -euo pipefail
root=/data/hohndor/gspa-v3/extracellular-products-v2
parent=/data/hohndor/gspa-v3/empty-input-v1
work=/scratch/$SLURM_JOB_ID/extracellular
out=$root/jobs/$SLURM_JOB_ID
mkdir -p "$work/fs" "$out"
export TMPDIR=/scratch/$SLURM_JOB_ID/tmp
mkdir -p "$TMPDIR"
finish() {
 code=$?
 trap - EXIT
 set +e
 if [ -d "$work/fs/home/leechuck/tmp/extracellular-results" ]; then
  rsync -a "$work/fs/home/leechuck/tmp/extracellular-results/" "$out/results/" || code=90
 fi
 printf '%s\n' "$code" > "$out/exit-code.txt"
 exit "$code"
}
trap finish EXIT
cd "$root"
sha256sum -c packages-v2.sha256 > "$out/package-check.log"
cp "$0" "$out/executed.sbatch"
hostname > "$out/hostname.txt"
zstd -dc "$parent/before.tar.zst" | tar -xf - -C "$work/fs"
zstd -dc "$parent/python-runtime.tar.zst" | tar -xf - -C "$work/fs"
tar -xzf "$root/source-v2.tar.gz" -C "$work/fs"
cd "$work/fs"
sha256sum -c "$root/source-v2.sha256" > "$out/source-check.log"
repo=$work/fs/home/leechuck/Public/software/gspa
mkdir -p "$work/fs/home/leechuck/tmp/extracellular-results"
container=(singularity exec --no-home --no-mount tmp --pwd /home/leechuck/Public/software/gspa
 --bind /usr,/lib,/lib64,/bin,/etc
 --bind "$work/fs/home/leechuck:/home/leechuck"
 --env PYTHONPATH=/home/leechuck/Public/software/gspa/scripts/tests:/home/leechuck/Public/software/gspa/benchmark/v3
 --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1
 --env XDG_CACHE_HOME=/home/leechuck/tmp/cache --env TMPDIR=/home/leechuck/tmp
 /data/hohndor/gspa-v3/observed-gap-context-modes-v1-relocation/container-root
 /home/leechuck/tmp/gspa-v3-testenv/bin/python)
cp "$repo/scripts/metabolic_experiments.py" "$work/after.py"
cp "$repo/scripts/masked_assignment_experiments.py" "$work/after-masked.py"
cp "$work/fs/home/leechuck/tmp/extracellular-before/scripts/metabolic_experiments.py" "$repo/scripts/metabolic_experiments.py"
cp "$work/fs/home/leechuck/tmp/extracellular-before/scripts/masked_assignment_experiments.py" "$repo/scripts/masked_assignment_experiments.py"
set +e
"${container[@]}" -m unittest -v scripts.tests.test_extracellular_products > "$out/before-tests.log" 2>&1
before=$?
set -e
printf '%s\n' "$before" > "$out/before-exit-code.txt"
test "$before" -ne 0
cp "$work/after.py" "$repo/scripts/metabolic_experiments.py"
cp "$work/after-masked.py" "$repo/scripts/masked_assignment_experiments.py"
/usr/bin/time -v -o "$out/resources.txt" "${container[@]}" -m unittest -v scripts.tests.test_extracellular_products scripts.tests.test_intracellular_products scripts.tests.test_metabolic_experiments scripts.tests.test_masked_assignment_experiments scripts.tests.test_absence_experiments scripts.tests.test_metabolic_planner_feedback scripts.tests.test_knockout_identity scripts.tests.test_community_knockout_identity scripts.tests.test_community_metabolism > "$out/tests.log" 2>&1
"${container[@]}" -m unittest -v benchmark.v3.tests.test_run_metabolic_benchmark benchmark.v3.tests.test_check_metabolic_cobra benchmark.v3.tests.test_masked_assignment_oracle benchmark.v3.tests.test_metabolic_protocols benchmark.v3.tests.test_planner_benchmark > "$out/benchmark-tests.log" 2>&1

"${container[@]}" benchmark/v3/check_extracellular_products.py --model /home/leechuck/tmp/extracellular-input/Sco-GEM.xml --medium /home/leechuck/tmp/extracellular-input/minimal_medium.csv --out /home/leechuck/tmp/extracellular-results/real-model > "$out/diagnostic.log" 2>&1
