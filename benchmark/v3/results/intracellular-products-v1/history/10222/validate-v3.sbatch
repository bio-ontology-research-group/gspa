#!/bin/bash
#SBATCH --job-name=gspa-intracellular-products
#SBATCH --partition=debug
#SBATCH --cpus-per-task=3
#SBATCH --nodelist=node002
#SBATCH --mem=16G
#SBATCH --time=00:40:00
#SBATCH --output=/data/hohndor/gspa-v3/intracellular-products-v1/slurm/%j.log
set -euo pipefail
package=/data/hohndor/gspa-v3/empty-input-v1
current=/data/hohndor/gspa-v3/intracellular-products-v1
work=/scratch/$SLURM_JOB_ID/gspa-intracellular-products
jobout=$current/jobs/$SLURM_JOB_ID
mkdir -p "$work/fs" "$jobout"
export TMPDIR=/scratch/$SLURM_JOB_ID/tmp
mkdir -p "$TMPDIR"
finish() {
 code=$?
 trap - EXIT
 if [ -d "$work/fs/home/leechuck/tmp/product-diagnostic" ]; then rsync -a "$work/fs/home/leechuck/tmp/product-diagnostic/" "$jobout/diagnostic/"; fi
 if [ -d "$work/fs/home/leechuck/tmp/product-parity" ]; then rsync -a "$work/fs/home/leechuck/tmp/product-parity/" "$jobout/parity/"; fi
 if [ -f "$work/fs/home/leechuck/Public/software/gspa/gspa-cli/build/libs/gspa-3.0.0.jar" ]; then cp "$work/fs/home/leechuck/Public/software/gspa/gspa-cli/build/libs/gspa-3.0.0.jar" "$jobout/gspa-3.0.0.jar"; fi
 printf '%s\n' "$code" > "$jobout/exit-code.txt"
 for module in gspa-core gspa-predictors gspa-cli; do
  if [ -d "$work/fs/home/leechuck/Public/software/gspa/$module/build/test-results" ]; then
   mkdir -p "$jobout/final-tests/$module"
   rsync -a "$work/fs/home/leechuck/Public/software/gspa/$module/build/test-results/" "$jobout/final-tests/$module/" || true
  fi
 done
 exit "$code"
}
trap finish EXIT
hostname > "$jobout/hostname.txt"
cp "$current/validate-v3.sbatch" "$jobout/"
sha256sum "$package/before.tar.zst" "$package/python-runtime.tar.zst" "$current/source-v3.tar.gz" > "$jobout/packages.sha256"
zstd -dc "$package/before.tar.zst" | tar -xf - -C "$work/fs"
zstd -dc "$package/python-runtime.tar.zst" | tar -xf - -C "$work/fs"
tar -xf "$current/source-v3.tar.gz" -C "$work/fs"
repo=$work/fs/home/leechuck/Public/software/gspa
mkdir -p "$work/fs/home/leechuck/tmp"
container=(singularity exec --no-home --no-mount tmp --pwd /home/leechuck/Public/software/gspa
 --bind "$HOME/.config/borg-llm/api-key:/run/gspa-llm-key:ro"
 --bind /usr,/lib,/lib64,/bin,/etc
 --bind "$work/fs/home/leechuck:/home/leechuck"
 --bind "$work/fs/usr/lib/jvm/java-21-openjdk-amd64:/usr/lib/jvm/java-21-openjdk-amd64"
 --env GSPA_LLM_URL=http://10.72.186.139:8000/v1 --env PYTHONPATH=/home/leechuck/Public/software/gspa/scripts/tests --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1
 --env TMPDIR=/home/leechuck/tmp --env JAVA_TOOL_OPTIONS=-Djava.io.tmpdir=/home/leechuck/tmp
 --env JAVA_HOME=/usr/lib/jvm/java-21-openjdk-amd64
 --env PATH=/home/leechuck/tmp/gspa-v3-testenv/bin:/usr/lib/jvm/java-21-openjdk-amd64/bin:/usr/bin:/bin
 --env GSPA_TEST_PYTHON=/home/leechuck/tmp/gspa-v3-testenv/bin/python
 --env GRADLE_USER_HOME=/home/leechuck/.gradle
 /data/hohndor/gspa-v3/observed-gap-context-modes-v1-relocation/container-root)
cp "$repo/scripts/metabolic_experiments.py" "$work/after.py"
cp "$repo/scripts/masked_assignment_experiments.py" "$work/after-masked.py"
cp "$work/fs/home/leechuck/tmp/before-products/scripts/metabolic_experiments.py" "$repo/scripts/metabolic_experiments.py"
cp "$work/fs/home/leechuck/tmp/before-products/scripts/masked_assignment_experiments.py" "$repo/scripts/masked_assignment_experiments.py"
set +e
"${container[@]}" python -m unittest -v scripts.tests.test_intracellular_products > "$jobout/before-tests.log" 2>&1
before_code=$?
set -e
printf '%s\n' "$before_code" > "$jobout/before-exit-code.txt"
cp "$work/after.py" "$repo/scripts/metabolic_experiments.py"
cp "$work/after-masked.py" "$repo/scripts/masked_assignment_experiments.py"
"${container[@]}" python -m unittest -v scripts.tests.test_intracellular_products scripts.tests.test_metabolic_experiments scripts.tests.test_masked_assignment_experiments scripts.tests.test_absence_experiments scripts.tests.test_metabolic_planner_feedback scripts.tests.test_knockout_identity scripts.tests.test_community_knockout_identity scripts.tests.test_community_metabolism > "$jobout/tests.log" 2>&1
"${container[@]}" python -m unittest -v benchmark.v3.tests.test_run_metabolic_benchmark benchmark.v3.tests.test_check_metabolic_cobra benchmark.v3.tests.test_masked_assignment_oracle benchmark.v3.tests.test_metabolic_protocols benchmark.v3.tests.test_planner_benchmark benchmark.v3.tests.test_run_community_benchmark benchmark.v3.tests.test_check_community_cobra benchmark.v3.tests.test_nutrient_community > "$jobout/benchmark-tests.log" 2>&1
mkdir -p "$repo/gspa-cli/build/libs"
cp "$current/jobs/10221/gspa-3.0.0.jar" "$repo/gspa-cli/build/libs/gspa-3.0.0.jar"
sha256sum "$current/jobs/10221/gspa-3.0.0.jar" "$repo/gspa-cli/build/libs/gspa-3.0.0.jar" > "$jobout/validated-jar.sha256"
"${container[@]}" python benchmark/v3/check_planner_feedback.py --before /home/leechuck/tmp/before-products/scripts/metabolic_experiments.py --out /home/leechuck/tmp/product-parity > "$jobout/parity.log" 2>&1
"${container[@]}" python benchmark/v3/check_intracellular_products.py --model /home/leechuck/tmp/product-input/Sco-GEM.xml --medium /home/leechuck/tmp/product-input/minimal_medium.csv --jar gspa-cli/build/libs/gspa-3.0.0.jar --out /home/leechuck/tmp/product-diagnostic --live-qwen --key-file /run/gspa-llm-key > "$jobout/diagnostic.log" 2>&1
