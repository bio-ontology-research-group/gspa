#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --exclude=nid00022,nid00497,nid01033
#SBATCH --time=12:00:00
#SBATCH --no-requeue
#SBATCH --job-name=gspa-price-stage2b-evaluate
#SBATCH --output=/scratch/project/k10018/gspa-v3/out/price-stage2b-rc1/slurm/evaluate-%j.log
# Score stage0 vs stage2b of the verified Stage 2B collection on the 148 Price proteins (same audit as the rc1 scoring).
set -euo pipefail
O=/scratch/project/k10018/gspa-v3/out/price-stage2b-rc1
REPO=/project/k10018/code/gspa-v3
AUDIT=/scratch/project/k10018/gspa-v3/out/price-reference-rc1/ec-audit.json
OUT=$O/evaluation-$SLURM_JOB_ID
module load singularity
mkdir -p "$OUT"
[ "$(sha256sum $AUDIT | cut -d' ' -f1)" = 5302494b16ecd3f76fd96d24e48e27c298206b27608be72c5c66f1d197f4f79f ]
sha256sum $O/collection/collection-verification.json $O/code/stage2b/evaluate_price_stage2b.py $AUDIT $REPO/benchmark/manifests/suite6_price149_manifest.tsv > "$OUT/inputs.sha256"
for P in ${PREDICTORS:-topec mapred enzbert dpz graphec clean}; do
  singularity exec --no-home --no-mount tmp \
    --bind /project/k10018/sw/envs/gspa-tools:/project/k10018/sw/envs/gspa-tools:ro --bind /project/k10018/sw/envs/gspa-v3:/project/k10018/sw/envs/gspa-v3:ro \
    --bind /project/k10018/gspa-v3:/data/hohndor/gspa-v3:ro --bind /project/k10018/code:/project/k10018/code:ro \
    --bind /scratch/project/k10018/gspa-v3/out:/scratch/project/k10018/gspa-v3/out --bind /scratch/project/k10018/gspa-v3/shared/arm-fs-d04c9e11ee64db47/usr/lib/jvm/java-21-openjdk-amd64:/usr/lib/jvm/java-21-openjdk-amd64:ro --bind /scratch/project/k10018/gspa-v3/shared/arm-fs-d04c9e11ee64db47/home/leechuck:/home/leechuck:ro --bind /usr/bin/time:/usr/bin/time:ro --bind /scratch/project/k10018/gspa-v3/work:/scratch/project/k10018/gspa-v3/work:ro \
    --env PYTHONDONTWRITEBYTECODE=1 \
    --env PYTHONPATH=$O/code/stage2b:$O/code:/data/hohndor/gspa-v3/price-production-v1/venv/lib/python3.12/site-packages \
    /project/k10018/sw/images/debian-trixie \
    /project/k10018/sw/envs/gspa-tools/bin/python3.12 $O/code/stage2b/evaluate_price_stage2b.py \
    --collection $O/collection --manifest $REPO/benchmark/manifests/suite6_price149_manifest.tsv \
    --ec-audit $AUDIT --predictors $P --out "$OUT/report-$P" > "$OUT/evaluation-$P.log" 2>&1 &
done
fail=0; for j in $(jobs -p); do wait $j || fail=$((fail+1)); done
echo "failed predictors: $fail"; [ $fail -eq 0 ]
echo "exit 0" > "$OUT/exit-code.txt"
