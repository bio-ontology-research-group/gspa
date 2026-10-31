#!/bin/bash
# WP5 masked assignment recovery on unimatrix01 (design tasks/v3-wp5-masked-recovery-design-2026-10-29.md §8-§9).
# Prepared by Task N; NOT submitted. Same recipe as nutrient_community.sbatch: an immutable package directory is
# copied to node-local scratch, verified against SHA256SUMS, run, and the whole work tree is copied back to
# <package>/jobs/<jobid>/ (also on failure).
#
#   sbatch benchmark/v3/wp5_recovery.sbatch /data/hohndor/gspa-v3/wp5-recovery-v1/package
#
# Package layout (all read-only; every file listed in SHA256SUMS):
#   benchmark/ scripts/                         source snapshot (the WP5 scripts + scripts/metabolic_experiments.py,
#                                               scripts/masked_assignment_experiments.py)
#   inputs/keio/cases_manuscript_regulatory.jsonl   prepare_metabolic_benchmarks.py derived Keio regulatory cases
#   inputs/nutrient/{inputs,observations}.jsonl.gz  benchmark/v3/results/nutrient-environments-v1/data/derived-v3/
#   inputs/sulfonate/{inputs,observations}.jsonl    benchmark/v3/results/sulfonate-diagnostic-v1/dataset/
#   inputs/models/iML1515.xml inputs/models/iYO844.xml   registered models (register_nutrient_models.py)
#   inputs/proteomes/Ec.faa inputs/proteomes/Bs.faa      NCBI translated-CDS FASTA (MG1655 NC_000913.3; 168 AL009126.3)
#   inputs/swissprot/uniprot_sprot.dat.gz                frozen Swiss-Prot release (hash recorded by the reference stage)
#   inputs/gspa-3.0.0.jar                                production JAR (hash-bound; design §9 names build 10526)
#
# Environment: STAGES (default below; "none" runs only the WP5 unit tests), JAVA, DIAMOND (binaries; must exist),
# QWEN_SHARDS (default 4 concurrent Qwen workers), COBRA_REST (default 0: replay arm-visited cells only; 0.1 = the
# §9 seeded 10 % sample of the remaining exhaustive cells). The Qwen key is read by metabolic_experiments.load_api_key()
# from BORG_LLM_API_KEY / OPENAI_API_KEY or ~/.config/borg-llm/api-key; never put it here. GSPA_LLM_URL defaults to
# the unimatrix01 endpoint. Do not run on node005/node008 (Qwen services) — excluded below.
# The audit is a gate: a failed masking audit stops the job before any LP (an invalid run, not a negative result).
# join/calibrate/evaluate/report (evaluate_wp5_recovery.py) are deliberately NOT part of this job: the policy is
# frozen and its SHA-256 reviewed before test labels are joined.
#SBATCH --job-name=gspa-wp5-recovery
#SBATCH --partition=debug
#SBATCH --exclude=node001,node005,node006,node007,node008,unimatrix
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --signal=B:TERM@300
set -euo pipefail
WP5_PACKAGE=${1:?immutable package directory required}
WP5_PY=/data/hohndor/gspa-v3/tools/metabolic-venv/bin/python
WP5_WORK=/scratch/$SLURM_JOB_ID/wp5-recovery
WP5_JAVA=${JAVA:-java}
WP5_DIAMOND=${DIAMOND:-diamond}
WP5_STAGES=${STAGES:-prepare reference query diamond claims construct audit protocol exhaustive deterministic qwen transcripts cobra}
WP5_QWEN_SHARDS=${QWEN_SHARDS:-4}
WP5_COBRA_REST=${COBRA_REST:-0}
WP5_PIDS=()
mkdir -p "$WP5_WORK"
cp -a "$WP5_PACKAGE/benchmark" "$WP5_PACKAGE/scripts" "$WP5_PACKAGE/inputs" "$WP5_PACKAGE/SHA256SUMS" "$WP5_WORK/"
# Optional resume: <package>/resume/{bundle,runs} from an earlier run of the same frozen source (e.g. the Shaheen job
# that ran prepare..deterministic), every file hashed in resume/SHA256SUMS with paths relative to the work tree.
# Run only the later stages then (STAGES="qwen transcripts"); prepare refuses to overwrite an existing bundle.
if [ -d "$WP5_PACKAGE/resume" ]; then
    cp -a "$WP5_PACKAGE/resume/bundle" "$WP5_PACKAGE/resume/runs" "$WP5_WORK/"
    cp "$WP5_PACKAGE/resume/SHA256SUMS" "$WP5_WORK/RESUME-SHA256SUMS"
fi
cd "$WP5_WORK"
export TMPDIR=/scratch/$SLURM_JOB_ID
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPATH=benchmark/v3
export GSPA_LLM_URL=${GSPA_LLM_URL:-http://unimatrix01.kaust.edu.sa:8000/v1}
finish_wp5_job() {
    result_code=$?
    trap - EXIT
    printf '%s\n' "$result_code" > exit-code.txt
    date --iso-8601=seconds > finished-at.txt
    mkdir -p "$WP5_PACKAGE/jobs/$SLURM_JOB_ID"
    cp -a . "$WP5_PACKAGE/jobs/$SLURM_JOB_ID/"
    exit "$result_code"
}
stop_wp5_job() {
    for pid in "${WP5_PIDS[@]}"; do kill -TERM "$pid" 2>/dev/null || true; done
    for pid in "${WP5_PIDS[@]}"; do wait "$pid" 2>/dev/null || true; done
    exit 143
}
trap finish_wp5_job EXIT
trap stop_wp5_job TERM
date --iso-8601=seconds > started-at.txt
hostname > hostname.txt
sha256sum -c SHA256SUMS > input-verification.log
[ -f RESUME-SHA256SUMS ] && sha256sum -c --quiet RESUME-SHA256SUMS > resume-verification.log
"$WP5_PY" - <<'PY' > packages.txt
from importlib.metadata import distributions
print('\n'.join(sorted(d.metadata['Name'] + '==' + d.version for d in distributions())))
PY
# The WP5 unit tests (synthetic fixtures and mocks only) must pass on this interpreter before any stage runs.
"$WP5_PY" -m unittest discover -s benchmark/v3/tests -p 'test_wp5_*.py' > tests.log 2>&1
[ "$WP5_STAGES" = none ] && { echo smoke-ok > smoke.txt; exit 0; }

run_parallel() {   # run_parallel <n> <command...>: n shard workers "i/n" appended as --shard i/n
    local n=$1; shift
    WP5_PIDS=()
    for i in $(seq 0 $((n - 1))); do
        "$@" --shard "$i/$n" > "logs/$(basename "$2" .py)-$3-shard$i.log" 2>&1 &
        WP5_PIDS+=($!)
    done
    local status=0
    for pid in "${WP5_PIDS[@]}"; do wait "$pid" || status=1; done
    WP5_PIDS=()
    return $status
}
mkdir -p logs oracle
B=benchmark/v3
for stage in $WP5_STAGES; do
    date --iso-8601=seconds > "logs/$stage.started"
    case "$stage" in
    prepare)
        [ -e bundle ] && { echo "bundle exists (resume); refusing to overwrite" >&2; exit 2; }
        "$WP5_PY" $B/prepare_wp5_cases.py --out bundle \
            --keio-cases inputs/keio/cases_manuscript_regulatory.jsonl --keio-model inputs/models/iML1515.xml \
            --ec-model inputs/models/iML1515.xml --bs-model inputs/models/iYO844.xml \
            --nutrient-inputs inputs/nutrient/inputs.jsonl.gz --nutrient-observations inputs/nutrient/observations.jsonl.gz \
            --sulfonate-inputs inputs/sulfonate/inputs.jsonl --sulfonate-observations inputs/sulfonate/observations.jsonl \
            > logs/prepare.log 2>&1 ;;
    reference)
        "$WP5_PY" $B/build_wp5_candidates.py reference --bundle bundle --swissprot inputs/swissprot/uniprot_sprot.dat.gz \
            > logs/reference.log 2>&1 ;;
    query)
        "$WP5_PY" $B/build_wp5_candidates.py query --bundle bundle \
            --proteome Ec=inputs/proteomes/Ec.faa --proteome Bs=inputs/proteomes/Bs.faa > logs/query.log 2>&1 ;;
    diamond)
        command -v "$WP5_DIAMOND" > /dev/null || { echo "DIAMOND binary not found: $WP5_DIAMOND" >&2; exit 2; }
        "$WP5_PY" $B/build_wp5_candidates.py diamond --bundle bundle --threads "$SLURM_CPUS_PER_TASK" \
            --diamond-bin "$WP5_DIAMOND" > logs/diamond.log 2>&1 ;;
    claims)
        "$WP5_PY" $B/build_wp5_candidates.py claims --bundle bundle > logs/claims.log 2>&1 ;;
    construct)
        command -v "$WP5_JAVA" > /dev/null || { echo "java not found: $WP5_JAVA" >&2; exit 2; }
        "$WP5_JAVA" -version > logs/java-version.txt 2>&1
        "$WP5_PY" $B/build_wp5_candidates.py construct --bundle bundle --jar inputs/gspa-3.0.0.jar \
            --java "$WP5_JAVA" --python "$WP5_PY" --sidecar scripts/metabolic_experiments.py > logs/construct.log 2>&1 ;;
    audit)
        # Gate: exit status 1 (failed audit) aborts the job here.
        "$WP5_PY" $B/audit_wp5_masking.py all --bundle bundle --jar inputs/gspa-3.0.0.jar --java "$WP5_JAVA" \
            --python "$WP5_PY" --sidecar scripts/metabolic_experiments.py --scratch "$TMPDIR/wp5-audit" > logs/audit.log 2>&1 ;;
    protocol)
        "$WP5_PY" $B/run_wp5_recovery.py protocol --public bundle/public --out runs > logs/protocol.log 2>&1 ;;
    exhaustive|deterministic)
        run_parallel "$SLURM_CPUS_PER_TASK" "$WP5_PY" $B/run_wp5_recovery.py "$stage" --public bundle/public --out runs ;;
    qwen)
        run_parallel "$WP5_QWEN_SHARDS" "$WP5_PY" $B/run_wp5_recovery.py qwen --public bundle/public --out runs ;;
    transcripts)
        "$WP5_PY" $B/audit_wp5_masking.py transcripts --bundle bundle --runs runs > logs/transcripts.log 2>&1 ;;
    cobra)
        "$WP5_PY" $B/check_wp5_cobra.py --public bundle/public --runs runs --out oracle/cobra.json \
            --rest-sample "$WP5_COBRA_REST" > logs/cobra.log 2>&1 ;;
    *)
        echo "unknown stage: $stage" >&2; exit 2 ;;
    esac
    date --iso-8601=seconds > "logs/$stage.finished"
done
echo done > done.txt
