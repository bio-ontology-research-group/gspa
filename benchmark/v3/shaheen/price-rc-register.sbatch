#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --job-name=gspa-price-rc-register
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --no-requeue
#SBATCH --exclude=nid00022
#SBATCH --output=/scratch/project/k10018/gspa-v3/out/price-production-rc1/slurm/register-%j.log
# Fresh registration of the full 132-cell Price production collection for a release-candidate JAR on Shaheen III.
# Equivalent of unimatrix price-production-v3/prepare.sbatch (job 10292) without its pilot: the same frozen engine
# code (verified against the original code-hashes.json), the same export plan, audit, layout, parts and reference,
# only the JAR (a parameter; its hash is recorded by the registration) and the run_price_local.py wrapper (adds the
# opt-in external launch used on Shaheen) differ. The frozen unimatrix collection is never read or written: its
# inputs come from the read-only /project mirror bound at their original /data paths.
# Writes: $OUT_ROOT/{code,gspa-3.0.0.jar,collection/{plan,registration}.json,register-receipt.json,jobs/<job>/}.
# Parameters (environment): OUT_ROOT, JAR, EXPECT_JAR_SHA256 / EXPECT_WRAPPER_SHA256 (optional), FROZEN, WRAPPER, WRAPPER_TEST, PYTHON,
# PYLIB, FSROOT, IMAGE. Defaults below; price-rc-submit.sh passes them explicitly.
set -euo pipefail
module load singularity
unset JAVA_TOOL_OPTIONS JDK_JAVA_OPTIONS _JAVA_OPTIONS PYTHONPATH PYTHONHOME
OUT_ROOT=${OUT_ROOT:-/scratch/project/k10018/gspa-v3/out/price-production-rc1}
JAR=${JAR:-/project/k10018/gspa-v3/release-candidate-1/gspa-3.0.0.jar}
FROZEN=${FROZEN:-/project/k10018/gspa-v3/price-production-v3}
WRAPPER=${WRAPPER:-/project/k10018/code/gspa-v3/benchmark/v3/run_price_local.py}
WRAPPER_TEST=${WRAPPER_TEST:-/project/k10018/code/gspa-v3/benchmark/v3/tests/test_price_local_external.py}
PYTHON=${PYTHON:-/project/k10018/sw/envs/gspa-tools/bin/python3.12}
PYLIB=${PYLIB:-/data/hohndor/gspa-v3/price-production-v1/venv/lib/python3.12/site-packages}
FSROOT=${FSROOT:-/scratch/project/k10018/gspa-v3/shared/arm-fs-d04c9e11ee64db47}
IMAGE=${IMAGE:-/project/k10018/sw/images/debian-trixie}
WORK_ROOT=${WORK_ROOT:-/scratch/project/k10018/gspa-v3/work}
logical_data=/data/hohndor/gspa-v3
# sha256 of the original price-production-v3/code-hashes.json (unimatrix, 2026-10-28); it lists the frozen engine.
FROZEN_CODE_HASHES_SHA256=cfb8ac69520d2f1df07849cd0825ffe0224e97b9f29a275cde1779490968a7c4
# Fingerprints of the original plan (2407147c…) and registration (454f466c…): sha256 of json.dumps(x, sort_keys=True)
# of plan['cells'], of the 162 data inputs (parts, layout, export plan, audit, reference) and of the 17 parser files.
ORIGINAL_CELLS=4bc45a18801dc7bc2d1ef4df46aad44f82fd267be12d2312d6ebf5a41aa486bb
ORIGINAL_DATA_INPUTS=bcadd3959fd0ba06896df6eef4c88ebeb8f446ae6d51b295e5a88fb1aa989937
ORIGINAL_PARSER_INPUTS=0142bccd6e5ecac1411aaa5385acf375640884fd9e00485b8d537b72babd93be
ENGINE=(run_price_matrix_capacity.py run_price_production_collection.py stage_ledger.py v3_metrics.py v3_runner.py
        verified_batch.py verify_candidate_universe.py requirements-streaming.txt)
# sha256 of the original price-production-v3/test-source-hashes.json (133 files of the frozen test bundle).
FROZEN_TEST_HASHES_SHA256=fe33c177067ce9713aaa231421782ea165fb8fc6eeb29b03a715a318901b45fd

jobout=$OUT_ROOT/jobs/$SLURM_JOB_ID
work=$WORK_ROOT/$SLURM_JOB_ID
[ ! -e "$OUT_ROOT/collection" ] || { echo "refusing: $OUT_ROOT/collection already exists" >&2; exit 2; }
[ ! -e "$OUT_ROOT/code" ] || { echo "refusing: $OUT_ROOT/code already exists" >&2; exit 2; }
[ -f "$FSROOT/.complete" ] || { echo "shared package tree missing: $FSROOT" >&2; exit 2; }
mkdir -p "$jobout" "$work/tmp" "$OUT_ROOT/code/tests"
finish() { code=$?; trap - EXIT; printf '%s\n' "$code" > "$jobout/exit-code.txt"; date --iso-8601=seconds > "$jobout/finished-at.txt"; exit "$code"; }
trap finish EXIT
date --iso-8601=seconds > "$jobout/started-at.txt"; hostname > "$jobout/hostname.txt"
cp "$0" "$jobout/executed.sbatch"
singularity --version > "$jobout/singularity-version.txt"

# Snapshot: frozen engine files byte-for-byte, the candidate wrapper, the JAR. Hash checks run in the container below.
echo "$FROZEN_CODE_HASHES_SHA256  $FROZEN/code-hashes.json" | sha256sum -c - > "$jobout/frozen-code-check.log"
for f in "${ENGINE[@]}"; do cp "$FROZEN/code/$f" "$OUT_ROOT/code/$f"; done
echo "$FROZEN_TEST_HASHES_SHA256  $FROZEN/test-source-hashes.json" | sha256sum -c - >> "$jobout/frozen-code-check.log"
cp "$FROZEN/code/tests/test_price_local.py" "$OUT_ROOT/code/tests/test_price_local.py"
cp "$FROZEN/code-hashes.json" "$OUT_ROOT/code/frozen-code-hashes.json"
cp "$WRAPPER" "$OUT_ROOT/code/run_price_local.py"
cp "$WRAPPER_TEST" "$OUT_ROOT/code/tests/test_price_local_external.py"
cp "$JAR" "$OUT_ROOT/gspa-3.0.0.jar"
chmod a-w "$OUT_ROOT"/code/* "$OUT_ROOT"/code/tests/* "$OUT_ROOT/gspa-3.0.0.jar"
( cd "$OUT_ROOT" && sha256sum gspa-3.0.0.jar code/*.py code/*.txt code/*.json code/tests/*.py ) > "$jobout/snapshot.sha256"
if [ -n "${EXPECT_JAR_SHA256:-}" ]; then echo "$EXPECT_JAR_SHA256  $OUT_ROOT/gspa-3.0.0.jar" | sha256sum -c - >> "$jobout/snapshot-check.log"; fi
if [ -n "${EXPECT_WRAPPER_SHA256:-}" ]; then echo "$EXPECT_WRAPPER_SHA256  $OUT_ROOT/code/run_price_local.py" | sha256sum -c - >> "$jobout/snapshot-check.log"; fi

container=(singularity exec --no-home --no-mount tmp --pwd "$OUT_ROOT"
  --bind /project/k10018/sw/envs/gspa-tools:/project/k10018/sw/envs/gspa-tools:ro
  --bind "/project/k10018/gspa-v3:$logical_data:ro" --bind "$OUT_ROOT:$OUT_ROOT" --bind "$work:$work"
  --bind "$FSROOT/usr/lib/jvm/java-21-openjdk-amd64:/usr/lib/jvm/java-21-openjdk-amd64:ro" --bind "$FSROOT/home/leechuck:/home/leechuck:ro"
  --bind /usr/bin/time:/usr/bin/time:ro
  --env PATH=/usr/lib/jvm/java-21-openjdk-amd64/bin:/usr/bin:/bin --env LC_ALL=C.UTF-8 --env PYTHONDONTWRITEBYTECODE=1
  --env "PYTHONPATH=$PYLIB" --env "TMPDIR=$work/tmp"
  "$IMAGE")
printf '%s\n' "${container[@]}" > "$jobout/container-command.txt"
# Tools the cells need inside the image: GNU time (the capacity runner wraps java in /usr/bin/time -v), sort, java.
"${container[@]}" /usr/bin/time -v -o "$work/tmp/time-probe.txt" java -version > "$jobout/java-version.txt" 2>&1
grep -q 'Maximum resident set size' "$work/tmp/time-probe.txt"
"${container[@]}" sort --version > "$jobout/sort-version.txt"

export OUT_ROOT FROZEN_LOGICAL=$logical_data/${FROZEN#/project/k10018/gspa-v3/} ORIGINAL_CELLS ORIGINAL_DATA_INPUTS ORIGINAL_PARSER_INPUTS PYLIB SLURM_JOB_ID
"${container[@]}" "$PYTHON" - > "$jobout/preflight.log" 2>&1 <<'PY'
import hashlib, json, os, sys
from pathlib import Path
out = Path(os.environ['OUT_ROOT'])
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
frozen = json.loads((out/'code/frozen-code-hashes.json').read_text())
for relative, digest in frozen.items():   # tests/* other than test_price_local.py live in the frozen test bundle
    if relative != 'run_price_local.py' and (relative.count('/') == 0 or relative == 'tests/test_price_local.py'):
        assert sha(out/'code'/relative) == digest, relative
assert sorted(p.name for p in (out/'code').glob('*.py')) == sorted(k for k in frozen if k.endswith('.py') and '/' not in k)
tests = Path(os.environ['FROZEN_LOGICAL'])
for relative, digest in json.loads((tests/'test-source-hashes.json').read_text()).items():
    assert sha(tests/'test-source'/relative) == digest, relative
import ijson
assert ijson.backend == 'yajl2_c', ('a native ijson backend is required for multi-GB ledgers', ijson.backend)
assert ijson.__version__ == '3.4.0.post0', ijson.__version__
print(json.dumps(dict(python=sys.version, executable=sys.executable, ijson=ijson.__file__, backend=ijson.backend,
                      wrapper_sha256=sha(out/'code/run_price_local.py')), indent=2))
PY

# Tests, without skips: (1) the original prepare.sbatch set (64 tests of the frozen test bundle, job 10292);
# (2) the wrapper actually used here (frozen test_price_local.py plus the external-launch tests) on the snapshot.
# The tests' temporary files go to node memory: files on Shaheen's Lustre report st_blksize 4 MiB, so Python's default
# file buffers alone exceed test_price_capacity's 3 MiB tracemalloc ceiling (8,463,296 B in job 15707970); /dev/shm
# reports 4 KiB, as on the hosts where the frozen tests were written and passed (unimatrix 10292, ws).
testtmp=/dev/shm/price-rc-tests-$SLURM_JOB_ID; mkdir -p "$testtmp"
"${container[@]}" env TMPDIR="$testtmp" "$PYTHON" - > "$jobout/tests.log" 2>&1 <<'PY'
import os, sys, unittest
frozen = os.environ['FROZEN_LOGICAL']+'/test-source/benchmark/v3'; code = os.environ['OUT_ROOT']+'/code'
def run(paths, names, expected=None):
    saved = list(sys.path); sys.path[:0] = paths
    for name in [m for m in sys.modules if m.startswith(('test_', 'run_price', 'verified_batch', 'stage_ledger',
                                                          'v3_', 'verify_candidate'))]:
        del sys.modules[name]
    suite = unittest.defaultTestLoader.loadTestsFromNames(names)
    result = unittest.TextTestRunner(verbosity=2).run(suite); sys.path[:] = saved
    print(names, 'run', result.testsRun, 'failures', len(result.failures), 'errors', len(result.errors), 'skipped', len(result.skipped))
    assert result.wasSuccessful() and not result.skipped and (expected is None or result.testsRun == expected)
    return result.testsRun
run([frozen, frozen+'/tests'], ['test_price_local', 'test_price_capacity', 'test_price_production_collection',
                                'test_verified_batch', 'test_streaming_ledger'], 64)
run([code, code+'/tests'], ['test_price_local', 'test_price_local_external'])
PY

/usr/bin/time -v -o "$jobout/prepare.resources.txt" "${container[@]}" "$PYTHON" "$OUT_ROOT/code/run_price_production_collection.py" prepare \
  --out "$OUT_ROOT/collection" --export "$logical_data/price-full-matrix-v2/plan.json" \
  --audit "$logical_data/price-full-collection-audit-v1/verification.json" --jar "$OUT_ROOT/gspa-3.0.0.jar" \
  --reference "$logical_data/price-production-v1/reference" > "$jobout/prepare.log" 2>&1

# Same plan cells and same data inputs as the original registration; record what differs (JAR, runtimes, wrapper).
"${container[@]}" "$PYTHON" - > "$jobout/compare-original.log" 2>&1 <<'PY'
import hashlib, json, os
from pathlib import Path
out = Path(os.environ['OUT_ROOT']); data = '/data/hohndor/gspa-v3/'
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
digest = lambda x: hashlib.sha256(json.dumps(x, sort_keys=True).encode()).hexdigest()
plan = json.loads((out/'collection/plan.json').read_text())
registration = json.loads((out/'collection/registration.json').read_text())
inputs = registration['inputs']
assert len(plan['cells']) == 132 and plan['expected_claims'] == 1526919081
assert digest(plan['cells']) == os.environ['ORIGINAL_CELLS'], 'plan cells differ from the original collection'
prefixes = tuple(data+p for p in ('price-full-matrix-v2/', 'price-full-collection-audit-v1/', 'price-production-v1/reference/'))
assert digest({k: v for k, v in inputs.items() if k.startswith(prefixes)}) == os.environ['ORIGINAL_DATA_INPUTS'], 'data inputs differ'
parser = {k: v for k, v in inputs.items() if k.startswith(data+'price-production-v1/venv/')}
same_parser = digest(parser) == os.environ['ORIGINAL_PARSER_INPUTS']
if os.environ['PYLIB'] == data+'price-production-v1/venv/lib/python3.12/site-packages':
    assert same_parser, 'parser bytes differ from the original registration'
jar = str(out/'gspa-3.0.0.jar')
assert plan['jar'] == jar and inputs.get(jar) == sha(jar), 'registration does not bind the snapshot JAR'
receipt = dict(schema='gspa.price-rc-registration/1', job=os.environ['SLURM_JOB_ID'],
    plan_sha256=registration['plan_sha256'], registration_sha256=sha(out/'collection/registration.json'),
    jar=plan['jar'], jar_sha256=inputs.get(jar), java=plan['java'], java_sha256=inputs.get(plan['java']), python=plan['python'],
    wrapper_sha256=inputs.get(str(out/'code/run_price_local.py')), original_plan_cells_identical=True,
    original_data_inputs_identical=True, original_parser_identical=same_parser,
    original=dict(plan_sha256='2407147c51dd0fe25cea0ea8193feefcac75d8cdbe46f0009d9001474339c4a4',
                  registration_sha256='454f466cb635b70a9bbb65cbb70d1022065568d1931ac944be95e27a754b02dd',
                  jar_sha256='0606bec0bf009bc6b9c1301a28ae8d21998730c36138c24cb5720a44eb31a50a'),
    interpretation='Fresh registration of the unchanged 132-cell Price plan for a candidate JAR; no cell has run.')
with (out/'register-receipt.json').open('x') as stream:
    json.dump(receipt, stream, indent=2, sort_keys=True); stream.write('\n')
print(json.dumps(receipt, indent=2, sort_keys=True))
PY
cp "$OUT_ROOT/register-receipt.json" "$jobout/"
rm -rf "$work"
echo "registered: $(cat "$OUT_ROOT/collection/registration.json" | sha256sum | cut -c1-64)"
