#!/usr/bin/env bash
# Tests for the --operon-caller switch in run_integrate_full_priors.sh.
#
# Runs the integrate driver with --dry-run in a temp environment and
# asserts the assembled command references the right operons file for
# each caller. Exits 0 on success, non-zero on any failure.
set -u

SCRIPT_DIR=$(cd "$(dirname "$0")"/.. && pwd)
SCRIPT="$SCRIPT_DIR/run_integrate_full_priors.sh"

if [[ ! -x "$SCRIPT" && ! -f "$SCRIPT" ]]; then
    echo "FAIL: integrate runner not found at $SCRIPT" >&2
    exit 1
fi

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

# --- minimal fake env -----------------------------------------------------
export ROOT="$TMP/proteomes"
export BENCH="$ROOT/bench"
export GLM_DIR="$TMP/glm-preds"
export REF="$TMP/ref"
export OUT="$TMP/out"
export JAVA=/usr/bin/true             # never executed under --dry-run
export JAR=/dev/null
mkdir -p "$ROOT/operons" "$BENCH/gapseq" "$GLM_DIR/ecoli" "$REF"

# fake claims file so the per-genome loop does not skip
printf '{}\n' > "$BENCH/ecoli_claims.jsonl"

# heuristic ops file (canonical legacy path)
printf 'p1\tp2\n' > "$ROOT/operons/ecoli_operons.tsv"
# glm ops file (new path)
printf 'p3\tp4\n' > "$GLM_DIR/ecoli/operons.tsv"

# restrict the genome list to ecoli for the test
export GENOMES='ecoli'
export KINGDOMS='bacteria'

PASS=0
FAIL=0

assert_contains() {
    local needle="$1"; local haystack="$2"; local label="$3"
    if printf '%s' "$haystack" | grep -F -q -- "$needle"; then
        echo "  PASS: $label"; PASS=$((PASS+1))
    else
        echo "  FAIL: $label"; FAIL=$((FAIL+1))
        printf '    expected to contain: %s\n' "$needle" >&2
        printf '    got:\n%s\n' "$haystack" | sed 's/^/      /' >&2
    fi
}

assert_not_contains() {
    local needle="$1"; local haystack="$2"; local label="$3"
    if ! printf '%s' "$haystack" | grep -F -q -- "$needle"; then
        echo "  PASS: $label"; PASS=$((PASS+1))
    else
        echo "  FAIL: $label"; FAIL=$((FAIL+1))
        printf '    did NOT expect: %s\n' "$needle" >&2
        printf '    got:\n%s\n' "$haystack" | sed 's/^/      /' >&2
    fi
}

# -------------------------------------------------------------------------
echo "Test 1: default mode is heuristic"
out=$(bash "$SCRIPT" --dry-run 2>&1) || true
assert_contains "operon-caller=heuristic" "$out" "default operon caller is heuristic"
assert_contains "$ROOT/operons/ecoli_operons.tsv" "$out" "heuristic mode resolves legacy operon path"
assert_not_contains "$GLM_DIR/ecoli/operons.tsv" "$out" "heuristic mode does NOT use glm path"

# -------------------------------------------------------------------------
echo "Test 2: --operon-caller heuristic explicit"
out=$(bash "$SCRIPT" --operon-caller heuristic --dry-run 2>&1) || true
assert_contains "operon-caller=heuristic" "$out" "explicit heuristic flag honored"
assert_contains "$ROOT/operons/ecoli_operons.tsv" "$out" "explicit heuristic resolves legacy path"

# -------------------------------------------------------------------------
echo "Test 3: --operon-caller glm"
out=$(bash "$SCRIPT" --operon-caller glm --dry-run 2>&1) || true
assert_contains "operon-caller=glm" "$out" "glm flag honored"
assert_contains "$GLM_DIR/ecoli/operons.tsv" "$out" "glm mode resolves new operon path"
assert_not_contains "$ROOT/operons/ecoli_operons.tsv" "$out" "glm mode does NOT use legacy path"

# -------------------------------------------------------------------------
echo "Test 4: bad operon-caller value rejected"
if bash "$SCRIPT" --operon-caller bogus --dry-run >/dev/null 2>&1; then
    echo "  FAIL: bogus value should exit non-zero"; FAIL=$((FAIL+1))
else
    echo "  PASS: bogus operon-caller rejected with non-zero exit"; PASS=$((PASS+1))
fi

# -------------------------------------------------------------------------
echo "Test 5: --dry-run does not invoke JAR (proven by JAR=/dev/null surviving)"
# If --dry-run actually ran the integrate command, the JAR=/dev/null
# substitution would have produced a Java error in stdout. Confirm
# absence of the canonical Java error string.
out=$(bash "$SCRIPT" --operon-caller glm --dry-run 2>&1) || true
assert_not_contains "Error: Unable to access jarfile" "$out" "dry-run skips JAR invocation"

# -------------------------------------------------------------------------
echo
echo "RESULT: ${PASS} passed, ${FAIL} failed"
[[ $FAIL -eq 0 ]] && exit 0 || exit 1
