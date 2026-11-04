#!/bin/bash
# Generate and print (default) or submit (--submit) the Shaheen III SLURM jobs for row5-masked-context.
#
# Usage:
#   row5-masked-context-submit.sh [--loci PATH] [--claim-sets PATH] [--max-nodes N] [--pilot] [-- tokens-out PATH] [--submit]
#
# Rules:
#   - Maskable loci from loci.tsv (or loci.json) x claim sets of their group
#     (Price: all 7 claim sets; cohort: diamond only).
#   - --pilot selects the 5 pilot loci of plan A2 in hash order of locus_id (2 cohort, 3 Price) and all arms.
#   - Packs <locus_id>:<claim_set> tokens into jobs of roughly equal total published claim count (--max-nodes).
#   - Always passes --exclude=nid00022,nid00497,nid01033 and explicit REPO_ROOT.
#   - Never queues more than 400 jobs.
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
REPO_ROOT=${REPO_ROOT:-$(cd "$here/../../.." && pwd)}
OUT_ROOT=${OUT_ROOT:-/scratch/project/k10018/gspa-v3/out/row5-masked-context-v1}
LOCI=${LOCI:-$OUT_ROOT/registry/loci.tsv}
CLAIM_SETS=${CLAIM_SETS:-$OUT_ROOT/registry/claim_sets.json}
PYTHON_BIN=${PYTHON_BIN:-$(command -v python3)}
MAX_NODES=""
PILOT=0
SUBMIT=0
INCLUDE_STAGE=${INCLUDE_STAGE:-0}
INCLUDE_EVAL=${INCLUDE_EVAL:-0}
TOKENS_OUT=""
CELL_TIME=${CELL_TIME:-12:00:00}
STAGE_TIME=${STAGE_TIME:-02:00:00}
EVAL_TIME=${EVAL_TIME:-01:00:00}
ARMS=${ARMS:-}

while [ "$#" -gt 0 ]; do
  case "$1" in
    --loci) LOCI=$2; shift 2 ;;
    --claim-sets) CLAIM_SETS=$2; shift 2 ;;
    --max-nodes) MAX_NODES=$2; shift 2 ;;
    --pilot) PILOT=1; shift ;;
    --submit) SUBMIT=1; shift ;;
    --include-stage) INCLUDE_STAGE=1; shift ;;
    --include-evaluate) INCLUDE_EVAL=1; shift ;;
    --tokens-out) TOKENS_OUT=$2; shift 2 ;;
    --arms) ARMS=$2; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

if [ -z "$MAX_NODES" ]; then
  if [ "$PILOT" -eq 1 ]; then
    MAX_NODES=1
  else
    MAX_NODES=16
  fi
fi

if [ "$MAX_NODES" -lt 1 ] || [ "$MAX_NODES" -gt 400 ]; then
  echo "ERROR: --max-nodes must be between 1 and 400 (got $MAX_NODES)" >&2
  exit 2
fi

common_flags="--account=k10018 --partition=workq --nodes=1 --exclude=nid00022,nid00497,nid01033 --no-requeue"
pass_vars=(
  REPO_ROOT OUT_ROOT RC1_PLAN PRICE_LAYOUT COHORT_LAYOUT_DIR COHORT_BATCH_MANIFEST
  COHORT_LABELS COHORT_SPLITS PRICE_SPLIT PRICE148_MANIFEST PRICE_SWISSPROT
  COHORT_ORGANISM_METADATA COHORT_METADATA COHORT_ORTHOGROUPS PRICE_COHORT_DIR
  PRICE_DIAMOND_DIR COHORT_DIAMOND_DIR MODELS_TAR COHORT_MODELS_TAR REFERENCE
  JAR EXPECT_JAR_SHA256 JAR_CTX EXPECT_JAR_CTX_SHA256 MODEL GSREF PATHWAY_GROUPS
  PYTHON PYLIB FSROOT IMAGE WORK_ROOT MEM_BUDGET_G HEAP_SMALL HEAP_MEDIUM HEAP_LARGE
  ALLOC_SMALL_G ALLOC_MEDIUM_G ALLOC_LARGE_G TIMEOUT_SMALL TIMEOUT_MEDIUM TIMEOUT_LARGE
  WINDOW_BP KEEP_OUTPUTS
 CELLS_DIR
  SMALL_MAX_CLAIMS MEDIUM_MAX_CLAIMS)
envs=("REPO_ROOT=$REPO_ROOT" "OUT_ROOT=$OUT_ROOT")
for v in "${pass_vars[@]}"; do
  if [ "$v" = "REPO_ROOT" ] || [ "$v" = "OUT_ROOT" ]; then
    continue
  fi
  [ -z "${!v:-}" ] || envs+=("$v=${!v}")
done

# Pack tokens into <= MAX_NODES balanced bins by claim count using Python
packed_json=$("$PYTHON_BIN" - "$LOCI" "$CLAIM_SETS" "$MAX_NODES" "$PILOT" "$TOKENS_OUT" "$OUT_ROOT" "${SKIP_EXISTING:-0}" "${CELLS_DIR:-}" <<'PY'
import csv
import heapq
import json
import sys
from pathlib import Path

loci_path = Path(sys.argv[1])
csets_path = Path(sys.argv[2])
max_nodes = int(sys.argv[3])
pilot = int(sys.argv[4]) == 1
tokens_out = sys.argv[5]
out_root = sys.argv[6] if len(sys.argv) > 6 else ""
skip_existing = len(sys.argv) > 7 and sys.argv[7] == "1"
cells_dir = sys.argv[8] if len(sys.argv) > 8 and sys.argv[8] else str(Path(out_root) / "cells")

PRICE_CLAIM_SETS = ("clean", "diamond", "dpz", "enzbert", "graphec", "mapred", "topec")

maskable = []
if loci_path.suffix == ".json":
    data = json.loads(loci_path.read_text(encoding="utf-8"))
    rows = data.get("all_loci") or list(data.get("loci", {}).values())
    for r in rows:
        if r.get("maskable"):
            maskable.append({
                "locus_id": str(r["locus_id"]),
                "genome": str(r["genome"]),
                "group": str(r["group"]),
            })
else:
    with loci_path.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if str(r.get("maskable", "")).strip().lower() == "true":
                maskable.append({
                    "locus_id": r["locus_id"].strip(),
                    "genome": r["genome"].strip(),
                    "group": r["group"].strip(),
                })

# Sort by locus_id (hash order of locus_id)
maskable.sort(key=lambda x: x["locus_id"])

if pilot:
    cohort_pilot = [x for x in maskable if x["group"] == "cohort"][:2]
    price_pilot = [x for x in maskable if x["group"] == "price"][:3]
    selected_loci = sorted(cohort_pilot + price_pilot, key=lambda x: x["locus_id"])
else:
    selected_loci = maskable

claim_counts = {}
csets_by_genome = {}
if csets_path.is_file():
    cdata = json.loads(csets_path.read_text(encoding="utf-8"))
    for e in cdata.get("entries", []):
        cs = e["claim_set"]
        g = e["genome"]
        claim_counts[(cs, g)] = int(e.get("published_claims") or 50000)
        csets_by_genome.setdefault(g, []).append(cs)

tokens = []
for loc in selected_loci:
    lid = loc["locus_id"]
    g = loc["genome"]
    grp = loc["group"]
    if grp == "price":
        csets = sorted(set(csets_by_genome.get(g, PRICE_CLAIM_SETS)))
    else:
        csets = ["diamond"]
    for cs in csets:
        default_w = 350000 if cs in ("clean", "dpz", "graphec") else 60000
        w = claim_counts.get((cs, g), default_w)
        # SKIP_EXISTING=1: leave out cells that already hold a receipt (for example the pilot cells)
        if skip_existing and (Path(cells_dir) / lid / cs / "receipt.json").is_file():
            continue
        tokens.append({"token": f"{lid}:{cs}", "claims": w, "group": grp, "claim_set": cs})

if tokens_out:
    out_p = Path(tokens_out)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    out_p.write_text("\n".join(t["token"] for t in tokens) + ("\n" if tokens else ""), encoding="utf-8")

n_bins = min(max_nodes, max(1, len(tokens)))
# Heap of (total_claims, bin_index)
heap = [(0, i) for i in range(n_bins)]
heapq.heapify(heap)
bins = [[] for _ in range(n_bins)]
bin_weights = [0] * n_bins

# Greedy LPT packing (largest claim count first, deterministic tie-break by token)
for item in sorted(tokens, key=lambda x: (-x["claims"], x["token"])):
    total_w, b_idx = heapq.heappop(heap)
    bins[b_idx].append(item["token"])
    new_w = total_w + item["claims"]
    bin_weights[b_idx] = new_w
    heapq.heappush(heap, (new_w, b_idx))

non_empty_bins = [
    {"index": i, "tokens": bins[i], "total_claims": bin_weights[i]}
    for i in range(n_bins)
    if bins[i]
]
print(json.dumps({
    "pilot": pilot,
    "loci_count": len(selected_loci),
    "token_count": len(tokens),
    "job_count": len(non_empty_bins),
    "jobs": non_empty_bins,
}))
PY
)

job_count=$("$PYTHON_BIN" -c 'import json,sys; print(json.loads(sys.argv[1])["job_count"])' "$packed_json")
token_count=$("$PYTHON_BIN" -c 'import json,sys; print(json.loads(sys.argv[1])["token_count"])' "$packed_json")
loci_count=$("$PYTHON_BIN" -c 'import json,sys; print(json.loads(sys.argv[1])["loci_count"])' "$packed_json")

extra_jobs=$(( INCLUDE_STAGE + (INCLUDE_EVAL * 2) ))
total_planned_jobs=$(( job_count + extra_jobs ))
if [ "$total_planned_jobs" -gt 400 ]; then
  echo "ERROR: total planned jobs ($total_planned_jobs) exceeds the 400-job queue limit" >&2
  exit 2
fi

q() { [ "$#" -eq 0 ] || printf '%q ' "$@"; }

lines=()
emit() { lines+=("$1"); }

emit "set -euo pipefail"
emit "mkdir -p $(q "$OUT_ROOT/slurm")"

dep_stage=""
if [ "$INCLUDE_STAGE" -eq 1 ]; then
  emit "STAGE_JOB=\$(env $(q "${envs[@]}")MODE=stage sbatch --parsable $common_flags --time=$STAGE_TIME --job-name=gspa-row5-stage --output=$(q "$OUT_ROOT/slurm/stage-%j.log") $(q "$here/row5-masked-context.sbatch"))"
  emit "echo stage=\$STAGE_JOB"
  dep_stage="--dependency=afterok:\$STAGE_JOB "
fi

emit "CELL_JOBS=''"
while IFS=$'\t' read -r b_idx b_claims b_tokens; do
  [ -n "$b_tokens" ] || continue
  arm_env=()
  if [ -n "$ARMS" ]; then
    arm_env=("ARMS=$ARMS")
  elif [ "$PILOT" -eq 1 ]; then
    arm_env=("ARMS=stage0 gaps unmasked_gaps sham_gaps gaps_learned universal pathway organism ctx_universal ctx_pathway ctx_organism sham_ctx_universal gaps_ctxjar")
  fi
  jname="gspa-row5-cells-$b_idx"
  [ "$PILOT" -eq 1 ] && jname="gspa-row5-pilot-$b_idx"
  emit "J=\$(env $(q "${envs[@]}")$(q "${arm_env[@]}")CELLS=$(q "$b_tokens") sbatch --parsable $common_flags --time=$CELL_TIME ${dep_stage}--job-name=$jname --output=$(q "$OUT_ROOT/slurm/cells-%j.log") $(q "$here/row5-masked-context.sbatch")); CELL_JOBS=\$CELL_JOBS:\$J; echo \"job $b_idx (claims=$b_claims) job=\$J\""
done < <("$PYTHON_BIN" -c '
import json, sys
for j in json.loads(sys.argv[1])["jobs"]:
    print(j["index"], j["total_claims"], " ".join(j["tokens"]), sep="\t")
' "$packed_json")

if [ "$INCLUDE_EVAL" -eq 1 ]; then
  emit "EVAL_DEV=\$(env $(q "${envs[@]}")MODE=evaluate SCOPE=development CHOOSE_WINDOW=1 sbatch --parsable $common_flags --time=$EVAL_TIME --dependency=afterok\$CELL_JOBS --job-name=gspa-row5-eval-dev --output=$(q "$OUT_ROOT/slurm/eval-dev-%j.log") $(q "$here/row5-masked-context.sbatch"))"
  emit "EVAL_HELD=\$(env $(q "${envs[@]}")MODE=evaluate SCOPE=heldout sbatch --parsable $common_flags --time=$EVAL_TIME --dependency=afterok:\$EVAL_DEV --job-name=gspa-row5-eval-heldout --output=$(q "$OUT_ROOT/slurm/eval-heldout-%j.log") $(q "$here/row5-masked-context.sbatch"))"
  emit "echo eval_dev=\$EVAL_DEV eval_heldout=\$EVAL_HELD"
fi

if [ "$SUBMIT" -eq 1 ]; then
  command -v sbatch >/dev/null || { echo "sbatch not found: run on a Shaheen login node" >&2; exit 2; }
  if command -v squeue >/dev/null; then
    current_q=$(squeue -u "${USER:-$(id -un)}" -h 2>/dev/null | wc -l || echo 0)
    if [ $((current_q + total_planned_jobs)) -gt 400 ]; then
      echo "ERROR: current queued jobs ($current_q) + planned ($total_planned_jobs) > 400" >&2
      exit 2
    fi
  fi
  printf '%s\n' "${lines[@]}" | bash
else
  echo "# pilot=$PILOT loci=$loci_count tokens=$token_count cell_jobs=$job_count total_jobs=$total_planned_jobs"
  printf '%s\n' "${lines[@]}"
fi
