# Community numerical audit and interrupted verification

This advances the [registered coculture evaluation](v3-nutrient-community.md).
It does not complete that evaluation or any whole audit requirement. The frozen
plans and production solver remain unchanged; no experimental outcome file is
read by this audit or the numerical diagnosis.

## Receipt and phenotype checks

`audit_nutrient_community.py` verifies the immutable package and each job's
copied code/models/plan, input-check log, production receipt, exact modulo-shard
membership, case order, solver-call ownership/budget, organism identities,
predicted endpoints and independent GPR/LP records. It recomputes agreement using
the registered tolerances, verifies independent source/package versions and
checks all declared counts/statuses against the terminal exit code. Input hashes
are rechecked after the audit. Missing or interrupted verification is incomplete,
never accepted because a comparison prefix or a success-looking log exists.

The audit separately compares growth classifications. Scalar agreement within
`1e-6` can straddle the registered `1e-6` growth threshold: a synthetic regression
uses `0.8e-6` versus `1.2e-6`, which passes the original scalar tolerance but must
fail endpoint agreement. This strengthens final verification without changing
the frozen verifier's original reports or production predictions.

Sixteen focused tests cover complete and partial coverage, resealed malformed
outputs, counts, identities, GPRs, primal validity, thresholds, input mutation,
solver budgets and replay history. The full benchmark gate passes **683 tests**
with thirteen optional skips. An initial test changed the prediction counts as
well as its intended endpoint fault; its fixture was corrected to reach the
endpoint check. A partial-prefix file handle was also closed explicitly. Initial
failures and successive auditor versions are retained.

## First terminal shards

Slurm task `10150_0` (job `10151`) finished all **1,320** plans and **6,600**
independent LP comparisons. All **2,640** GPR comparisons agree. There are
**6,599** scalar agreements and one production solver unknown; all 6,599
comparable growth classifications also agree. The terminal exit code is **1**,
matching the verifier's `disagreements` status. It is not recorded as a passed
numerical gate.

Task `10150_1` (job `10152`) completed production but reached its time limit during
verification. Its exit code is **143** and its preserved prefix contains **6,204**
comparison rows, all recorded as agreeing. With no completed verification
receipt, these rows remain an unverified prefix. The collection audit is
`incomplete`: the expected denominator is **39,600**, not the available 6,600.
Both terminal directories were collected intact and hashed against their frozen
package. The Slurm accounting service is unavailable; `scontrol` terminal job
states and persisted exit receipts establish termination.

Job **10158** reruns only the independent verification for shard 1, using the
unchanged complete production file. It has a three-hour allowance. The submission
script requires an interruption receipt and no completed oracle receipt; Slurm
terminal state was checked before submission. Original partial comparisons and the exit code are preserved in
the new job; copied hashes, submission identity and the reason for replay are
recorded. No partial comparison is reused as a new solve. The auditor compares
the completed replay against every complete prior-prefix row and retains any
differences, even if the new run otherwise succeeds. Tasks 2–3 are running with
three-hour allowances; tasks 4–5 remain queued at this checkpoint.

## Reproduced numerical exception

The exception is `NUTRIENT2026:Bs:AmSu + Gluc:R+C`, arm `printed:oxygen=0`,
policy `no_unsupplied_uptake`. Production HiGHS reports status 4 and retains an
`inconclusive`/`solver_error` result with no growth value. The independent SBML
oracle returns optimal growth approximately `1.97e-31`, with valid residuals.
Neither result is treated as an experimentally observed phenotype.

Fresh-process diagnostic job **10159** reproduces all five production prediction
states, including the same error. Its raw HiGHS message reports an unknown model
status and an infeasible primal status; that message does not establish LP
infeasibility. The exact sparse production matrices, objective and bounds for all
five calls were retained without changing the original calls.

Job **10160** then tests the exact failing production LP (2,211 rows, 2,959
variables, 11,398 nonzeros), rather than only the separately constructed oracle.
All five diagnostic variants are retained:

| Solver configuration | Result | Objective | Maximum balance residual |
|---|---|---:|---:|
| Original `highs` defaults | Status 4 / unknown | unavailable | unavailable |
| Explicit `highs-ds` | Status 4 / unknown | unavailable | unavailable |
| `highs-ipm` | Optimal | 0 | `1.35e-9` |
| `highs`, presolve disabled | Optimal | 0 | `6.50e-13` |
| GLPK 5.0 on the same matrix | Optimal | 0 | `6.99e-9` |

All successful primal solutions also satisfy bounds within `1e-6`. The diagnostic
environment is Python 3.12.3 and SciPy 1.15.3, matching the remote execution.
This locates a reproducible numerical failure in the default HiGHS configuration;
it does not justify relabelling the frozen unknown or choosing a solver based on
experimental outcomes. The frozen predictions remain unchanged. Diagnostic calls
are additional verification work, separate from each production plan's five-call
budget.

The [checkpoint archive](../../benchmark/v3/results/nutrient-community-audit-v1/README.md)
retains both terminal jobs, frozen inputs, audit versions and receipts, tests,
submission scripts, diagnostic matrices/primal solutions and logs. Full shard
verification and the frozen biological endpoint analysis remain pending.

## Completed independent replay

Replay job **10158** completed on 2026-10-27 at 12:49:43 Arabia Standard Time,
with exit code **0:0** after **1:31:48**. Its complete receipt verifies all
**6,600** numerical comparisons and **2,640** GPR comparisons. The independent
collection audit checks the full replay against the frozen inputs and every
preserved prefix row: all **6,204** original records are exactly equal and
numerically equivalent, with no differences or truncated trailing record.

The updated collection audit covers **13,200 of 39,600** comparisons. Of these,
**13,199** agree; all 13,199 comparable growth classifications also agree. The
original shard-0 unknown remains unchanged. Shards 2–5 are still outstanding, so
collection status remains `incomplete` and the biological endpoint analysis has
not been run. This is numerical verification, not empirical phenotype accuracy.

The [replay archive](../../benchmark/v3/results/nutrient-community-oracle-replay-v1/README.md)
retains the complete replay, original interruption history, updated audit and
bound inputs, terminal scheduler record, source and logs. The earlier checkpoint
archive remains unchanged.

## Three-shard checkpoint

Task **10150_2** (job **10156**) completed on 2026-10-27 at 13:18:15 Arabia
Standard Time, exit **0:0**, after **2:05:29**. Its complete receipt verifies
all **6,600** numerical and **2,640** GPR comparisons. The independent collection
audit now covers **19,800/39,600** comparisons: **19,799** agree, and all 19,799
comparable growth classifications agree. The original shard-0 unknown remains
unchanged; there are no replay-prefix or endpoint discrepancies.

Shards 3–5 remain outstanding at this checkpoint. The collection is still
`incomplete`; biological endpoint evaluation remains gated on complete collection.
The [three-shard checkpoint archive](../../benchmark/v3/results/nutrient-community-three-shards-v1/README.md)
retains the new terminal job, bound inputs, complete audit and scheduler evidence.
Earlier immutable checkpoint archives remain unchanged.

## Four-shard checkpoint

Task **10150_3** (job **10157**) completed on 2026-10-27 at 13:23:20 Arabia
Standard Time, exit **0:0**, after **2:10:14**. All **6,600** numerical and **2,640**
GPR comparisons agree. The independent collection audit now covers
**26,400/39,600** numerical comparisons: **26,399** agree, including every
comparable growth classification, and the original shard-0 unknown remains.
There are no new discrepancies. Shards 4–5 are still running, so the audit
correctly exits 2 with status `incomplete`; biological endpoint analysis is gated.
The current receipt is
`~/.cache/gspa-v3-data/nutrient-community-v1/audit-four-shards-v1.json`.
Earlier immutable checkpoint archives retain their original scope.

The [four-shard archive](../../benchmark/v3/results/nutrient-community-four-shards-v1/README.md)
preserves the completed shards, independent audit and terminal scheduler evidence.
It remains an incomplete collection checkpoint.

## Complete collection

The remaining two shards completed with exit zero. Both local and remote full
audits cover all **39,600** LP comparisons: **39,599 agreements**, the original
solver unknown, and no further numerical or comparable-endpoint discrepancies.
All **15,840** GPR comparisons agree; the replay retains all 6,204 exact prefix
matches. Audit status is `disagreements`, preserving the original inconclusive
result; it is not relabelled as a successful production solve.

Job 10167 then completed the frozen endpoint analysis. Its full output and
independent recomputation are described in the
[complete result record](v3-nutrient-community-results.md). This completes
numerical collection and endpoint scoring, not biological cross-feeding or
whole-item acceptance.
