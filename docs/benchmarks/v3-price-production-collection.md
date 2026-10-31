# Price production collection: registration and streaming verification

This records the original v1 execution checkpoint. The active full-scope run now
uses [local execution and checked publication](v3-price-local-execution.md), with
all 132 cells registered under the tested compact-audit JAR. Original v1 results
remain preserved and are not silently combined with that execution snapshot.

All 132 single-source parts (six predictors × 22 assemblies) are registered for
production integration without priors and with default priors. The frozen inputs
contain 1,526,919,081 claims. This is an ongoing inference collection, not a complete
accuracy benchmark. Ensemble, cohort, learned-model and biological comparisons
remain required. All twelve acceptance items remain open.

The existing complete-matrix export was independently audited in job 10109; see
[the export record](v3-price-streaming-export.md). Registration job 10227 checks
that audited plan, every publication receipt and genome layout, freezes the tested
job-10221 JAR and reference files, and hashes the runner sources, parser implementation,
bundled native library and runtime executables. No evaluation labels are opened.
The collection is retained at `/data/hohndor/gspa-v3/price-production-v1/`.

The first wave covers all 66 TopEC, MAPred and EnzBERT parts, with 13,682,870
claims: approximately 0.896% of the full claim population. Pilot job 10228 completed
index 0 successfully, preserving all 3,481 candidates in both arms (111.05 seconds
including verification). Array 10229 is now running indices 1–65 after pilot success. Both arms
must finish and preserve the exact candidate identities before a cell is accepted.
The two arm names, `stage0` and `stage1`, follow the existing capacity runner:
`stage0` still uses the production evidence combiner without priors. It is **not**
the native predictor score baseline.

Job 10235 completes both registered arms for DPZ index 66 (9,790,188 claims)
on node002 with a 64 GiB heap and 80 GiB allocation. Every candidate is retained,
and its default-prior TSV is byte-identical to the retained indexed-JAR result.
This is a complete single-source cell, not an accuracy comparison. Its receipts
also provide the baseline for the [compact audit replay](v3-compact-audit.md).

Job 10230 started after pilot success and separately profiles the largest CLEAN part:
45,290,880 claims for GCA_004402975.1. It uses a 168 GiB Java heap inside a 192 GiB
allocation on node005, respecting the node's 48 GiB service reservation. It retains
bounded Java Flight Recorder output, periodic memory readings, logs and failure
records. It runs only the arm without priors, so even a successful diagnostic
cannot count as a completed two-arm collection cell. Its JVM now exits successfully
after 35:39.63 with 180,227,072 KiB peak RSS, after an observed page-writeback pause
while closing the TSV. The full diagnostic subsequently verifies all 45,290,880
candidate identities in 5,087.89 seconds total. This remains a single-arm capacity
result, not a complete two-arm cell or an accuracy result.

## Verification repair and measurements

Previously, `verified_batch.output_hashes` used `json.loads` on the complete stage
ledger. One retained DPZ ledger is 2,683,657,262 bytes, mostly individual candidate
changes; loading it into Python could consume much more memory than the file.
The verifier now parses the entire JSON stream while retaining only the schema,
genome identity and stage metadata used for acceptance. It still hashes the full
ledger and prediction file, rejects incomplete/failed/duplicate stages, and checks
learned-model identity when applicable. The small-file fallback applies the same
metadata rules and rejects nonfinite JSON. Without the optional parser, ledgers
over 8 MiB fail explicitly. Learned summary reports have an 8 MiB guard too.

The dependency is pinned in `benchmark/v3/requirements-streaming.txt`.
The [ijson event API](https://github.com/ICRAR/ijson#lower-level-interfaces) avoids
constructing a whole stage's change list. Production replay used ijson
3.4.0.post0 with `yajl2_c`; both the extension and bundled libyajl are hash-bound.
Memory remains dependent on individual JSON token size and retained metadata;
this is a measurement for the actual ledger shape, not a bound for arbitrary JSON.

Final replay 10226 passes 24 focused tests and rechecks the retained 9,790,188-row
DPZ prediction and its ledger in 53.75 seconds, with 23,424 KiB peak process RSS
(22.88 MiB). Both SHA-256 values exactly match the original verified receipt:

* Prediction: `64701c3dd55cd6dc753139df7b93eb10e49f1efad361ea943f19f42f7aa6136f`.
* Ledger: `f491cd3ce5b450b525ba450fb0c375702587f56d0d3736d3c4f2f9d572e1b3ae`.

The broader local focused gate passes 44 tests without running models or JVMs.
Tests cover malformed/truncated skipped payloads, trailing data, duplicate metadata,
small/streaming acceptance parity, model-hash mismatches, candidate preservation,
failed/incomplete/duplicate attempts and changed output/receipt rejection. The
Python-memory fixture checks accumulation; the real remote replay measures total
process RSS, including native allocations. No old-parser peak RSS was measured,
so there is no numerical before/after speedup or memory-reduction claim.

`run_price_production_collection.py finalize` requires exactly one completed
attempt for every registered cell, both requested arms, matching populations,
valid chained receipts, unchanged predictions/ledgers and unchanged inference
inputs. It refuses partial collections and ambiguous duplicate completions.
Started markers alone do not establish liveness or success; scheduler termination
can leave an unfinished attempt without a Python failure record.

The portable checkpoint is in
[`benchmark/v3/results/price-production-checkpoint-v1/`](../../benchmark/v3/results/price-production-checkpoint-v1/).
It retains registration, code, tests, review adjudication and measured replay
receipts. The large raw inputs and ongoing production outputs remain at the remote
locations bound by those receipts. This checkpoint does not establish improved EC
prediction, learned calibration or completion of the full benchmark.

The separate caller `v3_runner.stage_artifacts` was subsequently connected to the
shared parser and verified in [job 10244](v3-generic-ledger.md). This running
collection retains its original frozen verifier; the later extraction does not
change its registered source or outputs.
