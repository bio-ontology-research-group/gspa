# Full Price execution with local output and checked publication

The complete 132-cell comparison is registered at
`/data/hohndor/gspa-v3/price-production-v3/`: six predictors, 22 assemblies and
1,526,919,081 input claims, with both no-prior and default-prior production arms.
Job 10292 verifies the first cell. Every remaining cell has a submitted resource
wave; none is dropped or selected using evaluation labels. This is an execution
checkpoint, not a completed accuracy benchmark. Ensemble, learned/context and
biological comparisons remain separate requirements.

The previous largest-CLEAN diagnostic, job 10230, now completes independent
candidate verification: all 45,290,880 input identities appear in its output.
Its JVM takes 35:39.63 with 180,227,072 KiB peak RSS; the complete diagnostic,
including input and output verification, takes 5,087.89 seconds. It executes only
the no-prior arm. The observed writeback waits and successful final result are
recorded in the [compact-audit investigation](v3-compact-audit.md); no two-arm
largest-cell or accuracy result follows from that diagnostic.

## Compatibility with the EC identity correction

The ongoing collection retains its original JAR. A separate complete-cell replay,
job 10405, uses the tested JAR from 10402 on all 3,481 TopEC candidates under both
prior settings. Its initial byte-equality assertion fails because output EC IDs
now carry `EC:`. Independent full-output checking confirms that this prefix is
the only TSV change in either arm: scores, sources and ordering are unchanged.
Both complete stage ledgers agree after explicit EC-spelling and timing
normalization, and all 34 shared receipt input hashes match per arm. The original
failed assertion remains recorded. See the [portable compatibility evidence](../../benchmark/v3/results/price-ec-compatibility-v1/).

This one-cell result does not certify the other Price cells or improved accuracy.
The final benchmark must identify its JAR and account for the new representation;
older full-collection results cannot silently become current-version evidence.

## Storage and receipt preservation

`run_price_local.py` runs the existing registered collection engine inside a
Singularity namespace. A private directory on job scratch is mounted at the exact
logical output path under `/data`. The engine therefore writes local files while
its commands, manifests and receipts retain their final published paths. A check
inside that same namespace verifies that both paths identify the same directory
before starting inference. An ignored or failed bind cannot silently run the
engine against shared storage.

After execution, the wrapper hashes every private artifact, copies the tree to
staging outside `collection/attempts`, and checks both source and copied bytes.
It adds a separate publication receipt, flushes files and directories, and renames
the complete tree into place. Staging and destination must share a filesystem,
checked before inference. A completed marker becomes visible to the collection
finalizer only with the checked tree. Existing attempts cannot be overwritten;
failures preserve private/staged data and durable node/path/return-code records.
Successful private copies are removed only after publication.

No original production receipt is edited. The finalizer reads its original
absolute paths, checks the complete output and input hashes, and additionally
binds publication identity, inventory and hashes. Java receives an explicit
`-Djava.io.tmpdir` under private scratch; setting `TMPDIR` alone did not move the
JVM's temporary files. The container receives `TMPDIR` explicitly too.

This uses the existing empty container-directory skeleton and explicit host
filesystem bindings. It is not a general sandbox/SIF launcher. Runtime executable,
source and parser hashes are retained; they do not attest every host library or
promise durability beyond the filesystem's `fsync` semantics. Inputs still come
from shared storage. Full-scale runtime and memory gains remain to be measured.

## Validation and schedule

The missing-Java-scratch regression fails before the command change. Final local
tests run 64 cases with eight optional streaming-parser skips; job 10292 runs all
64 without skips. Tests cover copy corruption, source mutation, failed or
contradictory completion, same-attempt overwrite, cross-filesystem staging,
flush failure, symlinks, an inactive bind, early engine failure and finalization
using unchanged logical receipts. The finalizer rejects changed publication
identity, inventory and hashes.

Job 10287 fails before registration because the first test bundle omits a module
dependency. Its failure is retained. Job 10289 verifies the initial publication
pilot; adding the explicit bind check is then frozen into the final execution
snapshot. Job 10292 registers all 132 cells and completes the final pilot. Both
arms preserve all 3,481 candidates, every TSV byte and all ledger content except
stage durations relative to the original retained pilot. Host-side verification
rechecks published files and confirms the private output has been removed.

| Wave | Indices | Cells | Input claims | Submitted job | Allocation / heap |
|---|---|---:|---:|---|---|
| Pilot | 0 | 1 | 3,481 | 10292, complete | 36 GiB / 24 GiB |
| Small | 1–65 | 65 | 13,679,389 | 10296, concurrency 2 | 36 GiB / 24 GiB |
| Medium | 66–100 | 35 | 579,802,657 | 10297, concurrency 2 | 112 GiB / 96 GiB |
| Largest | 131 | 1 | 45,290,880 | 10295 | 192 GiB / 168 GiB |
| Large | 101–130 | 30 | 888,142,674 | 10298, concurrency 1 | 192 GiB / 168 GiB |

The small, medium and largest waves require successful pilot 10292 and complete
compact-audit DPZ parity job 10279. The large wave additionally requires successful
completion of the largest two-arm cell. These are scheduler dependencies, not
completed results. Small cells use node004; medium cells can use node002–004;
large cells use node005 while respecting its service reservation. Unrelated jobs
are untouched. Node008's probe lacks the required Java/container tools, so this
schedule does not use it. No cluster service or driver is changed.

The finalizer still requires exactly one completed attempt for every registered
cell and both arms, plus unchanged inputs, outputs, publication and exact candidate
populations. The earlier v1 collection and v2 publication pilot remain preserved;
their receipts are not rewritten or silently merged into this snapshot.

## Evidence and recovery

The [portable checkpoint](../../benchmark/v3/results/price-local-execution-v1/)
contains 96 files, the 14-file production bundle and 133-file test-source bundle,
the failed-package log, review/adjudication, final registration, scheduler
submissions, whole old/new pilot outputs, and the completed largest-baseline
receipts. Run:

```sh
python3 benchmark/v3/results/price-local-execution-v1/verify_archive.py
```

This checks the stored source hashes, tests, full small-output parity, publication
inventory, every scheduled cell and claim count, and the largest-baseline receipt
chain. The large remote outputs are not shipped or re-read by this portable check.

On publication failure, first confirm that the original job is terminal. Its
`local-executions/<index>-<attempt>.json` records the node, private parent, logical
parent and command; the adjacent `.result.json` binds that locator and records
the process return code. Retain the recorded node's scratch until recovery. On
that node, verify the locator and frozen input/runtime hashes, then use
`run_price_local.publish` with the recorded local attempt, logical destination,
collection staging directory and execution record including that return code.
The function refuses an existing destination or inconsistent success. If hashes
changed, retain the failure and use a newly registered execution instead of
rewriting receipts. An interrupted copy can leave a staging directory; it cannot
be mistaken for a completed collection cell.

All twelve whole acceptance items remain open. No improved enzyme accuracy or
readiness for release is claimed by this storage and scheduling checkpoint.

## Largest-cell completion and final collection audit

Job 10295 completes both arms of the largest CLEAN cell, preserving all
45,290,880 candidates per arm, and publishes its checked outputs successfully.
At this checkpoint 71 of 132 cells have completed and large wave 10298 has
started. The [retained metadata checkpoint](../../benchmark/v3/results/price-finalization-checkpoint-v1/)
checks its small receipt chain; large outputs remain remote.

Finalizer 10410 is queued after successful arrays 10297 and 10298. It uses the
original registered finalizer to check all 132 cells, both arms, complete inputs,
outputs and publication before writing `collection-verification.json`. Its
read-only preflight verifies the original registration, plan and exact eight-file
Python source set. Initial wrapper 10409 was cancelled before execution because
it assumed 14 Python files; the corrected v2 snapshot is retained separately.
This is completion of the largest execution cell and scheduling of the final
audit, not a completed collection, latest-JAR comparison or accuracy finding.
