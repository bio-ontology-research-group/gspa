# Reaction-signature availability repair, 2026-10-27

Verified engineering repair; no whole audit item accepted. See the
[protocol, results and limitations](../../../../docs/benchmarks/v3-signature-availability.md).

Untestable signatures no longer create negative recurrence evidence. Family
presence, copy observations and layout remain available; unsupported local
context is withheld from learned inputs. Eight added JVM regressions cover the
repair, including three initial cohort failures and an intermediate bridge
failure. Genuine negative neighborhoods and family absence remain distinct.

- `cli/initial/`: retained first CLI diagnostic failure, caused by a mistaken
  test expectation for a second one-neighbor reaction under a two-anchor minimum.
- `cli/final/`: ten successful before/after production calls, preserving all
  candidates and functional TSV bytes across five controlled scenarios.
- `migration/before/` and `migration/after/`: complete matched eighteen-genome
  exports, exact commands, input/output bindings, supervision and model refit.
  All 484,202 candidate identities and non-cohort values survive; 4,780 candidates'
  context features change. Every functional prediction TSV stays byte-identical.
- `models/`: the two refitted models, also preserved inside the after-migration
  outputs. All 4,262 original supervision records and splits remain unchanged.
  Null layers are untrained; the target is reference support, not experimental
  enzyme correctness. Current benchmark runners default to these models.
- `controls/`: 118 independently verified whole-genome outputs. For each model,
  all 346,269 numeric probabilities remain exactly invariant to duplicate strains,
  with 137,933 unsupported pairs and no lost identities or support transitions.
  Frozen policies, test reports and all-genome comparisons preserve the failed
  calibration transfer and the lack of evidence for a neighborhood benefit.
- `source/`: production/test changes and frozen scripts bound by the export.
  `source/current-defaults/` separately records the later default-model selection
  and Nextflow fixture updates; historical receipts retain the source actually used.
- `verification/`: regression failures, full JVM XML/logs, live HMMER checks,
  benchmark/Nextflow gates, independent real-cohort comparison and archive scripts.

The JVM gate passes 1,229 tests with six optional skips. All four skipped live
HMMER/AntiFam cases pass separately, leaving two unrelated optional checks
unexecuted. The benchmark gate passes 600 tests, three optional skips, including
all nine live Nextflow checks; the final gate uses the refitted model.

The logistic whole-genome EC policy selects 0.45 with 25/26 calibration calls
supported by reference annotations; test calls are 6 supported/6 unverified.
Shuffle yields 7/5, net +1 with a paired genome-bootstrap interval of −3 to +3.
Logistic GO and both tree policies abstain. These findings do not certify
biological prediction accuracy, honest null inference or context benefit.

Production engines are retained in the local immutable engine cache:

- Before: `1f0a7da89894896b05626d386879be494fd1e8d1ecbacba3b048e600086befd6`.
- After: `6757088967769a6d9204f50b34a7e96ebf42a177e4ae50588db27a562cdc3dd7`.

`storage-index.json` maps every payload to its original path, uncompressed size
and SHA-256. Payloads above 512 KiB use deterministic gzip. All archived payloads
are decompressed and independently rehashed; `SHA256SUMS` covers every stored
file except itself. The large original inputs remain in their pinned caches or
earlier archives, with their hashes in the execution/migration receipts.
