# Shared ledger validation for dense production benchmarks

The generic benchmark runner now uses the same streaming stage-ledger reader as
`verified_batch`. Previously, `v3_runner.stage_artifacts` loaded the complete JSON
ledger into memory. It also accepted nonfinite JSON constants in ignored trailing
fields. Two regression tests reproduce these failures before the change.

The existing parser was extracted into `benchmark/v3/stage_ledger.py` without
changing the verified-batch acceptance logic or its artifact hashes. The old
helper names remain available through `verified_batch`. The generic runner now
validates the entire JSON stream, retains only the metadata used for acceptance,
and continues to hash all emitted sidecars. Failed or malformed evidence remains
an execution failure even when the engine exits successfully.

Both runner identities bind the shared parser source. The generic runner also
records the parser version/backend and hashes its package and bundled native
library. Tests show that changing the parser library, backend or shared source
invalidates reuse; an old receipt without the new source binding is rejected.
This intentionally changes cache identity. Historical runs must be audited through
read-only verification, rather than passed to a runner that could rerun them.
The [runtime instructions](v3-protocol.md#ledger-verification-runtime) describe the
required module and pinned dependency for files larger than 8 MiB.

## Evidence

Local validation and final cluster job **10244** each exercise 106 cases:
105 pass and the optional real-JVM smoke test is skipped. No application models or
JVMs were run on the laptop. Tests cover the existing runner's label isolation,
input/engine changes, stale outputs, stage failures, sidecar integrity and metrics,
as well as the shared parser and full-collection completion checks. These are
execution-contract tests, not biological validation.

Job 10244 then replays both new validation entry points on the retained
2,683,657,262-byte DPZ ledger and 9,790,188-row prediction. The verified-batch check
takes 69.59 seconds; the generic sidecar check takes 58.87 seconds. Peak RSS for the
whole replay process is 24,400 KiB (23.83 MiB), including native allocations.
The timings cover different artifact sets and are not a speed comparison.
Both entry points reproduce the historical artifact hashes exactly:

* Prediction: `64701c3dd55cd6dc753139df7b93eb10e49f1efad361ea943f19f42f7aa6136f`.
* Ledger: `f491cd3ce5b450b525ba450fb0c375702587f56d0d3736d3c4f2f9d572e1b3ae`.

The original receipt and output files are unchanged. The replay reads no labels,
performs no inference or model fitting, and writes a separate audit. The running
Price collection keeps its independently frozen earlier verifier; its registration
and source files were not overwritten.

The first cluster attempt, **10241**, failed before replay because the source
bundle omitted the existing `benchmark/neural/make_temporal_split.py` dependency
of a test helper. Its log and source bundle are retained. A fresh remote directory
and complete bundle were used for 10244; the failed attempt is not counted as a
successful validation.

## Scope and remaining work

This closes the two identified full-ledger reads in the benchmark entry points.
Memory still depends on individual JSON token size and retained metadata; this
measurement applies to the real ledger shape. The parser hashes cover the bundled
native library used here, not arbitrary system-linked libraries or an entire
Python environment. Parsing and subsequent artifact hashing are separate reads,
so verification assumes finished, stable output files.

The portable evidence is in
[`benchmark/v3/results/generic-ledger-v1/`](../../benchmark/v3/results/generic-ledger-v1/).
It contains both source bundles, the failed attempt, final tests, source-bound
replay receipt and Sonnet review disposition. Large raw outputs remain at the
remote paths in the receipts. This work does not establish improved EC prediction,
complete the [Price collection](v3-price-production-collection.md), or close any of
the twelve whole acceptance items.
