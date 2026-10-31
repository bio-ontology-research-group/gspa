# Local reaction context without peers: verification archive

The production bridge previously discarded local observations when no independent
peer genome existed. The repair retains local/intrinsic features in a `localOnly`
payload with peer statistics null and recurrence availability zero.

`cli/verification.json` records eight real production CLI runs on 29,079 unchanged
target homology claims and 50 verified pathway-gap hypotheses. Singleton local
features now cover 1,182 candidates, including 20 with nonzero anchors. Missing
coordinates withhold local windows. With two independent peers, every old/new
feature value is identical. Sequence features, candidate identities and baseline
prediction TSVs match in every arm. Both current learned models abstain on the
1,182 local-context candidates in the singleton arm.

`failed-cli/` preserves the initial invocation rejected for supplying metadata
for genomes outside the chosen cohort. `verification/` contains the four
failing-before tests, focused/final checks, intermediate export-reader and
fixture failures, the engine receipt, full test logs and prior source snapshots.
The final JVM gate has 1,172 tests, two skips, zero failures/errors. The benchmark
gate has 572 tests with twelve optional skips; all nine live Nextflow cases pass
separately. `test-results/` contains complete final JVM XML results.

Files larger than 512 KiB are stored with gzip compression. `storage-index.json`
maps each original archive-relative path to its stored file, uncompressed size,
SHA-256 and original execution path. Decompress before comparing against CLI
receipt hashes. `sha256.json` covers all stored archive bytes except itself.
`source-files.json` separately binds the final code and repair document.

`inputs/` snapshots the bound data, models and reference files; engines and the
Java executable remain identified by their hashes and original paths in the
receipt. The genuine external-search provenance is preserved in
`../gapsmith-observed-v2/`; these runs use its derived gap hypotheses, not truth
labels. The two reference genomes are used for software compatibility checks,
not a claim that this is an optimal biological cohort.

See [the repair record](../../../../docs/benchmarks/v3-local-reaction-context.md)
for semantics, limits and remaining work. No whole audit item is accepted.
