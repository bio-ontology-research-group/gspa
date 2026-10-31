# Automatic pathway-gap handoff: verification archive

`after-cli/verification.json` verifies ten actual production CLI invocations
replaying the frozen observed gapsmith output. Both source tables, all 3,469 raw
claims, and exactly 50 independently audited gap hypotheses survive shared
integration. Annotation, observation-file integration and explicit-gap replay
have identical feature vectors and predictions. Zero, missing and malformed
inputs remain distinct; both learned models abstain on the unsupported source.

`before-cli/` and its log retain the old engine's failing replay: it loses the
pathway table. `verification/comparison.json` binds before/after feature and raw
claim files. Only gap availability/listing/guessed-listing columns change; raw
claims and baseline TSVs are byte-identical. The corresponding comparison script
is preserved alongside the result.

`verification/` also contains the initial failing pipeline test, the reproduced
input-cleanup failure, the reproduced stale-receipt failure, intermediate test
failures and the complete final logs. Final gates: 1,211 JVM tests (two optional
skips), 573 benchmark tests (twelve optional discovery skips), and nine live
Nextflow tests passing separately. `test-results/` preserves all final JVM XML.
`verification/engine.json` identifies the tested content-addressed JAR.

Files larger than 512 KiB are gzip-compressed. `storage-index.json` maps original
archive-relative paths to stored files, uncompressed sizes/hashes and original
execution paths. Decompress before checking against CLI receipt hashes.
`sha256.json` covers every stored archive file except itself.
`source-files.json` binds final implementation, tests, checkers and documentation.
`inputs/` contains the bound observed tables, gap oracle, mapping, proteome and
models; engines and Java remain identified by their hashes and execution paths.

The actual external search and reference inventory are preserved in
`../gapsmith-observed-v2/`. This archive does not certify another search or
biological correctness of the gap hypotheses. The new source-qualified pathway
IDs deliberately replace the historical converter's blanket MetaCyc prefix.
See [the repair record](../../../../docs/benchmarks/v3-pathway-gap-bridge.md) for
behavior, exact scope and outstanding acceptance requirements.
