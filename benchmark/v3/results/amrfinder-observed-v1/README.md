# AMRFinder observed adapter repair

Actual AMRFinderPlus 4.2.7 execution, database 2026-08-07.1, and official version-3/4
fixtures expose the old fixed-column parser and missing pipeline registration.
The repaired engine preserves all 14 observed calls and raw provenance. Seven live
CLI checks verify success, zero calls, invalid input/database, missing executable,
and annotation/direct integration parity. Twenty-one focused cases and 1,327 JVM
tests pass (six optional skips). Biological accuracy and the full 12-item plan
remain unaccepted. See [the report](../../../../docs/benchmarks/v3-amrfinder-observed.md).

`index.json` maps source paths to content-addressed gzip objects. Both stored and
expanded SHA-256 hashes are verified. Before sources, failed trials, original
downloads, all bound execution inputs, current sources and final JVM XML are kept.
Route v2 predates the final stderr-tail reporting repair; route v3 is authoritative.
`SHA256SUMS` covers every archive file except itself.
