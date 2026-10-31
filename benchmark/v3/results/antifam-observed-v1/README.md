# AntiFam adapter execution regression

The adapter invokes the same installed HMMER 3.4 used by the preceding
`hmmer-observed-v1` archive, with the same explicitly synthetic sequences/model.
The synthetic profile is not a genuine AntiFam family and these results do not
establish spurious-ORF detection accuracy. The model has no GA threshold:
successful execution uses E-value gating; the default `--cut_ga` failure is
explicitly tested and must be FAILED rather than a zero-hit result.

All 12 new cases failed before the fix and pass afterwards. They cover database
availability, the supported help/version probe, malformed statistics/accession,
live matching/no-hit output, and the real missing-GA tool failure. Thread defaults
and overrides are also checked through existing CLI and pipeline cases.
The full JVM suite passes 1,003 tests, with zero failures/errors and two skips.
All eight live Nextflow cases pass against the same repaired JAR.

`cli/` contains old/new real annotation CLI replay: the old build produces
UNAVAILABLE and no annotation; the repaired build produces SUCCESS and one
PSEUDOGENE annotation on the target, no negative-control hit, and uses one CPU.
Both runs emit zero raw integration claims. AntiFam is absent from the current
source-to-evidence map; functional rescoring from this flag remains unverified.
It would be incorrect to report a functional-confidence improvement here.

`hmmer-regression/` repeats the four-run direct-domain annotation/integration
check after extracting the shared HMMER executable probe. The two integration
routes still produce byte-identical predictions. Input/output hashes were
rechecked after both replays. The observed tool artifacts remain in the preceding
HMMER archive; `code/` freezes the changed sources and both replay scripts.

Repaired JAR SHA-256:
`0812564ad19764a65234f418bf2c539eb16578fb66ff62034421fd0ad39e4caa`.
AntiFam before-fix JAR SHA-256:
`8a9dedbd9358b1e1249dba3c3a2e1006e23c71e601e4d68beeb36fd00d27e65a`.
