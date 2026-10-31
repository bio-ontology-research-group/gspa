# Candidate-count scoring regression

This archive retains the failing core test, passing regressions, and old/new
production CLI outputs. Inputs are explicitly synthetic and isolate a software
invariant; they do not establish biological accuracy.

The old build switched a weak candidate's context ceiling on at 100 candidates.
The repaired build applies the same conservative rule at every input size. Actual
CLI outputs show the old 99-to-100 score jump (0.96095 to 0.52434) disappearing.
The complete 100/101-candidate predictions are byte-identical between builds.

`cli/verification.json` binds the eight genome/build runs to their content-verified
batch receipts. `engine.json` identifies both JARs; the repaired JAR is cached
immutably, although the original after-run receipt names its build output path.
`code/` records the exact changed source, tests and replay script. `tests/` includes
the full JVM and benchmark logs: 972 JVM tests (two skips), 523 benchmark tests
(three optional skips); all eight live Nextflow tests execute. No failures.

See docs/benchmarks/v3-candidate-count.md and tasks/v3-acceptance.md. This repairs
one remaining requirement within item 4; broader training, null/role data,
calibration and independent biological validation are still required.
