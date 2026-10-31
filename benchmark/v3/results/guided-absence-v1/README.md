# Guided absence mechanism: verified scope

Automatic construction connects run gaps and retained raw candidates to separate
absence comparisons with shared budgets. These are synthetic mechanism checks,
not biological absence classification or gene-function validation.

Final CLI evidence is in `verification/cli-replay/`: six packaged CLI runs,
including a baseline, and 100 independent COBRApy/GLPK phenotype comparisons.
All agree, while annotation and learned-feature exports remain byte-identical.
The full JVM gate has 1,265 tests, no failures and six optional skips; 50 focused
tests and 89 Python mechanism tests pass. Exact commands and package versions
are retained. Four optional JVM skips require the installed-HMMER environment
variable; two are existing optional core tests.

`focused-initial` retains the mapping callback's set/list failure.
`focused-second` retains a test assertion using nonexistent annotation fields.
`focused-third.log` and `focused-final-xml/` are the passing focused results.
`full-jvm-xml/` and `full-jvm.log` are the passing full-suite results.

Source snapshots, inputs, outputs, child logs and construction provenance are
retained. The immutable engine location and SHA-256 are in
`verification/engine.json`; the JAR is cached locally, not duplicated here.
The replay runner and receipts retain original workspace paths: remap them when
restoring on another host. The source report documents the input contract and
reproduction command. Large files use deterministic gzip; `storage-index.json`
records their uncompressed hashes and sizes. Every stored payload was expanded
and checked against both its recorded hash and the original source.

Run `sha256sum -c SHA256SUMS` to check stored files. No whole audit item is accepted.
