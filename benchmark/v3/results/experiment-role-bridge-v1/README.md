# Experiment role handoff: verification evidence

Raw candidate roles now survive masked assignment construction. Complex-subunit
observations require supported partner clauses. Noncatalytic-only candidates
cannot become enzyme assignments, and malformed role metadata fails validation.
Explicit and guided absence comparisons also reject a standalone activity from
mixed catalytic/subunit observations. No scientific accuracy claim is made.

Twelve failing-before regression cases are retained in `verification/before-xml/`.
All 105 focused cases pass in `after-xml/`. The full JVM suite passes 1,278 tests
with six optional skips, recorded in `full-jvm-xml/` and `full-jvm.log`.

Final CLI evidence is `verification/cli-replay-final/`: nine matched old/new
cases, 18 invocations and 90 independent GLPK comparisons. Five invalid masked
assignments now abstain; one mixed-role absence request fails validation without
publishing a result. Eight successfully completed pairs preserve annotation
bytes. Valid default/supported cases retain identical phenotype cells; masked
public problems are unchanged. Proposed hypotheses retain role provenance.

`verification/cli-replay/` is an incomplete initial harness run: its unknown
assembly fixture omitted the required nullable completeness_fraction field.
The initial harness and failure are retained; use the final replay for results.

Sources, exact commands, package versions, inputs, logs, test XML and old/new CLI
outputs are included. Engine hashes and immutable local cache paths are in the
receipts; JARs are not duplicated. Replays preserve original workspace paths;
remap paths on another host. The source report documents the scope and runner.
Large files use deterministic gzip. The storage index records original sizes
and hashes; every payload was checked after expansion against its original file.
Run `sha256sum -c SHA256SUMS` to verify stored files.

All role, partner and model fixtures are synthetic. Biological role accuracy,
independent partner evidence and calibrated recovery remain open. No whole
audit item is accepted.
