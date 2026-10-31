# Reaction-local missing-evidence regression

Three actual-suggester regressions fail on the old code and pass on the repair.
`before-tests.xml` retains the failures. `cli/verification.json` verifies ten
production CLI runs, five interventions per frozen engine, on explicitly
synthetic inputs. Functional TSVs remain byte-identical; suggestion features
no longer treat absent GO evidence as probability 0.5.

`engine.json` identifies the repaired JAR and all 1,113 JVM tests (two skips).
The nine live Nextflow checks pass. `source/` preserves the implementation,
regression and replay scripts; `source-files.json` records their repository paths
and hashes. `sha256.json` covers every archive file except itself. Absolute paths
inside original receipts refer to the original run locations and are preserved.

See [the repair record](../../../../docs/benchmarks/v3-rlc-missing-evidence.md)
for numerical results, reproduction and limitations. No biological accuracy or
whole-item acceptance is claimed.
