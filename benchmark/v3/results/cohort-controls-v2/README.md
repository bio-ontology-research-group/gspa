# Verified development-cohort controls

See [the protocol](../../../../docs/benchmarks/v3-cohort-controls.md).

Three production arms contain 59 genome outputs; only the original four test
genomes enter evaluation. The frozen policy abstains for GO and selects EC 0.44.
EC top-1 reference agreement remains 14/37 in every arm. Duplication changes
numeric learned probabilities by up to 0.0017367852562789055: score invariance
is not established. Positive-only reference support does not measure false
assignments, and this previously examined cohort is not a fresh blind test.

`check_inputs.py` and `check_outputs.py` independently rehash sources/artifacts
and check conservation/stage coverage without importing the production verifier.
`source/` preserves the runner/evaluator versions used by the receipts. `tests/`
contains the three failing before-fix cases, ten passing focused tests, eleven
passing sequence-loss regressions and the 538-test benchmark gate (three skips).
`sha256.json` inventories this archive; large predictions stay in the cache paths
recorded in the receipts. Future source changes invalidate reuse of these runs;
the snapshots and original receipts preserve what was actually executed.
