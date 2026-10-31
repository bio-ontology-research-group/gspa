# Raw predictor role handoff: reproducible repair evidence

The old engine ignores raw-claim `candidate_role` during learned candidate
construction. `before-tests.xml` records twelve failing cases out of thirteen;
`before-cli/` and `before-cli.log` retain the full real-output replay that fails
the independent role oracle. The final fourteen role cases pass.

`after-cli/verification.json` verifies seven production CLI invocations on the
same frozen observed gapsmith output. Both learned models retain 829 subunit,
2,207 catalytic and 433 unspecified roles. All remain unsupported by the current
models. `comparison.json` binds the previous successful replay and the new one:
829 role identities/feature digests change; 2,640 digests do not change for each
model. Raw claims and functional prediction TSVs remain byte-identical.

`engine.json` identifies the final content-addressed JAR and test totals.
`test-results/` contains the full 1,161-test JVM XML gate, with two optional skips.
Nine live Nextflow tests pass. `nextflow-old-launcher-failure.log` preserves an
initial invocation using an older PATH launcher that rejected the installed Java;
the successful run uses the established cached Nextflow 26.04.6 launcher.

`source/` and `source-files.json` bind the implementation, regression tests, oracle
and repair record. `before-ContextFeatureExtractor.groovy` preserves the prior
implementation. `sha256.json` covers every archive file except itself. Absolute
paths in receipts are the original execution paths and have been preserved.

The external search inputs, full reference inventory, executable hashes and
Slurm receipt are in `../gapsmith-observed-v2/`; these are frozen-output CLI
replays, not fresh external searches. No independent role truth, role calibration
or whole-item acceptance is claimed. See
[the repair record](../../../../docs/benchmarks/v3-claim-role-bridge.md).
