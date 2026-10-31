# Genuine gapsmith observation and production adapter verification

Job 10111 completed with exit 0. `observation/` contains the full observed tables,
receipt, 41,274-file reference inventory and execution logs/resources. The exact
training input is `input.faa`; no held-out function labels were supplied.

`final-cli-replay/verification.json` checks seven real CLI invocations using those
frozen outputs. This is a replay of an actual external search, not seven new
searches. It verifies 3,469 candidates, all 13,147 raw rows, raw-claim preservation,
annotation/direct TSV parity, both models' unsupported-source abstention and
empty/failed states. `gap-audit/verification.json` verifies 50 new candidate-gap
records across 46 eligible reactions, against 234 old records.

`before-parser.json` and `after-parser.json` retain the full parser comparison:
zero versus 3,469 candidates, for the same observed table. `parser-comparison.json`
identifies both engines. `engine.json` identifies the final JAR and the complete
1,147-test JVM gate (two skips). Nine live Nextflow cases pass. The separate
fixture-order gate checks the synthetic fixture after placing its nonempty status
column last; this preserves its named-column semantics and avoids stripping TSV
trailing fields. The unmodified observed fixture retains all original fields.

All intermediate compilation, inherited-dispatch and complete-table failures
are retained, including the failed first route replay. The initial external
staging failure is preserved separately under `../gapsmith-observed-v1/`.
`source/` and `source-files.json` bind the final implementation, checks and tests.
`inspected-gapseq_find.sh` is the installed shell source used to verify gapseq
flags; this archive does not certify a full gapseq execution.

`sha256.json` covers every archive file except itself. Absolute paths inside
receipts refer to original execution locations and are preserved. See
[the repair record](../../../../docs/benchmarks/v3-gapsmith-observed.md) for scope,
reproduction and remaining requirements. No biological accuracy or whole-item
acceptance is claimed.

`cli-replay/` retains the earlier successful replay before stale-file cleanup;
`final-cli-replay/` uses the final engine. `test-results/` preserves the complete
final JVM XML results.
