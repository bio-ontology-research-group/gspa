# Gap-table contract regression

`before-parser.py` preserves the original converter. The original twelve tests
produced eleven failures; `before-tests.log` retains them. The final seventeen
CLI-based tests pass after reaction-level aggregation and strict table validation.
`benchmark-tests.log` records 569 tests, no failures, 12 optional skips.

All inputs to these tests are explicitly synthetic tables. No biological accuracy
is claimed here. The subsequent genuine external-output validation is recorded
in `../gapsmith-observed-v2/gap-audit/`. See
[the repair record](../../../../docs/benchmarks/v3-gap-table-contract.md).
`source-files.json` binds the repaired parser and final tests to their repository
paths; `sha256.json` covers all archive files except itself.
