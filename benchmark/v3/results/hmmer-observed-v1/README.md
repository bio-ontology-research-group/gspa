# HMMER observed-output regression evidence

Real HMMER 3.4 executes on explicitly synthetic sequences and a synthetic HMM.
This is not a biological accuracy or real-Pfam-family benchmark.

`observed-tool/` contains generated inputs, model, tool output and provenance.
`cli/` contains the four-run old/new annotation and direct-integration replay,
with content hashes and explicit assertions recorded in `verification.json`.
`cli-before-thread-fix/` retains the first CLI reproduction, which exposed
the ignored thread budget. `code/` freezes the changed sources and replay code.
`tests/` preserves failed before-fix cases and the final full-suite results.

The old build reports HMMER unavailable and emits no claim; the repaired build
emits exactly one expected hit, omits the negative control, respects one CPU,
and produces identical predictions through annotation and direct integration.
No whole audit item is signed off by this fixture. See
`docs/benchmarks/v3-hmmer-observed.md` and `tasks/v3-acceptance.md`.
