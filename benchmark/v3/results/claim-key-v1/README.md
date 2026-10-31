# Protein identifiers containing pipes

The old typed key parser splits at the first two pipes and mistakes part of the
protein ID for the annotation type. Targeted essentiality rescue consequently
crashes on a UniProt-shaped identifier. The repaired parser shares integration's
last-two-delimiters rule. Three new round-trip regressions fail before repair.

`initial-reproduction/` retains the first actual CLI crash and ordinary-ID control.
`cli/` contains a four-run replay on the first repaired build; `cli-indexed/`
repeats it on the subsequent essentiality-index build. Old pipe case fails,
repaired pipe case succeeds, ordinary-ID scores remain unchanged, and repaired
predictions agree except for the deliberately different protein identifier.
Input/output hashes were rechecked after both replays.

The original exploratory essentiality replay crash is retained in `tests/`.
`code/` contains the parser, regression cases and replay scripts. Inputs are
synthetic; this does not establish biological accuracy. See
`docs/benchmarks/v3-essentiality-index.md` and `tasks/v3-acceptance.md`.
