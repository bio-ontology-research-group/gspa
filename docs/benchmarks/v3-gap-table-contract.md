# Reaction-level gap preparation

`benchmark/parse_gapseq_gaps.py` previously emitted a gap for every row whose
status differed from `good_blast`. A reaction with both a strong and a weak
alignment therefore became a gap. Missing reference sequences, spontaneous
reactions and missing complex subunits also became whole-enzyme gaps. The
`gapseq_guessed` flag could depend on which duplicate row appeared first.

The converter now groups all rows by pathway and reaction before deciding
eligibility. It emits a candidate gap only in a predicted pathway, with exclusively
`bad_blast`/`no_blast` evidence, and without spontaneous or complex annotations.
Any strong hit or unknown reference coverage suppresses that gap. Any weak hit
sets `gapseq_guessed`, independent of row order. All EC and mapped GO alternatives
remain explicit; absent identifiers stay null. These are candidate gaps supported
by an unsuccessful sequence search, not proof that an enzyme is absent.

Named columns, row widths, statuses, identifiers and booleans are checked.
Malformed input fails before opening the output, preserving an existing output
file. A valid header with no rows remains a valid empty table. EC prefixes are
normalized before mapping; invalid eligible EC strings fail.

The regression suite invokes the real converter CLI on synthetic named-column
tables. Eleven of the first twelve checks fail before repair; all seventeen final
checks pass. The broader benchmark suite runs 569 tests without failures, with
12 optional skips (nine Nextflow checks separately passed on the unchanged JVM
engine in the reaction-local repair). Evidence is in
`benchmark/v3/results/gap-table-contract-v1/`.

The status interpretation was checked against local gapsmith source in
`crates/gapsmith-find/src/types.rs`, `runner.rs` and `output.rs`. This does not
establish that an installed binary matches that checkout. The subsequent [complete observed-table audit](v3-gapsmith-observed.md) now
verifies this converter against a successful actual search, job 10111. It yields
50 candidate records over 46 eligible reactions, versus 234 old records; strong
hits, complexes, spontaneous reactions and unknown reference coverage account
for the excluded evidence. Both the initial staging failure and successful
repeat are preserved. This remains an engineering audit, not biological truth.

Complexes require role/subunit-aware inference and are deliberately outside this
whole-enzyme converter's current scope. Its historical `MetaCyc:` pathway prefix
is retained; mixed pathway databases still need explicit namespace handling.
The JVM adapter now passes the separate observed-output repair; full
reaction-context ablations and biological validation remain open.
