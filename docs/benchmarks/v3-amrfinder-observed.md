# AMRFinder output and execution contract

The previous adapter parsed AMRFinderPlus 4.2.7's real header as an annotation
for a fake protein named `Protein id`. It assigned zero identity to all 14 real
calls and shifted their classifications. Its synthetic test used an invented
column order that concealed the defect. Separately, requesting `amrfinder` in
`predictors.enable` registered no predictor and exported no claims.

## Repair

The parser now resolves named columns from the published version-3 and version-4
formats, including optional coordinate columns, and preserves the original row
fields. Missing/ambiguous headers, truncated rows, invalid percentages, and
missing required values raise a typed output-contract failure. A complete header
without data means a successful search with no reported calls.

`predictors.amrfinder` supplies `database`, `executablePath`, `organism`, and
`threads`; `predictors.enable: [amrfinder]` requests execution. `--threads`
overrides the configured count. The command uses protein input and `--plus`.
AMR, stress, and virulence classifications retain their distinct values and
metadata, while using the existing AMR annotation channel.

Reported identity divided by 100 remains an explicitly uncalibrated similarity
score. It is not a probability of resistance. For HMM calls lacking BLAST
statistics, `identity` and `coverage` remain null and `score_available` is false.
The primitive score field uses a zero storage placeholder; claim extraction,
JSONL round trips, evidence combination, refinement, learned feature extraction,
and cohort baseline construction do not treat that placeholder as functional
evidence. Such observations remain in raw claims for audit. Other scored calls
still enter the existing legacy calibration transform; no new calibration or
biological accuracy is claimed.

NCBI documents the [report columns and missing alignment fields](https://github.com/ncbi/amr/wiki/Running-AMRFinderPlus#fields)
and explains that AMRFinderPlus [does not predict resistance phenotypes](https://github.com/ncbi/amr/wiki/Interpreting-results).

## Evidence and checks

The [official 4.2.7 binary release](https://github.com/ncbi/amr/releases/tag/amrfinder_v4.2.7)
was downloaded and its SHA-256 verified against the release API:
`68045a8bccdbe3c5dcdf941bebe2352ed419758a9914c41f48f0bbbd6fbade56`.
The tool ran on its supplied 17-protein test FASTA against database
`2026-08-07.1`, with `--plus --organism Escherichia --threads 1`, yielding 14
calls. The independently captured output, log and input are stored under
`gspa-predictors/src/test/resources/observed-output/amrfinder/v4/`.
Official coordinate-bearing expected outputs from releases 4.2.7 and
[3.12.8](https://github.com/ncbi/amr/blob/amrfinder_v3.12.8/test_prot.expected)
provide additional parser fixtures. Version 3 was not locally executed.

All 16 initial regressions fail against the old implementation. The final
contract tests also cover process-level failure recording, version provenance,
missing-score inference exclusion, and pipeline configuration. Malformed and
missing-statistic variants are explicitly synthetic fault tests, not additional
observed tool outputs. The fabricated old fixture was removed.

`benchmark/v3/check_amrfinder_route.py` binds the full tool and database trees,
input fixtures, engine, Java executable, BLAST/HMMER executables and checker.
It compares real annotation claims with the independently generated reference,
checks the actual thread/database arguments, and compares raw provenance and
integrated output across annotation and direct integration. Missing executable,
invalid database, invalid input and genuine header-only output exercise distinct
completion statuses. The initial failed no-hit test is retained: its artificial
sequence was rejected as nucleotide-like, correctly producing FAILED. The final
no-hit case uses the three upstream proteins unreported under this organism
setting; this selection establishes a parser/run contract, not biological
negative labels.

The final JVM gate passes **1,327 tests with no failures and six optional skips**;
the AMRFinder-specific specification contains 21 passing cases. All seven final
CLI executions pass their assertions: two successful annotation routes, genuine
zero calls, invalid input, invalid database, missing executable, and direct
integration. The two annotation routes preserve all 14 reference calls and raw
provenance; annotation-integrated and direct-integrated predictions are byte-identical.
The status manifest now retains the stderr tail, so long banners cannot hide the
actual invalid-input or database diagnostic. The CLI's overall exit code remains
zero for recorded per-predictor failures; consumers must inspect the manifest.

The previous engine is
`f80982650f49738235e80e0570ce3dd416270f64818126c456cabe4c48d3c345`;
the repaired engine is
`74d9e8fc8c170adcb8d222cc0fb2ba82199b877e63a51efffca40f8dac512405`.

Sources, downloaded release, database, failed and successful runs, old parser
output, source snapshots, and test results are preserved in
`benchmark/v3/results/amrfinder-observed-v1/` with stored and decompressed hashes.

This closes the demonstrated AMRFinder adapter and registration defects. It does
not establish AMR biological accuracy, learned model eligibility, nucleotide
search support, or completion of the full 12-item acceptance plan.
