# InterProScan statistic contracts and observed execution

This work advances audit items 2 and 3. It does not complete predictor coverage,
per-method calibration or the biological evaluation required by items 1 and 4.

## Defects and changes

The previous parser interpreted every TSV statistic as an E-value. Valid HAMAP
and ProSite profile scores above one therefore produced negative annotation
scores; malformed text, negative E-values and missing values could instead become
0.9. A reported zero E-value received less support than a tiny nonzero E-value.
Raw statistics were discarded, including the originating analysis on InterPro
entry annotations. The annotation pipeline also omitted InterProScan's configured
thread count.

The [official TSV contract](https://interproscan-docs.readthedocs.io/en/v5/OutputFormats.html)
defines a method-dependent score column. The
[upstream writer at revision 26256cbb](https://github.com/ebi-pf-team/interproscan/blob/26256cbb70357846da3c8275df59397c537f1aca/core/io/src/main/java/uk/ac/ebi/interpro/scan/io/match/writer/ProteinMatchesTSVResultWriter.java)
selects E-values for HMM/BLAST/PRINTS matches, profile scores for profile scans,
and a missing-value marker for other methods. The source is downloaded and hashed
in the evidence archive. The [InterProScan FAQ](https://interproscan-docs.readthedocs.io/en/v5/FAQ.html)
also cautions against treating different methods' E-values as a common scale.

The repaired parser:

- Checks numeric statistics for finiteness and rejects negative/missing E-values.
  A corrupt row fails the batch; preceding rows are not reported as a successful
  partial annotation.
- Keeps the existing `-log10(evalue)/20` transform for recognized E-value methods,
  bounded to [0,1], with zero mapped to one. Ordinary E-values in (0,1] are unchanged.
- Uses the existing match-only 0.9 fallback for profile/uninterpreted/missing
  non-E-value statistics. This fixes the sign/scale error but remains a **legacy
  heuristic**, not independent calibration or a measured probability. Unknown
  methods are not silently interpreted as E-values.
- Preserves `rawStatistic`, `statisticValue`, `statisticType`, `scoreModel`,
  `scoreCalibrated: false`, originating `analysis` and `via` signature on every
  emitted annotation. Raw claim export retains this metadata.
- Honors the domain thread setting, including the CLI `--threads` override.

`predictors.domains.disablePrecalc: true` now forwards the adapter's existing
`-dp` option, allowing a benchmark to use only installed database files. It defaults
to false to preserve ordinary lookup behavior. The observed replay sets it true;
the remote lookup service cannot silently supply a different release's matches.
InterProScan 5.72's `--version` probe succeeds; no availability repair was needed.

## Independent tool execution and regression proof

Slurm job 10095 completed on node003, exit 0, using installed InterProScan
5.72-103.0. It ran Pfam 37.1, HAMAP 2023_05, ProSiteProfiles 2023_05 and
ProSitePatterns 2023_05 on two public protein sequences (P18843 and P0ADB7) and an
explicitly synthetic 80-alanine control. This is a tool-format fixture, not an
accuracy or sensitivity benchmark. Source FASTA identity and selected accessions
are recorded. The real output has four signature matches and no control match.

The generator hashes 3,951 launcher, Java, main-JAR, library, selected database,
binary, property, generator and input files before and after execution. This is
the documented identity scope, not a claim to hash every installed dependency.
Commands, return codes, versions, stdout/stderr, input, TSV and identities are
under `gspa-predictors/src/test/resources/observed-output/interproscan/`.

The old production JAR parses that same output into nine annotations, three with
negative scores derived from the HAMAP profile score 47.316956. The repaired JAR
preserves all nine identities, fixes those three scores, and preserves all six
E-value-derived scores exactly. The ProSite profile row lacks an InterPro/GO map
and this adapter does not expose a separate ProSite signature type; no extra
function is inferred from it or the pathway column. EC remains absent from the
advertised and observed output types.

Seventeen failures were reproduced in the initial focused run. The final 19-case
statistic/observed-output suite and the pipeline thread/local-lookup regression
pass. Full JVM validation passes 1,046 tests, zero failures/errors and two skips.
The CLI Python test defaults were corrected from a hardcoded RAM-temporary path
to `~/tmp/gspa-v3-testenv`; the six initially skipped checks then ran successfully.
The final benchmark suite passes 528 tests with three optional skips, including
all eight live Nextflow checks on the exact repaired JAR.

JAR before: `006ecc5339125b5845c1a27cceb175077a4f24975dc02c9f2285f414455664cc`.
JAR after: `c0def1656fc17170869b1d8f25b80e1bf433ecbc583794f63c27aff23b1bc821`.

## Reproduction and limits

Run the capture on a compute allocation with `TMPDIR` pointing to allocated scratch:

```bash
python3 benchmark/v3/generate_interpro_fixture.py \
  --tool /installed/interproscan-5.72-103.0/interproscan.sh \
  --fasta input.faa --out observed-new-run
```

The output directory must be new. `check_interpro_route.py` rehashes the captured
dependencies, executes actual annotation twice (with and without integration),
checks the raw claims and invokes direct integration. It also verifies that the
actual command has one CPU and disables remote lookup. Slurm job 10096 completed with exit 0. Both actual annotation arms return SUCCESS
and nine identical raw claims. Annotation-plus-integration and direct integration
produce byte-identical prediction TSVs. All captured dependencies are rehashed
afterward; all 24 output hashes are independently checked after archive transfer.
This certifies this route/input combination, not all predictors or biological accuracy.

`replay_interpro_corpus.py` additionally compares old/new parsers on the seven
existing full-proteome TSVs. All 379,740 emitted annotation identities are preserved.
The old parser emits 33,758 negative scores. All repaired scores lie in [0,1];
40,273 scores change in total (including reported-zero E-values). All ordinary
E-value-derived scores in (0,1] remain exactly unchanged. The new provenance
classifies 333,215 annotations as E-value-derived, 32,285 as profile-score-derived
and 14,240 as lacking a reported numeric statistic. These are annotation counts,
not independent proteins, matches or correct predictions. The full before/after
JSONL files remain in the content-bound local cache, with summary/hashes archived.
Those historical TSVs have no independently checked
execution receipts here; parser replay does not upgrade their provenance or turn
their annotation content into experimental truth.

Evidence is archived under `benchmark/v3/results/interpro-observed-v1/`. Other
adapters, full upstream coverage loss, InterProScan 6 compatibility and independently
learned method-specific score calibration remain unverified. Legacy integration
can still transform the heuristic annotation scores; the metadata flag does not
make the provisional calibration table independently calibrated.
