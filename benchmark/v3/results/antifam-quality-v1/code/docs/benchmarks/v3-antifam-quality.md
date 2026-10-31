# AntiFam quality evidence through integration

The quality-evidence loss identified in the twelve-item review is repaired in the
working tree. This is an evidence-preservation and scoring-eligibility repair,
not a validated AntiFam accuracy result or a learned functional-confidence penalty.
The broader requirements for items 2–4 remain open.

## Evidence semantics

AntiFam detects recurring spurious ORFs. Some matches instead identify a local
boundary error: its Spurious_ORF_21 model addresses incorrectly reconstructed
prfB frameshift extensions. A match therefore does not establish whole-protein
inactivity. This distinction follows the
[original AntiFam paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC3308159/)
(Eberhardt et al., 2012; DOI 10.1093/database/bas003).

The adapter keeps the legacy `PSEUDOGENE` annotation namespace and now emits
`GENE_MODEL_QUALITY` evidence. Raw JSONL retains the observation, its raw score,
model accession, threshold mode and domain statistics. `calibrated_prob` is null:
a GA threshold-pass indicator or normalized E-value is not a probability of
functional inactivity. Legacy records with misleading positive probabilities or
functional evidence types are normalized to this quality representation on read.

One legacy flag per protein/family is retained, with **all** reported domains in
`metadata.domainsJson`. This JSON string preserves regional records through the
existing scalar metadata contract. Coordinates are protein-relative, one-based,
inclusive; model, alignment and envelope bounds are checked. Invalid coordinates
fail the adapter rather than producing a misleading flag.

Quality observations do not enter Noisy-OR, the functional candidate universe or
the cohort EC-anchor snapshot. Integration preserves their original annotations
and metadata when replacing functional annotations. Existing heuristic functional
scores remain unchanged; neither whole-gene deletion nor an arbitrary numerical
penalty is imposed.

## Learned scoring and reporting

The current models have no supported quality-evidence input. Functional candidates
on flagged proteins are explicitly `unsupported`, with reason
`gene_model_quality_not_supported_by_model` and null learned probability/margin.
Their existing heuristic predictions are retained and identified as such; this
is abstention by the learned stage, not a retraction of the underlying candidate.
An incomplete locus competition is not ranked. Unflagged candidates follow the
existing model coverage rules.

The numeric feature contract remains `gspa.context-features/3`: no feature values,
columns or fitted coefficients changed. Quality observations are retained in
candidate metadata, learned reports and feature exports. Exports without these
observations retain their previous bytes. A later model that actually learns from
quality evidence will require an explicit feature contract and suitable labelled
training/calibration data. The current training-data join and JSONL/TSV trainer
refuse flagged rows rather than silently discarding the observations. Python
probability/locus APIs accept `quality_evidence` and apply the same guard as the JVM;
callers using vectors directly must pass the corresponding observation metadata.

Every integration output now has `<out>.quality.json`. It records the original
quality observations and AntiFam scan coverage separately. Actual annotation runs
supply their own outcome; standalone integration defaults to `UNKNOWN` unless an
external-context `predictorStatus.antifam` is supplied. A reused annotation
pipeline cannot assign a previous genome's scan status to a new genome.
`SUCCESS_ZERO` is a completed scan with no hits; `FAILED`, `UNAVAILABLE`, `PARTIAL`
and unknown coverage remain distinct. Even a completed no-hit scan does not prove
that every gene model is correct. CLI replay with matched raw evidence and scan
coverage yields identical quality reports between annotation and direct integration.

## Verification

The [archive](../../benchmark/v3/results/antifam-quality-v1/) preserves sources,
regression logs, actual CLI outputs, engine identities and a SHA-256 inventory.

- Three initial regressions failed before implementation: dropped raw evidence,
  positive integration of legacy quality records, and loss during annotation
  projection. Twenty-eight JVM cases were added overall, including coverage,
  region-bound validation, model eligibility and legacy serialization checks.
- The JVM gate passes 1,110 tests with zero failures/errors and two optional
  real-data demo skips. All installed HMMER/AntiFam cases ran; their inputs and HMM
  are explicitly synthetic, not a biological AntiFam benchmark.
- Python verification passes 552 benchmark tests (12 skips in discovery), with
  all nine opt-in Nextflow cases then executed successfully against the repaired
  JAR. Three optional benchmark skips remain. All 41 trainer and nine input-support
  tests pass. These totals include four new Python regression tests.
- Seven actual CLI invocations reproduce the old loss, preserve both the HMMER
  domain claim and AntiFam observation after repair, retain the flag after
  projection, and exercise direct integration with known/unknown coverage,
  successful no-hit scanning and a genuine missing-GA tool failure.
- The mixed fixture's raw claims and final annotations increase from one to two.
  Its functional integration TSV remains byte-identical to the old engine; the
  added record is the preserved quality observation, not an extra function.
- An initial CLI checker incorrectly required identical record order across
  concurrent predictor runs. Its failure is retained; the corrected comparison
  checks the complete multiset of records, preserving duplicate multiplicity.
- Two intermediate fixture tests exposed an unnecessary empty export field. The
  writer now emits quality metadata only when present, preserving existing
  no-quality fixture bytes. These failures are also retained.

Before engine: `94db839ff7c8eaf1bcb03b9d28d710d2cdbe45342f5941dfc520fe3809b0469c`.
Repaired engine: `cb74fc910024c3d199a38e94a1a5d24466cb9740fc11fd88ff72071b67579122`.
Both frozen JARs remain in the content-addressed local engine cache. Older engines
may drop this new evidence; use the repaired engine to replay the new raw claims.
No production model was refitted or replaced in this repair.
