# Predictor identity and observed-range eligibility

This implements the identity/distribution-support gate following the
[verified DIAMOND reconstruction](v3-predictor-contract-replay.md). It does not
establish biological calibration or complete any whole audit requirement.

## Contract and inference

Every usable baseline claim may carry `metadata.predictor_contract` with exactly
eight string fields: `schema` (`gspa.predictor-contract/1`), `tool`, `tool_version`,
`tool_sha256`, `reference_sha256`, `parameters_sha256`, `adapter_sha256`, and
`score_model`. Hashes are lowercase SHA-256 strings. The reconstruction binds
these declarations to real execution receipts; a declaration alone is not a
cryptographic attestation or an accuracy guarantee. Other adapters need their own
verified identities and corresponding training support.

The extractor retains each normalized source/contract pair in
`sequence_contracts`, independently of the strongest-per-source numeric features.
Identical pairs are deduplicated. A weaker usable claim with a missing or malformed
contract retains an explicit null entry and prevents scoring. Derived context
claims and unusable sequence scores do not create predictor coverage. Candidates
with no sequence evidence have an empty identity list and still require their own
supported profile. Exports, training tables and learned reports retain this field.

`gspa.model-input-support/3` profiles bind exact source identities, predictor
contracts, function/role/availability indicators and context mode. Each profile
needs both classes in both training and calibration. For every non-indicator
`seq.*` feature, inclusive lower/upper bounds are the extrema observed in training
and calibration for that same profile. Test rows cannot widen them. Out-of-range
inputs abstain; no clipping, extrapolated probability or partial-locus ranking is
substituted. Unknown predictor identity, unsupported profile and out-of-range
scores have distinct reasons. Quality-evidence exclusions remain in force.

Legacy support schemas 1 and 2 remain readable but abstain with
`model_predictor_contract_support_unavailable`. Missing input identities produce
`candidate_predictor_contract_unavailable`; a changed contract on an otherwise
known profile produces `unseen_predictor_contract`. Range violations produce
`sequence_score_outside_training_calibration_range`.

The numerical feature contract is unchanged. Trainer version 4 binds predictor
metadata into the training digest. The gate is only a necessary coverage test:
matching hashes and marginal ranges do not establish a matching joint
distribution, taxonomic transfer, calibration within the retained subset, or
correct biological function. Synthetic mechanism models remain explicitly exempt
from production eligibility and cannot be used as calibration evidence.

## Reproduction

`benchmark/v3/replay_predictor_support.py` has separate `export`, `train` and
`evaluate` stages. Export requires a verified reversible claims replay and checks
every candidate, numerical feature and legacy prediction. The sole normalization
of existing metadata makes an omitted historical context mode explicitly
`UNIVERSAL`; every other change fails. Training preserves all supervision and
splits, refits the original settings, and requires exactly unchanged normalization,
coefficients/trees, calibration and display thresholds before freezing models.
Evaluation reads those frozen files. Historical models and data remain intact.

`benchmark/v3/check_predictor_support.py` runs packaged old/new CLIs on every EC
claim from the selected real genome and deliberate version/reference/score-model,
missing-identity, weaker-unidentified-evidence and range faults. It also checks
legacy models against the new engine. These are mechanism interventions, not new
biological labels. `check_context_mode_models.py` compares JVM and Python coverage
and probabilities over complete migrated supervision tables.

The first export check rejected the historical omitted-mode field. Its failed
receipt is retained; the explicit-default conversion has a regression test. The
CLI verifier also retains its initial option-name error and its overly strict
assumption that every real candidate had layout support. Two unsupported real
candidates are retained in the corrected complete-universe comparison.

## CLI verification

All **40** packaged old/new/baseline CLI calls pass over **1,921** real EC claims.
The old engine scores 1,919 candidates, including after identity changes. Two
candidates lack a represented layout-availability profile and already abstain.
With verified identities, the new engine scores 1,895 candidates and additionally
withholds 24 outside the frozen sequence-feature ranges. Every retained score is
unchanged; maximum Python/JVM probability error is below `1e-12`.

For both logistic and tree models, changing the version, reference hash or score
transform with identical vectors withholds all 1,919 previously supported cases.
The two existing profile gaps remain unsupported. Missing identity and an added
weaker unidentified claim withhold all 1,921. Setting every raw score to zero
withholds 1,919 by range plus the two profile gaps. Loading each legacy model in
the new engine withholds all 1,921 by unsupported contract schema. Unsupported
records have no numeric learned probability or locus probability; their baseline
annotation rows remain identical to the matched baseline. No candidate is dropped.

Engines are frozen as SHA-256
`705b865d9591746227c1e8d46475240dcb78426f49f5aa54f1153a695ee944f5`
(before) and
`f80982650f49738235e80e0570ce3dd416270f64818126c456cabe4c48d3c345`
(after).

## Complete cohort migration

Four separate exports (historical default, explicit universal, pathway and
organism modes) each preserve **484,202** candidate identities, every numeric
feature and all legacy annotation bytes. Eight refits preserve **4,262**
supervision records each, including splits and locus membership. All eight retain
exactly unchanged scoring parameters, calibration, display thresholds and complete
locus-model records. Every null layer remains untrained. All models were frozen
before this replay's test evaluation.

| Export mode | Supported table rows per model | Supported test rows | Test abstentions |
|---|---:|---:|---|
| Historical default | 3,175 / 4,262 | 462 / 608 | 146 insufficient class coverage |
| Explicit universal | 3,175 / 4,262 | 462 / 608 | 146 insufficient class coverage |
| Pathway | 3,179 / 4,262 | 461 / 608 | 145 insufficient class coverage; 2 range violations |
| Organism | 3,158 / 4,262 | 460 / 608 | 147 insufficient class coverage; 1 range violation |

All **34,096** JVM/Python candidate comparisons agree on eligibility. Maximum
probability difference is `6.11e-16`. The default and explicit universal runs use
their respective historical context settings; they are preserved separately.
Metrics remain conditional on supported rows and reference-support labels. Changed
coverage is not an accuracy gain, and this previously examined cohort is not a
new blind test.

The full JVM suite passes **1,307** tests with six optional skips. The Python
scripts suite passes **286** with sixteen optional skips. The support-specific
Python tests cover all six changed identity components, range boundaries, missing
contracts, per-version class coverage, mutable aliasing, serialization, corrupted
manifests and partial loci. JVM tests additionally exercise extraction, weaker
claims and public reports. The initial full JVM run correctly rejected two stale
export snapshots; their only added data are explicit unknown predictor contracts.

The benchmark gate passes **667** tests with thirteen optional skips; all ten
optional live Nextflow cases then pass separately with the new engine/default
model. They verify feature/report parity, context-mode abstention, cache
invalidation and failure handling. The isolated migration-test invocation also
passes after explicitly setting its local module search path.

## Defaults and archive

Current ablation defaults select the plain model copies under
`benchmark/v3/results/predictor-support-v1/models/default/` and the verified claims
under `~/.cache/gspa-v3-data/predictor-contract-v1/claims/`. The claims are produced
by `replay_predictor_contract.py`, not by inventing identities for old files.
`GSPA_COHORT_MODEL_DIR` and `GSPA_COHORT_CLAIMS_DIR` provide explicit overrides.
Original claim files without identities remain unsupported. The family-excluded
transfer has a different verified parameter identity and no matching learned
training here. The strict ablation runner rejects a wholly unsupported learned
stage instead of evaluating its fallback predictions as a learned method.

The [archive](../../benchmark/v3/results/predictor-support-v1/README.md) retains
both engines, all four exports/refits, all old/new CLI results, parity outputs,
source snapshots, models, tests, logs and initial failed checks. Migration source
snapshots are taken before changing default paths, preserving their exact input
hashes; current-default sources are stored separately. Every gzip payload has an
original hash and a stored hash, and all eight plain model copies are checked
against the frozen refits. The archive's verification report and checksum list
record the final storage checks.

Python visualization tests regenerate the core JVM XML using one focused writer
test. The full JVM gate was therefore rerun after all Python/Nextflow checks and
its complete 1,307-case XML captured under `final-full-jvm:` archive entries.
The initial 434-case XML inventory is retained as an intermediate snapshot; it
does not replace the final complete-suite report. The final JAR hash still
matches the engine used by every migration and CLI check.

Predictor families lacking training support, complete experimental role/locus
truth, calibrated null learning, independence from reference labels and biological
generalization remain open. No whole audit requirement is accepted.
