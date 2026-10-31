# Matched joint-cohort inference with observed search gaps

Status: all four export arms, frozen model fitting, test evaluation, table
reader parity, full production replay and matched production-score analysis are verified.
This advances items 1 and 5; it does not establish biological benefit or complete
any of the twelve acceptance requirements.

The main historical ablation runner writes no gap inputs, so its requested RLC
stage has no eligible work. The complete [pathway replay](v3-cohort-pathway-replay.md)
now supplies verified observations for all 18 genomes. Its local-context checks
are extended here to matched **joint cohort** runs and both transfer regimes.
`run_observed_gap_cohort.py` is the explicit registered entry point for this
comparison. The historical ablation command is not certified by these results.

## Registered design

The four arms pair homology and family-excluded claims, each with and without
the **521 observed gap records**. Both arms of each pair use the same complete
candidate universe, sequence evidence, predictor identities, genome layouts,
orthogroups, clade holdout, reaction network, reference files and integration
rules. All 18 genomes are retained. The export experiment reads no truth labels and
fits no models; separately registered fitting follows verified export. Registration binds all inputs before execution; the first version covered 1,223 files.

The verified gap files come from complete proteome searches and preserve
reaction/pathway identities, weak-hit flags and EC-to-GO mappings. The preparer
checks all historical replay input/output hashes and recomputes reaction-group
eligibility from complete observed tables. Incomplete acquisitions, changed
files or invented gap groups fail. The observed-gap arm requests RLC explicitly;
its ledger must record executed reaction context for every genome.

Each arm exports every candidate. Comparison requires exact candidate identities
and predictor provenance, unchanged non-gap features, correct gap availability,
unchanged sequence likelihoods and exact retention of the observed gaps.
All functional score changes are retained, including decreases and descriptive
threshold crossings; none is labelled correct or harmful without independent truth.
RLC feature rows and suggestions are counted separately; suggestions are not
silently treated as revised functional predictions. Fresh execution receipts
bind commands, inputs, outputs and the shared production engine:
`74d9e8fc8c170adcb8d222cc0fb2ba82199b877e63a51efffca40f8dac512405`.

## Interpretation and outstanding work

The searches use overlapping references. In particular, “family-excluded” names
the primary homology transfer regime; it does **not** mean that pathway-search
references exclude those families. These outputs cannot establish unseen-family
performance. Missing search hits remain hypotheses, not experimentally proven
absence. EC-to-GO mapping correctness and biological pathway activity are not
established by numerical execution checks.

Existing learned artifacts were fitted without these gap-availability profiles.
They must retain their unsupported-input behavior. Separate registered refitting and calibration-only policy selection are now
complete. Full production model validation and matched held-out reference-support
comparisons are complete for the four registered arms. Independent biological
calibration, roles and null outcomes remain separate unmet requirements.

Five focused regression tests cover complete acquisition, altered outputs,
invented gap groups despite resealed hashes, candidate coverage, predictor
identity, unintended feature changes and four-arm registration. The full benchmark
suite passes **716 tests**, with **13 optional skips**. No JVM source changed.

## Reproduction

Use fresh output directories and the registered immutable engine. The template
commands provide existing cohort settings; each run binds its current actual
inputs anew. No old fitted model or evaluation result enters this export.

```bash
export TMPDIR="$HOME/tmp"
export JAVA_TOOL_OPTIONS="-Djava.io.tmpdir=$TMPDIR"
PYTHON="$HOME/tmp/gspa-v3-testenv/bin/python"
CACHE="$HOME/.cache/gspa-v3-data"
"$PYTHON" benchmark/v3/run_observed_gap_cohort.py prepare \
  --homology-template "$CACHE/predictor-support-v2/default" \
  --family-template "$CACHE/family-excluded-training-v1/default" \
  --claims-replay "$CACHE/predictor-contract-v1" \
  --pathway-replay "$CACHE/cohort-pathway-replay-v2" \
  --jar "$CACHE/engines/74d9e8fc8c170adcb8d222cc0fb2ba82199b877e63a51efffca40f8dac512405.jar" \
  --out "$HOME/tmp/observed-gap-cohort-new"
"$PYTHON" benchmark/v3/run_observed_gap_cohort.py run \
  --out "$HOME/tmp/observed-gap-cohort-new"
```

The corrected working run is `~/.cache/gspa-v3-data/observed-gap-cohort-v2/`.
The first run, `observed-gap-cohort-v1/`, is preserved separately.
A completed `verification.json` is required before accepting production results.

## Preserved first verifier failure

Both homology batches completed successfully, but the first verifier incorrectly
required unchanged functional prediction TSVs. The full workflow's existing
`gap_filling` prior consumes the supplied gaps, unlike the earlier local replay
with priors disabled. The failure is preserved with its source and registration;
no production outputs or gap inputs were changed to satisfy the assertion.

A separate diagnostic compares all 484,202 original homology predictions. It
finds **202** score changes: **198 increases and four decreases**, with **197**
upward crossings of 0.5 and none downward. The largest absolute change is
**0.343139**. Candidate identities and sequence likelihoods remain identical.
These are heuristic changes, not demonstrated corrections or calibrated calls.
The full changed-row inventory is retained. The corrected verifier measures
these effects explicitly and repeats the four registered production arms in a
fresh directory.

Three added comparison tests cover both directions and crossings, unchanged
sequence evidence, identity loss, and nonfinite/inconsistent probabilities.
`fit_observed_gap_models.py` separately registers frozen supervision and original
selection rules, preserves every supervision row and split, and freezes both
model kinds for both transfer modes before any test evaluation. Four lifecycle
tests cover train/calibration-only fitting, missing artifacts, changed inputs,
changed runtime and unverified templates. An initial test used the wrong mock
argument index; the corrected fixture passes, and the failure log is retained.
Real fitting and table-reader verification subsequently completed; the results
below supersede that pending state. The failure archive remains unchanged.

The [failure checkpoint archive](../../benchmark/v3/results/observed-gap-cohort-checkpoint-v1/README.md)
preserves the failed verifier, both completed first-run production batches,
every bound input, the full score-change diagnostic, corrected registration,
source snapshots and test logs. It does not contain a completed corrected
comparison or certify the pending fitted models.

## Verified four-arm export and frozen models

The corrected comparison completed all **72 genome/arm executions**. Both transfer
regimes retain all 521 gaps and exact candidate/provenance/sequence identities.
Only the three declared gap feature columns change. Each regime has identical
candidate coverage with and without gaps:

| Transfer regime | Candidates per arm | RLC rows | Suggestions | Score increases | Decreases | Upward crossings of 0.5 |
|---|---:|---:|---:|---:|---:|---:|
| Homology | 484,202 | 7,408 | 371 | 198 | 4 | 197 |
| Family-excluded homology | 272,366 | 2,241 | 116 | 95 | 11 | 95 |

Neither regime has downward crossings of 0.5. These are descriptive effects of
the existing rules, without an accuracy interpretation. Genomes with no RLC rows
remain included and their stage execution is checked.

`fit_observed_gap_models.py` preserves all **4,262 homology** and **2,444
family-excluded** supervision rows and their original splits. It freezes all four
models before evaluating test rows. Both readers agree across **13,412** table
row/model comparisons; maximum absolute error is **4.44e-16**. Test results are
conditional on supported inputs and target reference-annotation support:

| Regime | Model | Supported/total test rows | Brier | AUROC | Precision | Recall over all test positives |
|---|---|---:|---:|---:|---:|---:|
| Homology | Logistic | 462/608 | 0.122885 | 0.759259 | 0.451613 | 0.186667 |
| Homology | Tree | 462/608 | 0.127387 | 0.719533 | 0.368421 | 0.186667 |
| Family-excluded | Logistic | 299/416 | 0.090293 | 0.742802 | 0.500000 | 0.029412 |
| Family-excluded | Tree | 299/416 | 0.088369 | 0.751697 | 0.285714 | 0.352941 |

Precision uses each model's original calibration-selected threshold. The apparent
50% precision of the family-excluded logistic model means **one supported positive
among two calls**. No superiority or biological calibration follows from these
point estimates. All four models still lack an honest trained null outcome.

The current artifacts are `~/.cache/gspa-v3-data/observed-gap-models-v1/` and the
two `observed-gap-models-v1-parity-*` directories. These are working artifacts;
the earlier failure archive does not certify them. The complete production check
uses `check_observed_gap_models.py`, registers all four commands before execution,
and compares every candidate with the matching gap-enabled rule baseline. It
checks unsupported fallback, taxon projection, preserved gap inputs and Python
reader parity. Two focused tests cover all four command plans, source bindings,
correct baseline/model selection, unchanged original plans and invalid genome
coverage. Its fresh output is `observed-gap-model-replay-v1/` and is now verified.

## Matched production-score analysis

`analyze_observed_gap_models.py` is prepared for the verified production replay.
It compares no-gap rules, gap-enabled rules and both gap-enabled learned models
on the complete frozen supervision table within each transfer regime. Final
production probabilities are used, including unsupported-input fallback and
taxon enforcement. Missing or duplicate supervised functions fail; EC prefix
normalization preserves protein identifiers containing pipes.

Separate max-F1 cutoffs are selected from all calibration rows, with conservative
ties, and frozen before test metrics. A common 0.5 cutoff is also reported. The
analysis retains every test row, reference-label basis, Brier/ECE, confusion
counts and reference-relative corrections/harm. Paired genome resampling uses
10,000 draws and seed 20261027 for the mean per-genome Brier difference; these
four previously examined genomes provide descriptive sensitivity, not a
calibrated population uncertainty guarantee. Unconfirmed reference rows remain
explicit and are not experimentally proven negatives. No promotion or novel-family
claim is permitted.

Four tests verify calibration/test separation, conservative threshold ties,
complete identity matching, nonfinite rejection, correction/harm direction,
equal genome weighting, undefined single-genome intervals and retention of all
reference labels. An initial threshold-tie fixture was not actually tied; that
failure is retained and the corrected fixture passes. The complete benchmark
gate passes 716 tests with 13 optional skips. Analysis execution required the verified complete production replay and has
now completed with all registered test rows retained.

## Supervision does not adequately cover observed gaps

The complete production exports contain gap-linked candidates, but the existing
supervision table provides very little evidence about them. Both gap-linked
indicator columns have the following nonzero coverage (the availability column
is one for every row):

| Regime | Split | Gap-linked positive labels | Gap-linked unconfirmed labels |
|---|---|---:|---:|
| Homology | Train | 0 | 2 |
| Homology | Calibration | 3 | 1 |
| Homology | Test | 0 | 2 |
| Family-excluded | Train | 0 | 2 |
| Family-excluded | Calibration | 1 | 1 |
| Family-excluded | Test | 0 | 0 |

Consequently, successful refitting and numerical replay cannot establish learned
recovery of real missing enzymes. In particular, the family-excluded table-test
metrics contain **no gap-linked examples**. The two training examples are not
experimentally confirmed negatives. Model selection and thresholds are unchanged;
this coverage audit does not promote any model or redefine the acceptance scope.
The complete inventory is `~/tmp/gspa-observed-gap-models-v1/supervised-gap-feature-coverage.json`.

The [verified export/training archive](../../benchmark/v3/results/observed-gap-model-training-v1/README.md)
now preserves the complete corrected export, all bound inputs, four frozen models,
table-reader parity and tests. It predates the completed full production/model
comparison documented below. Its archiver initially mishandled relative parity-output paths; the
failed attempt is retained, and the resumed archive verifies all 1,375 checksums.

## Complete production replay and matched test results

All **1,513,136 candidate/model comparisons** pass across four models and 18 genomes.
Each homology model scores 343,944 candidates: 339,118 are projected and 4,826
remain constrained by taxonomy. Each family-excluded model scores 196,621:
193,883 are projected and 2,738 remain taxon-constrained. Unsupported candidates
retain exact rule fallback (**140,258 homology**, **75,745 family-excluded**, per
model). Maximum probability error is **5.56e-16**; all supplied gap records remain
unchanged. This replay uses the frozen 74d9e8fc engine above, before the separate
community-result validator repair; these runs do not request that community stage.

The registered final-score analysis includes **all 608 homology test rows** and
**all 416 family-excluded test rows**, with unsupported fallback and taxon effects
included. There are four test genomes per regime. Its calibration-selected
cutoffs apply to these complete final-score populations and can differ from the
supported-only model advisory cutoffs reported earlier.

| Regime | Final-score arm | Brier | Precision | Recall | F1 |
|---|---|---:|---:|---:|---:|
| Homology | No-gap rules | 0.107801 | 0.297297 | 0.440000 | 0.354839 |
| Homology | Gap rules | 0.108086 | 0.273810 | 0.306667 | 0.289308 |
| Homology | Gap logistic | 0.105613 | 0.341463 | 0.186667 | 0.241379 |
| Homology | Gap tree | 0.109034 | 0.250000 | 0.266667 | 0.258065 |
| Family-excluded | No-gap rules | 0.080812 | 0.060000 | 0.176471 | 0.089552 |
| Family-excluded | Gap rules | 0.080812 | 0.060000 | 0.176471 | 0.089552 |
| Family-excluded | Gap logistic | 0.075744 | 0.333333 | 0.029412 | 0.054054 |
| Family-excluded | Gap tree | 0.074362 | 0.272727 | 0.352941 | 0.307692 |

The homology no-gap rules have the highest F1 in this comparison. Relative to
those rules, logistic/tree classification changes yield respectively **56/24**
and **38/33** reference-relative corrections/harms at calibration-selected cuts.
Family-excluded counts are **92/5** and **91/23**. These counts use the declared
closed-world policy: **533** homology and **382** family-excluded test labels are
unconfirmed, not experimentally established negatives.

The paired mean-genome Brier deltas versus no-gap rules are:

| Regime | Arm | Mean delta | Descriptive 95% genome-resampling range |
|---|---|---:|---|
| Homology | Gap logistic | -0.001697 | [-0.007863, 0.005080] |
| Homology | Gap tree | 0.003421 | [-0.003357, 0.015463] |
| Family-excluded | Gap logistic | -0.004890 | [-0.008046, -0.001441] |
| Family-excluded | Gap tree | -0.006613 | [-0.010525, -0.003327] |

Both homology ranges include zero. The family-excluded improvements concern a
test population with **no gap-linked examples**, so they do not establish a gap
benefit. At the common 0.5 cutoff, homology logistic has nine corrections and
11 harms relative to no-gap rules; homology tree makes no calls. Both
family-excluded models recover zero test positives at that cutoff. All results,
including these unfavorable outcomes, remain reported. No model is promoted.

Complete receipts and scored rows are in `observed-gap-model-replay-v1/` and
`observed-gap-model-analysis-v1/` under the data cache. Neither full-system
biological acceptance nor unseen-family generalization follows from this reused,
sparse reference-support cohort.

The [complete production/analysis archive](../../benchmark/v3/results/observed-gap-model-production-v1/README.md)
preserves all bound inputs, candidate comparisons, final-score policies and
results, including the sparse gap-supervision audit. All 1,753 checksums pass.
