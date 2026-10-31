# Context modes with observed gaps

Status: all eight export arms, twelve model fits, 216 production runs and final
analysis are complete. Independent job 10223 verifies all 6,981 bound files,
reconstructs the scores and statistics, and evaluates all 37 experimental-reference
EC proteins, including missing candidates. See the [complete results and retained
archive](v3-context-comparison-results.md). No broad context benefit or superiority
over CLEAN is established. No model is promoted and all twelve whole audit items
remain open.

This extends the [observed-gap comparison](v3-observed-gap-cohort.md) to all
three context modes. The earlier [context comparison](v3-context-modes.md)
did not include these observed gaps and used the earlier predictor identity.
The present comparison retains both transfer regimes, all 18 genomes, and all
521 observed gaps. Family exclusion applies to homology transfer; overlapping
pathway-search references still preclude an unseen-family claim.

## Registered execution

`benchmark/v3/run_observed_gap_context_modes.py` registers eight export arms:
a default compatibility control and UNIVERSAL, PATHWAY and ORGANISM context
for each transfer regime. All three explicit modes receive the same verified
pathway groups and organism metadata. Only the context mode differs. The
default arm checks the new engine against the preceding observed-gap outputs.

The registration binds commands, manifests, raw claims, gap records, context
inputs, engine, Python sources and runtime. Each arm preserves candidate and
sequence-source identity, predictor contracts, roles, quality evidence, sequence
likelihoods and non-cohort numerical features. Gap availability must remain one
and supplied gap records must match exactly. The default control additionally
requires identical numerical features and functional prediction bytes.
Functional changes in the other arms are measured and retained; gaps can change
rule scores, so unchanged functional predictions are not assumed.

Refits must keep all 4,262 homology and 2,444 family-excluded supervision rows' labels,
splits, identities and provenance. Logistic and tree models are fitted separately
for each mode using the preceding train/calibration settings. All twelve models
must be frozen before the test evaluation stage begins.

Execution root: `~/.cache/gspa-v3-data/observed-gap-context-modes-v1`.
The engine is the immutable JAR with SHA-256
`b23a7f6d94eb06c375e63ad344cfbf4c15d7644e33c1848766a36c4287cfd090`.
The pipeline records export, training, evaluation and both transfer-regime JVM
table-parity stages separately. All four homology arms contain 484,202 candidates
each; all four family-excluded
arms contain 272,366 each, across all 18 genomes. Default functional predictions
remain unchanged. All 521 observed gaps are retained in each applicable arm.

The original local PATHWAY execution was interrupted solely to move computation
off the laptop. Job 10162 then failed before script execution on node007. Job
10163 completed all exports on node002, but its train stage could not find
`benchmark/v3/registry/training_cohort_registry.json`. No model was fitted: the
failure occurred while loading the registry, before a table or model was written.
The complete export tree and original failure logs remain preserved at
`/data/hohndor/gspa-v3/observed-gap-context-modes-v1-relocation/jobs/10163/`.

The original registration bound prior training artifacts but omitted direct
supervision-file dependencies from the relocation file list. A supplemental
registration binds the seven actual inputs before training. Six match their
prior training-manifest SHA-256 values. The registry hash is absent from those
table manifests, but it **is already bound by the parent registration**
(`observed-gap-models-v1/registration.json`, SHA-256
`c1628d8bf03ad81c95daf0a7a36446ffac8f29d78086fffdcc53b948123ddcd9`).
Its registry hash is
`ba86d3c44ea7771e5edc702e5a2c158a0266c68733bc7e627d4c23dbb7b00e38`.
The previous wording suggesting no prior registry binding was incorrect.
Both original table manifests are also bound. The supplement is verified before
and after every resumed stage; source, recipe, splits and model settings are
unchanged. The final analysis registration binds the parent registration rather
than directly listing these seven files; the independent collection audit checks
the full chain and all seven files explicitly.

Job 10176 restores the original immutable runtime/package and the verified
exports, adds the supervision supplement, and resumes at training. It does not
rerun exports. Persistent resumed records are under
`/data/hohndor/gspa-v3/context-resume-v1/`. Its saved training tree contains only
new/changed files relative to the job-10163 export base, named in
`unchanged-export-base.txt`; reconstruction requires that base plus the new tree.
This preserves the original evidence and avoids duplicating unchanged exports.
No model execution runs on the laptop.

## Production and final-score checks

`check_observed_gap_context_models.py` prepares all twelve replay commands
before executing any. Each model is paired with its own context baseline.
The intended scope is 216 genome/model runs and 4,539,408 candidate/model
comparisons. Every candidate is checked for Python/JVM probability agreement,
support or abstention, rule fallback, taxonomy effects, features and provenance.
Every run must retain its supplied gaps. At most two independent arms execute
concurrently, each with four audit workers. Job 10176 completed all 216 runs and
4,539,408 candidate/model comparisons. Its two table-reader checks also completed
all 40,236 comparisons across the twelve models. The original workflow records
every stage as completed with exit zero; the retained job finished at
2026-10-27T21:25:23+03:00. Independent collection and expanded EC evaluation are complete in the
[final audit](v3-context-comparison-results.md); execution and numerical agreement
do not establish biological accuracy.

`analyze_observed_gap_context_models.py` compares nine final-score arms per
transfer regime: three contexts, each with rules, logistic and tree scoring.
It requires complete verified production artifacts and the same supervision
identities across contexts. It keeps every frozen test row, including unsupported
learned candidates that fall back to rules. No common-support intersection is
used. Each PATHWAY/ORGANISM arm is compared against UNIVERSAL with the same
scorer; each learned arm is also compared against rules in its own context.

All max-F1 thresholds are selected from calibration rows and saved before the
first test metric is calculated. Results also use a common threshold of 0.5.
The report includes confusion counts, Brier score, ten-bin ECE, reference-relative
corrections and harms, and paired genome Brier differences. The seeded 10,000
resamples are descriptive uncertainty across four previously examined test
genomes, not calibrated population confidence. Both transfer regimes remain
separate. The complete analysis was written by job 10176, with 608 test rows
in homology and 416 in family-excluded transfer. These are selected supervision
rows, including only 15 and five EC-positive pairs respectively. The independent
collection audit additionally evaluates all 39 experimental-evidence EC assertions
across all 37 reference proteins in the four test genomes, including unemitted
functions, using the existing frozen thresholds without retuning.

The completed historical run uses comparison schema 1. The [schema-2 reporting
repair](v3-ec-candidate-coverage.md) now includes full-reference EC populations
in routine reports and directly binds supervision dependencies. Its complete
replay preserves all historical scores, cutoffs and original summaries.

## Validation and limitations

Six export/training tests cover complete planning, invariant inputs, finite
features, preserved supervision and freezing all models before test evaluation.
Two replay-plan tests cover all twelve baseline/model/context pairings, source
bindings, genome coverage and overwrite rejection. Four final-score analysis
tests cover complete populations, missing/nonfinite inputs, calibration-only
thresholds, matched contrasts, harms, genome weighting and saving all thresholds
before test metrics. The complete benchmark gate passes 730 tests with 13 optional skips.

The table-parity checker previously allowed NaN to evade its tolerance comparison.
It now requires both Python and JVM probabilities to be finite real numbers,
excluding booleans, within [0,1]. Two tests reproduce five invalid cases before
the repair and pass afterward. The completed table-parity runs use this repaired
checker; earlier receipts retain their original source identity.

This is a reference-support comparison on a previously examined cohort.
Unconfirmed rows are not experimentally proven negative labels. Each transfer
regime has only two gap-linked training rows, both unconfirmed; the family-excluded
test table has no gap-linked rows. Models retain an untrained null layer.
Numerical agreement or improved aggregate reference scores cannot establish
missing-enzyme recovery, biological calibration or independent accuracy.
