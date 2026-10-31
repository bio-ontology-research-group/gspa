# Context modes with observed gaps

Status: local execution was stopped at Robert’s request. Continuation job
10163 is running on Unimatrix node002 (12 CPUs, 48 GB, eight-hour limit). The replay and final-score comparison
implementations pass their focused tests; their full production results are
pending. No model is promoted and no whole audit item is accepted.

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

All 4,262 homology and 2,444 family-excluded supervision rows keep their labels,
splits, identities and provenance. Logistic and tree models are fitted separately
for each mode using the preceding train/calibration settings. All twelve models
must be frozen before the test evaluation stage begins.

Execution root: `~/.cache/gspa-v3-data/observed-gap-context-modes-v1`.
The engine is the immutable JAR with SHA-256
`b23a7f6d94eb06c375e63ad344cfbf4c15d7644e33c1848766a36c4287cfd090`.
The pipeline records export, training, evaluation and both transfer-regime JVM
table-parity stages separately. The homology default and UNIVERSAL controls have completed, each retaining
484,202 candidates. The default control preserves all functional predictions.
The PATHWAY JVM was interrupted with exit 143 solely to move computation off
the laptop; its partial outputs are retained and are not certified. The remaining
stages now continue remotely under the same registered file paths and inputs.
A filesystem view binds the transferred Python/Java runtimes and input tree
inside node-local scratch. No further model execution is scheduled on the laptop.
Job 10162 failed before script execution on node007; its scheduler receipt is
retained. The transferred JDK has a pre-existing dangling optional `src.zip`
link; runtime files and all registered inputs still require verification.
Persistent package/results root:
`/data/hohndor/gspa-v3/observed-gap-context-modes-v1-relocation`.

## Production and final-score checks

`check_observed_gap_context_models.py` prepares all twelve replay commands
before executing any. Each model is paired with its own context baseline.
The intended scope is 216 genome/model runs and 4,539,408 candidate/model
comparisons. Every candidate is checked for Python/JVM probability agreement,
support or abstention, rule fallback, taxonomy effects, features and provenance.
Every run must retain its supplied gaps. At most two independent arms execute
concurrently, each with four audit workers. The replay has not yet run.

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
separate. This analysis has not yet run.

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
the repair and pass afterward. The upcoming table-parity runs use this repaired
checker; earlier receipts retain their original source identity.

This is a reference-support comparison on a previously examined cohort.
Unconfirmed rows are not experimentally proven negative labels. Each transfer
regime has only two gap-linked training rows, both unconfirmed; the family-excluded
test table has no gap-linked rows. Models retain an untrained null layer.
Numerical agreement or improved aggregate reference scores cannot establish
missing-enzyme recovery, biological calibration or independent accuracy.
