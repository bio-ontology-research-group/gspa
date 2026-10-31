# Family-excluded predictor training and production verification

The verified family-excluded homology transfer has its own predictor parameter
identity. The current homology-trained models correctly reject that identity, so
they cannot supply a learned comparison for this arm. This run fits separate
logistic and tree artifacts using actual family-excluded production features.
It does not silently reuse the homology models or weaken the identity gate.

This is a reference-support experiment on an already examined, reference-overlapping
cohort. Its closed-world negative rows are unconfirmed annotations, not measured
losses of function. It does not establish an untouched test, biological calibration,
role accuracy, or a trained null outcome. Neither model is promoted on the basis
of these results.

## Frozen recipe

`benchmark/v3/run_family_excluded_training.py` separates registration, production
export, fitting, and test evaluation. Registration binds the reconstructed claims,
exact predictor identity, engine, all reference/context inputs, supervision policy,
splits, source modules and training runtime. It copies the previous selection
rules: group weighting, five grouped folds, seed 42, and calibration-selected
maximum-F1 display threshold. It does not copy fitted coefficients or thresholds.
Both artifacts are frozen before either test evaluation. Existing paths are not
overwritten; changed registered inputs invalidate subsequent stages.

The predictor is DIAMOND 2.1.11, with parameter SHA-256
`910a3d2f261916bdfa56ba7b0509ac034d16383c5f0bac523eda25cd0b071912`.
Its full version/reference/adapter identity comes from the
[verified predictor reconstruction](v3-predictor-contract-replay.md).
The production engine is
`74d9e8fc8c170adcb8d222cc0fb2ba82199b877e63a51efffca40f8dac512405`.
The existing genomic/taxonomic split, family masking, label policy, locus
certificates and cohort context settings remain fixed.

Production export preserves all **272,366** unique claim identities across all
**18 genomes**, with the correct family-excluded predictor contract on every
candidate. The supervision policy yields 2,444 rows:

| Split | Rows | Positive reference labels | Explicit negatives | Unconfirmed closed-world negatives |
|---|---:|---:|---:|---:|
| Training | 1,683 | 183 | 1 | 1,499 |
| Calibration | 345 | 29 | 0 | 316 |
| Test | 416 | 34 | 0 | 382 |

These are the eligible rows of the new candidate set; they are not a claimed
preservation of the larger homology-supervision table. Candidate losses and truth
coverage remain in `training_table_manifest.json`. The only complete locus group
is still a training competition; there are no complete null groups and no
calibration/test locus groups. Both artifacts therefore keep the null layer
explicitly untrained.

## Evaluation and limits

Both models support 299 of 416 test rows. They abstain on 115 rows whose profiles
lack both classes in training and calibration and on two outside the frozen score
ranges. The following metrics describe only supported reference-support rows;
these are not comparable to a different method's complete test universe.

| Model | Brier | AUROC | Log loss | Calibration threshold | Predicted positives | Precision | Recall over all 34 test positives |
|---|---:|---:|---:|---:|---:|---:|---:|
| Logistic | 0.09082 | 0.74292 | 0.31702 | 0.29819 | 4 | 0.25000 | 0.02941 |
| Tree | 0.08837 | 0.75170 | 0.30472 | 0.21581 | 42 | 0.28571 | 0.35294 |

The low precision and limited recall do not support an accuracy or calibration
claim. Matching predictor identities and observed feature ranges are necessary
eligibility checks, not proof of calibrated biological probabilities. Complete
paired stage comparisons, genome uncertainty, experimental negatives/roles/nulls,
and independent family/genome validation remain required. No threshold or model
was selected from these test results.

## Production checks

The existing `check_context_mode_models.py` verifies all 2,444 table rows against
both model readers: **4,888** JVM/Python comparisons, identical eligibility, and
maximum probability error `4.996e-16`. Each model supports 1,779 table rows.

`benchmark/v3/check_family_excluded_models.py` then invokes full production batch
inference for both models on all 18 genomes. It checks every candidate against the
frozen export and Python model reader, including input features, predictor
provenance, null probabilities for unsupported records, unchanged fallback rows,
and agreement of final TSV scores with the learned report when projection is permitted.
Taxon-forbidden scores remain reported but cannot overwrite the prior annotation;
the verifier requires that explicit disposition and unchanged output. The first
verifier incorrectly required every supported score to be projected. Its failed
run and source are retained; a focused regression now covers this policy.

Both full production runs are verified: **544,732 candidate/model comparisons**
across 36 genome/model runs. Each model scores 196,621 candidates, projects 193,883
scores, and retains 2,738 scores without projection because of taxon constraints.
The remaining 75,745 candidates abstain: 70,536 lack both classes in their
training/calibration profile, 3,538 have an unseen profile, and 1,671 are outside
frozen score ranges. Every fallback is unchanged. Maximum absolute probability
error is `6.106e-16` for logistic and `1.110e-16` for the tree.

The ablation runner previously selected the homology artifact for both claim
regimes. It now selects the corresponding frozen model directory inside each
regime, retaining both logistic and tree comparisons. An explicit
`GSPA_COHORT_MODEL_DIR` override remains supported; predictor eligibility still
applies. Three added tests cover arm construction, overrides and invalid modes.
This wiring enables a valid family-excluded comparison but does not execute or
certify the entire multi-arm evaluation.

Ten focused lifecycle/auditor tests cover incorrect source identity, missing
or duplicate candidates, changed registration inputs, test-row exclusion from
training, model freezing before evaluation, altered report provenance, incorrect
final predictions, and modified unsupported fallbacks. The complete benchmark
suite passes **696 tests with 13 optional skips**. No JVM implementation changed
in this run; the bound engine is the previously validated 1,327-test build.

## Artifacts and reproduction

The working run is `~/.cache/gspa-v3-data/family-excluded-training-v1/`.
The registration and protocol sit at its root; `default/` contains the complete
export, training table, two frozen models and evaluation. The independent parity
run is `family-excluded-training-v1-parity/`; the complete production replay is
`family-excluded-model-replay-v2/` under the same cache root.

Use the Python environment with the registered runtime and a fresh output path:

```bash
export TMPDIR="$HOME/tmp"
export JAVA_TOOL_OPTIONS="-Djava.io.tmpdir=$TMPDIR"
PYTHON="$HOME/tmp/gspa-v3-testenv/bin/python"
"$PYTHON" benchmark/v3/run_family_excluded_training.py prepare \
  --template "$HOME/.cache/gspa-v3-data/predictor-support-v2/default" \
  --claims-replay "$HOME/.cache/gspa-v3-data/predictor-contract-v1" \
  --jar "$HOME/.cache/gspa-v3-data/engines/74d9e8fc8c170adcb8d222cc0fb2ba82199b877e63a51efffca40f8dac512405.jar" \
  --out "$HOME/tmp/family-excluded-new"
"$PYTHON" benchmark/v3/run_family_excluded_training.py export --out "$HOME/tmp/family-excluded-new"
"$PYTHON" benchmark/v3/run_family_excluded_training.py train --out "$HOME/tmp/family-excluded-new"
"$PYTHON" benchmark/v3/run_family_excluded_training.py evaluate --out "$HOME/tmp/family-excluded-new"
```

The [immutable archive](../../benchmark/v3/results/family-excluded-training-v1/README.md)
contains both frozen models, registrations, complete exports, fitting/evaluation
receipts, both-reader parity, full production reports, the failed first verifier,
source snapshots and test logs. Every stored and decompressed payload is hashed.
The original run binds the pre-default-wiring runner; that exact source snapshot
is preserved. Current benchmark defaults are recorded separately. Replaying the
original receipt requires its bound snapshot; a fresh registration binds current
sources. Changed source inputs must never be silently accepted against an old
receipt.

The overall twelve-item acceptance plan remains open. These artifacts enable
further matched comparisons; they do not certify those comparisons or satisfy
missing biological evidence by themselves.

## Report-auditor hardening

A later review reproduced a verifier defect: tolerance comparisons could accept
NaN probabilities, and supported candidates could carry `score_status=unsupported`.
`check_family_excluded_models.py` now requires finite real probabilities within
[0,1], excludes booleans, and checks the supported score status. The check covers
Python predictions, reports, exported features, projection outcomes and both
baseline/final TSVs. Two added tests produce five failing subcases before the
fix and pass afterward; the complete benchmark gate passes 712 tests with 13
optional skips.

A separate scan rehashes the exact previously archived outputs before and after
checking all **544,732** candidate/model records across **36** reports. All
**393,242** supported probabilities, feature probabilities and numeric outcomes
are finite and valid; **151,490** unsupported records retain explicit nulls.
All prediction/feature rows and statuses pass. This repairs the verifier without
changing the original production outputs or their scientific interpretation.
The working receipt is
`~/tmp/gspa-observed-gap-models-v1/family-excluded-finite-audit.json`.

The [immutable verifier addendum](../../benchmark/v3/results/family-excluded-finite-audit-v1/README.md)
retains the old/new checkers, failing/passing tests, complete scan receipt and
bindings to the original output archive.
