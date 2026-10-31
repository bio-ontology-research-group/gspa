# Cohort feature contract version 3

The [duplicate-strain replay](v3-cohort-permutation.md) showed that removing the
permutation defect did not remove confidence changes from copied evidence. The
fitted model still consumed raw peer and background-window counts. Two metadata
coverage fractions also used raw genome counts, so copying a known strain could
change apparent metadata coverage without adding independent evidence.

## Contract and compatibility

`gspa.context-features/3` removes these eight columns:

- `cohort.available`
- `cohort.log_n_genomes`
- `cohort.frac_with_orthogroup`
- `cohort.frac_with_function_seq`
- `cohort.frac_with_context`
- `cohort.n_genomes_incomplete`
- `cohort.w.n_peer_genomes`
- `cohort.w.bg_windows`

The weighted cohort block retains its own availability flags. Unknown weighted
context remains distinct from a measured zero. Legacy unweighted summaries remain
readable as provenance but cannot influence the new feature vector. Raw counts
remain in cohort reports. The phylogeny-known and organism-reaction-set-known
fractions now use the same clade/strain weights as other peer summaries.

The new contract fingerprint is
`d9537b61c3e41251e6a9cafba9acb376385b66c58db0ede758f52a4994e58803`.
The loader rejects version-2 artifacts. It does not silently pad, drop or reinterpret
model coefficients. The shared exporter and JVM/Python mechanism fixtures have
been updated, including the context-availability report reader.

Forward migration consists of a production feature export, rebuilding the table
with unchanged supervision and splits, refitting train/calibration models, and
then evaluating frozen artifacts. Preserve old directories. For rollback, use the
old version-2 engine and models together; the previous engine is preserved at
`~/.cache/gspa-v3-data/engines/a792849b501e7ad7a9989488ea3b2ce2bd85242de44b0f279b213272bd5a8547.jar`.
Old results describe their recorded engines and are not recertified against the
new implementation. Never relabel a version-2 artifact as version 3.

## Re-export and training protocol

`migrate_cohort_feature_contract.py` has separate `export`, `train` and `evaluate`
stages. Each refuses overwrite. Production export is verified through the batch
receipt and additional hashes for every exported feature file. Input identity
includes the previous table, original models, labels, splits, policies, groups,
source modules and both feature contracts. Changed inputs invalidate reuse.

The migration compares every old/new candidate identity, genome and source profile.
Exactly the eight declared columns may be removed; no columns may be added, and
no non-cohort feature may change. Changed cohort columns are counted explicitly.
The rebuilt training table must preserve every field except features, including
positive/negative polarity, evidence/references, locus groups and original splits.
Both model kinds use the original training settings, group-wise cross-validation
and calibration rules. Both artifacts are frozen before test evaluation. Test
results cannot choose model settings or thresholds.

This remains the previously examined development cohort. Re-export and refitting
do not supply missing null/role truth or establish generalization. In particular,
one reviewed complete training locus does not support fitting a null layer.

## Current verification

Three new regressions failed before repair: learned vectors changed with raw
counts, legacy unweighted summaries entered scoring, and duplicate known strains
changed metadata coverage. They pass after correction. Contract tests also reject
the genuine frozen version-2 fixture. The full JVM suite passes 1,055 tests with
two optional skips; the Python trainer passes 40 tests and migration integrity
checks pass four tests. JAR:
`ec9be608a606c05fa5d543fddec2434e9d7ce51f185465f7d151a30def33e136`.

Production re-export verifies all 484,202 candidate identities, source profiles and
non-cohort feature values. Exactly eight columns are removed; permutation excess
changes in 5,124 rows and permutation p-value in 896. Both models are refitted with
all 4,262 supervision rows unchanged. Both null layers remain untrained: there are
zero fully-negative and one single-positive training locus groups.

Frozen test evaluation covers 462 of 608 rows; 146 input profiles remain
unsupported. Logistic/tree AUROC values are 0.7576/0.7154, Brier scores
0.12449/0.12771 and log loss 0.41960/0.42774. At their separate frozen display
thresholds, reference-support fractions are 15/32 (46.9%) and 21/59 (35.6%).
These labels include declared closed-world unconfirmed rows; the target is
reference support, not experimentally established functional correctness.

The full benchmark gate passes 543 tests with three optional skips, including
nine live Nextflow checks. The added route test verifies actual rejection of a
version-2 model from the task logs; Nextflow's displayed stderr tail was too short
to expose the underlying contract error. Original failed assertion logs are kept.
Both logistic and tree whole-genome controls completed: all 118 genome/arm outputs
pass independent hash and stage checks. Across all 18 original genomes, the same
346,269 candidates have numeric learned probabilities and 137,933 remain
unsupported. No candidate or support status is lost. Duplicate-arm probability
change is at most 2.22e-16 for logistic and exactly zero for the tree, meeting the
declared absolute tolerance of 1e-12 for this control.

The separate whole-genome call policy requires ten calibration calls and 90%
reference support. Logistic selects EC 0.44 (25/27 supported in calibration) and
abstains for GO. The tree has no qualifying threshold for either type and abstains
for both. These policies differ from the models' display thresholds above.

| Test measure | Logistic intact | Logistic shuffled | Logistic duplicated | Tree, each arm |
|---|---:|---:|---:|---:|
| EC top-1 reference matches | 14/37 | 14/37 | 14/37 | 13/37 |
| Supported EC calls under the whole-genome policy | 6 | 7 | 6 | 0 |
| Unverified EC calls under that policy | 6 | 5 | 6 | 0 |
| GO calls under that policy | 0 | 0 | 0 | 0 |

Calibration support still does not transfer to the logistic test calls. The shuffle
adds two supported EC calls and loses one; these few observations do not establish
an intact-neighborhood benefit. Only three test genomes have EC reference labels.

## Remaining mixed-cluster failure

Subsequent work: the [observed-input identity repair](v3-cohort-evidence-identity.md)
now fixes the archived mechanism failure and adds distinct-member preservation
tests. Fresh model refitting/control validation is recorded separately. The
following paragraph describes the failure in this archive's original engine.

Passing the real control is not general duplicate invariance. A broader synthetic
replay places two genomes with different context evidence in the same strain
cluster, then makes exact copies of one member. Effective N remains 1, but the
weighted, smoothed support fraction rises from 0.5 to 0.678571 with five copies.
The fitted logistic margin's feature contribution changes by 1.59597. The unit's
raw average still counts identical copies repeatedly and shifts its composition.

The failing replay, inputs and output are preserved under
`benchmark/v3/results/cohort-contract-v3/mixed-cluster-failure/`. It is a mechanism
finding, not a calibrated biological prediction. Deduplicating identical evidence
within heterogeneous strain clusters and validating distinct-member preservation
remain required. No whole acceptance item is complete.

## Reproduction

Use a fresh output directory and the recorded input cache:

```bash
export TMPDIR="$HOME/tmp"
export JAVA_TOOL_OPTIONS="-Djava.io.tmpdir=$TMPDIR"
PY="$HOME/tmp/gspa-v3-testenv/bin/python"
MIGRATION="$HOME/.cache/gspa-v3-data/cohort-contract-v3-replay"
$PY benchmark/v3/migrate_cohort_feature_contract.py export \
  --original-manifest "$HOME/.cache/gspa-v3-data/model-input-support-v2-final/table/training_table_manifest.json" \
  --original-models benchmark/v3/results/model-input-support-v1/models \
  --jar "$HOME/.cache/gspa-v3-data/engines/ec9be608a606c05fa5d543fddec2434e9d7ce51f185465f7d151a30def33e136.jar" \
  --out "$MIGRATION"
$PY benchmark/v3/migrate_cohort_feature_contract.py train --out "$MIGRATION"
$PY benchmark/v3/migrate_cohort_feature_contract.py evaluate --out "$MIGRATION"
```

This archive's models remain in `benchmark/v3/results/cohort-contract-v3/models`.
Current local cohort/sequence-loss/agent runners use the subsequent
`benchmark/v3/results/cohort-evidence-dedup-v1/models` refit;
`GSPA_COHORT_MODEL_DIR` can select another explicitly compatible model directory.
Remote frozen jobs retain their recorded JARs and models.
