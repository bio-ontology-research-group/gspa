# Exact-copy evidence within heterogeneous strain clusters

The [contract-v3 control](v3-cohort-contract.md) passed for its 18 real genomes,
but a broader mechanism replay failed: copying one of two distinct observations
in a strain cluster raised smoothed support from 0.5 to 0.678571. Capping the
cluster's total weight did not stop its internal raw average being reweighted.

## Implementation

`CohortEvidenceIdentity` fingerprints each genome's complete observed input.
Genome, protein and contig identifiers do not define its identity. Coordinates,
strand, contig boundaries and multiplicity, assembly uncertainty, family/subfamily
membership and weights, all candidate scores/mappings/provenance, intrinsic role
evidence, family-search coverage and organism reaction sets do. Unknown and known
empty inputs remain distinct. Equal input representations are not a claim of
identical sequences or independent ancestry.

`CladeModel.units` can now use these identities to select one representative per
observed input **within each existing strain cluster**. It never collapses
different clusters or clades. Representatives are selected from each eligible
subset, so a filtered/held-out member cannot suppress an available equivalent
member. Distinct observations still share the cluster's capped weight equally.
Every weighted support, coverage, background, copy-number, joint-context and
permutation calculation uses the same policy.

Candidate records and descriptive raw genome/background-window counts are
preserved. The singleton flag uses distinct informative observations rather than
counting copies. No protein or alternative function is removed from output.
The cohort engine records `cohort-v3` and the policy
`complete-observed-input-within-strain/1` in provenance.

The boolean permutation pools are sampled by success count and size, with a
canonical member order. This is the same per-pool sampling distribution, but
renaming or reordering contigs/proteins cannot change the random realization.
Members of a strain cluster still share a randomization stream. Realized
permutation features can differ from earlier engines, so models are refitted
after a fresh export rather than reusing old coefficients as updated results.
The feature schema remains `gspa.context-features/3`: no columns were added or
removed, and their intended meanings are unchanged by the bug repair.

## Regression and production proof

Four new heterogeneous-cluster cases failed before the repair (1, 3 or 5 copies,
two seeds). They now pass. Additional tests cover arbitrary identifier renaming,
input order, multi-contig layouts, paralogs, unknown metadata, preservation of all
candidate records, filtered representatives and independent strain units. Nineteen
input perturbations verify that genuinely different observations are not collapsed.

The original archived replay run through the repaired production JAR returns
support 0.5 at 0, 1, 3 and 5 copies. Every extracted model column is identical;
the duplicate-induced logistic-margin change is zero. This is mechanism evidence,
not a calibrated biological prediction.

JAR SHA-256:
`94db839ff7c8eaf1bcb03b9d28d710d2cdbe45342f5941dfc520fe3809b0469c`.
Fresh 18-genome production export preserves all 484,202 candidate identities,
source profiles and non-cohort feature values. No columns change. Permutation
excess changes in 5,043 records and p-values in 893; all other exported values
remain identical on this cohort.

Both models have been refitted with all 4,262 supervision rows unchanged.
Frozen table evaluation covers 462 of 608 test rows; 146 profiles remain
unsupported. Logistic/tree AUROC is 0.75764/0.71544 and Brier score is
0.124491/0.127707. Reference-support fractions at their frozen display thresholds
remain 15/32 and 21/59. Both null layers remain untrained because the data still
contain zero fully-negative and one single-positive training competition.
These are reference-support metrics, not experimentally established accuracy.

Both models completed all three production arms: 118 genome/arm outputs in total.
Independent receipt/hash checks pass. Across all 18 original genomes, each model
has 346,269 supported numeric candidates and 137,933 unsupported candidates.
Every supported probability is exactly unchanged under duplication; no candidate
or support status changes. This complements the heterogeneous-cluster regression
and its actual-JAR replay; it does not replace them.

The separate whole-genome calibration policy still selects logistic EC 0.44
(25/27 supported calibration calls), abstains for logistic GO, and abstains for
both tree outputs. Both policies were frozen before test scoring. Held-out EC
top-1 reference agreement is 14/37 for logistic and 13/37 for the tree in each
arm. Logistic intact/duplicated arms each make six supported and six unverified
calls; the shuffle makes seven supported and five unverified calls (two new
supported calls, one lost). Thus calibration support still does not transfer,
and the control does not establish an intact-neighborhood advantage. Unverified
calls are not experimentally demonstrated false assignments.

The full JVM gate passes 1,082 tests with six optional skips. The benchmark gate
passes 551 tests with three optional skips, including nine executed live Nextflow
checks against the new engine and refitted model. The migration
runner now also accepts a same-contract replay, requiring exactly the same
columns and unchanged non-cohort values; five migration-integrity tests pass.

## Scope and reproduction

The archive is `benchmark/v3/results/cohort-evidence-dedup-v1/`. The original
engine, models and failing replay remain in `cohort-contract-v3/`, whose completed
114-file checksum inventory has been verified. Re-run the archived Groovy
replay with either pinned engine to compare the mechanism directly.

Current cohort runners default to this archive's refitted models. Use
`GSPA_COHORT_MODEL_DIR` for an explicitly selected alternative. Historical results
remain tied to their recorded engine, source snapshots and models. The production
re-export/refit is reproducible with `migrate_cohort_feature_contract.py`'s separate
`export`, `train` and `evaluate` stages, using the previous
`cohort-contract-v3/table/training_table_manifest.json`, previous archived models,
the pinned new JAR and a fresh output directory.

The real cohort is previously examined development data. Exact observed-input
deduplication does not identify near-duplicate assemblies with changed coordinates,
scores or annotations, nor infer missing strain metadata. It does not establish
biological validity, null/role training, complementary-gap recovery or all of
recommendation 6. Those acceptance requirements remain open.
