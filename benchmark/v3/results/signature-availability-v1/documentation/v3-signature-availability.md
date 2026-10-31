# Untestable reaction signatures remain unknown

This repair addresses an inference defect relevant to requirements 5 and 6.
It does not complete either requirement or establish biological accuracy.
The [evidence archive](../../benchmark/v3/results/signature-availability-v1/README.md)
preserves failing regressions, both production engines' outputs, exact commands,
feature comparisons, model refits and test logs.

## Defect and corrected behavior

The cohort engine called a present family `NEGATIVE` when its observed genome
window had no matching anchors, even if the reaction had no network or interaction
neighbors to test. This produced effective negative evidence and an adjusted
support estimate from missing reaction knowledge. A configured minimum of two
anchors caused the same problem for a signature with only one possible anchor.
The learned bridge also exposed these cases as available local context.

A signature now requires at least `minAnchors` possible distinct reactions to be
testable. A present family with an untestable signature has `UNKNOWN_CONTEXT`,
zero informative effective N, and no adjusted interval, background enrichment or
permutation evidence. Family membership weights and copy observations survive.
Certified absence and unobserved absence remain independently determined from
the family search and assembly metadata.

The raw feature report retains descriptive layout and anchor counts, records
`signature_testable`, and omits local positivity/enrichment when untestable.
The shared learned bridge retains layout availability and withholds its
local-context group. A testable signature with a fully observed window and no
anchors still supplies negative context; independently supplied interaction
edges can make an otherwise isolated reaction testable. No protein/function
candidate or raw sequence score is removed.

Provenance records `signatureAvailabilityPolicy` as
`minimum-distinct-anchors-or-unknown/1`, and raw candidate states distinguish
`reaction_signature: testable` from `untestable`. The existing feature-contract
version 3 columns remain unchanged; existing availability indicators carry the
corrected missingness. Models must be refitted for changed feature values.

## Regression and production verification

Seven new cohort cases cover empty and insufficient signatures, positive and
negative independently supplied interaction contexts, weighted family presence,
copy counts, and certified versus uncertified family absence. Three cases fail
before the repair. An eighth workflow test reproduces the local-availability
handoff defect. Its intermediate failure is retained. The existing fixture's
count of locally testable candidates changes from four to three because its
fourth reaction is isolated; its candidate count stays four.

The full JVM gate passes **1,229 tests**, zero failures/errors, six optional
skips. Separate live HMMER/AntiFam suites execute all four previously skipped
external-tool checks successfully (39 cases total, zero skips), leaving the
two unrelated optional JVM cases unexecuted. The benchmark gate passes **600
tests**, with three optional skips; all nine live Nextflow cases execute on
the repaired JAR.

Ten actual before/after CLI calls cover empty signatures, curated positive and
negative controls, an unattainable anchor minimum, and a single genome without
peers. They retain all candidate identities and byte-identical functional TSVs.
The connected controls remain identical. The initial CLI verifier incorrectly
treated a second one-neighbor reaction as testable under a two-anchor minimum;
that failed diagnostic and its source are preserved before the corrected replay.

Engines:

- Before: `1f0a7da89894896b05626d386879be494fd1e8d1ecbacba3b048e600086befd6`.
- After: `6757088967769a6d9204f50b34a7e96ebf42a177e4ae50588db27a562cdc3dd7`.

## Complete cohort effect and model refit

Fresh exports from both engines cover all eighteen fixed genomes. The latest
before-engine export exactly reproduces the preceding model-training features.
The matched before/after comparison preserves **484,202 candidate identities**,
source identities, all non-cohort feature values and every functional prediction
TSV. Feature values change for **4,780 candidates**. Local context becomes
unavailable for **4,779**; informative peer evidence becomes unavailable for
**4,064**. These are corrected availability decisions, not lost proteins or calls.

Alternative reactions remain candidates. The existing rule selecting a
representative by weighted support can choose a different reaction when false
negative evidence becomes unknown. This also changes some other cohort columns,
including local anchors and the representative's reaction-expansion count.
The bounded production ambiguity examples contain 178 changed representatives;
this is not an exhaustive count of changed alternatives. The complete feature
comparison is exhaustive.

Both logistic and tree models are refitted with all **4,262 supervision records**,
source identities, locus groups and original splits preserved. Training settings
and calibration-only threshold rules are unchanged. No threshold is selected
from these test results. The target remains reference-annotation support under
the declared closed-world policy, not experimentally established enzyme function.

| Model | Supported test rows | AUROC | Brier score | Reference-supported display calls |
|---|---:|---:|---:|---:|
| Logistic | 462 / 608 | 0.755645 | 0.124982 | 14 / 31 |
| Tree | 462 / 608 | 0.719533 | 0.127387 | 14 / 38 |

Both null layers remain untrained: the training table has zero complete
no-correct-candidate groups and one complete single-positive group. The 146
unsupported test rows are withheld. This previously examined cohort is not a
new blind test, and these metrics do not establish improved biological accuracy.

All **118 whole-genome control outputs** now complete and pass independent
verification for both refitted models. Across all eighteen original genomes,
each model has 346,269 numeric probabilities exactly invariant to adding five
copies of a known strain; 137,933 candidates remain unsupported, with no missing
identities or support transitions.

The whole-genome logistic EC policy selects 0.45 on intact calibration data:
25/26 calls have reference support. Test intact/duplicate calls remain six
supported and six unverified; shuffled neighborhoods yield seven and five.
The shuffle's net supported-call change is +1, with a paired genome-bootstrap
interval of −3 to +3. Logistic GO and both tree policies abstain. EC top-1
reference agreement remains 14/37 for logistic and 13/37 for the tree in every
arm. The failure of calibration support to transfer and the lack of evidence
for a neighborhood benefit remain visible; no accuracy acceptance follows.

Current cohort benchmark runners default to the refitted models in this archive.
`GSPA_COHORT_MODEL_DIR` selects an explicit alternative. Pair these models with
the repaired engine using `GSPA_JAR`. Historical outputs remain bound to their
original models and source snapshots, including the source before this default
selection change. Full input verification of an old run requires its frozen
sources, rather than assuming that the live checkout will remain unchanged.

## Scope still open

This makes the existing reaction-context feature semantics more faithful; it is
not the requested pathway-only/universal/organism-specific ablation. Complete
archaeal pathway acquisition, the reference-overlap and search-reproducibility
issues, independently labelled loci/roles/nulls, matched context comparisons and
biological acceptance remain open. No whole audit item is accepted.
