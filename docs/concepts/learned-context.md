# Learned context integration (foundations)

Status: **mechanism and training/inference code only. There is no trained production
model in the repository, and nothing here is evidence of improved annotation accuracy.**
The shared workflow (`IntegrationWorkflow`) now calls these classes: the `learned_integration` stage runs
last over the original raw claims, the `cohort_context`, `metabolic_experiments` and `community_metabolism`
stages fill `ExternalContext`, and `--feature-export-dir` writes production rows from the same extractor. What is
still missing is the real training run (see [What real training needs](#what-real-training-needs)); until a
trained artifact exists `scoring_mode` stays `legacy_heuristic`. The stage and the training-row export are
described in [learned-stage.md](learned-stage.md); stage order, ledger and open gaps in
[context-workflow.md](context-workflow.md).

This replaces hand-set boosts, floors and count-dependent caps with a model whose weights
are learned from labelled candidates, calibrated on a separate split, and able to move a
candidate's probability **down as well as up**.

| Piece | Location |
| --- | --- |
| Feature contract (typed, versioned, documented) | `gspa.integration.learned.FeatureContract` |
| Feature extraction from claims + state + optional context | `ContextFeatureExtractor`, `ExternalContext` |
| Artifact loader and inference (logistic or trees), null outcome | `LearnedContextModel`, `CandidateScore`, `LocusAssignment` |
| Trainer / exporter / reference predictor | `scripts/train_context_model.py` |
| Tests | `gspa-core/src/test/groovy/gspa/integration/learned/`, `scripts/tests/test_context_model_training.py` |

## 1. Feature contract

`FeatureContract.CURRENT` (`gspa.context-features/1`) lists every column with its group and
description; `toJson()` writes the file the trainer reads. An artifact is bound to the contract
by id, fingerprint (over ordered names + groups) and the name list, so extractor, trainer and
inference cannot silently drift.

Groups: `FUNCTION_IDENTITY`, `ROLE`, `BASELINE_SEQUENCE`, `RUN_COVERAGE`,
`CONTEXT_NEIGHBORHOOD`, `CONTEXT_COHORT`, `ASSEMBLY`, `FBA_CONDITIONAL`. **Baseline sequence
features are separate from context features.**

Rules that apply to every feature:

* **Missing is not absent.** Anything that can be unavailable is a value column plus an
  availability/reported column; unavailable values are 0 and the flag is 0. "No neighbourhood
  was measured" and "neighbours were measured and none support" are different vectors.
  Predictor failures are visible as coverage (`seq.cov.*`, from the predictor status manifest),
  not as evidence against a function.
* **No double counting.** Predictor evidence enters as the strongest raw score per source and
  per correlation group (`homology`, `ml_protein_seq`, `structure`, `motif`, `domain_specific`);
  repeating the same hit, or adding weaker hits from one tool, changes nothing. Counts are of
  distinct sources/groups, never of hits.
* **No posteriors as evidence.** Nothing is read from posterior log-odds, likelihoods, pinned
  floors, orthogroup consensus posteriors, or the table-calibrated probabilities. Claims whose
  evidence type is context-derived (operon, metabolic, dark-matter, reaction-local, cross-genome,
  ranker, localisation, region/viral) create a candidate but never enter the baseline block.
* **Gap information is an indicator**, not a Bayes factor: `ctx.gap.function_listed` /
  `ctx.gap.gapseq_guessed_listed` say that an open gap lists the function; what that is worth is
  learned. The gap prior is never relabelled as an RLGC likelihood ratio.
* **FBA is conditional.** `fba.*` values describe behaviour under the tested conditions only and
  are *suppressed* (flagged `fba.suppressed_objective_forced`) when the objective itself
  required the candidate, so choosing such an objective cannot count as confirmation.
* **Only what exists is used.** Neighbourhood (reaction-graph) support, cohort recurrence, FBA,
  predictor statuses and assembly completeness arrive through `ExternalContext` from real
  analyses. If they are not supplied the model sees "unavailable"; nothing is fabricated.

Local context is bounded and candidate-specific: the same operon plus genes within ±10 kb.
Adding unrelated candidates elsewhere cannot change a candidate's vector (tested).

### Candidates and roles

Candidate identity is `protein | type | function | role`. Roles (`CATALYTIC`, `ACCESSORY`,
`TRANSPORTER`, `REGULATOR`, `COMPLEX_SUBUNIT`, `UNSPECIFIED`) are separate candidates scored
separately. Roles are derived conservatively (EC claim → catalytic; GO complex term →
complex subunit) unless structured evidence supplies them. **Complex membership is never
equated with catalytic function.** Every claimed function yields a candidate; no score
threshold removes candidates.

## 2. Two different probabilities

1. **Calibrated per-candidate correctness** (`CandidateScore.probability`):
   `sigmoid(a·margin + b)` with `margin` from the logistic model or the tree sum and `(a, b)`
   fitted by Platt scaling on the calibration split. It answers "is this function correct for this
   protein in this role?". It is a *marginal, independent* quantity: several candidates can all be
   high (paralogues, isoenzymes, multi-reaction proteins). It is never normalised across proteins,
   functions or roles. `aboveDisplayThreshold` is advisory (threshold chosen on calibration) and
   does not gate retention.
2. **Conditional locus ranking with an explicit null** (`LocusAssignment`): for the candidates of
   *one* reaction-locus group (same function and role, distinct proteins),
   `P(i) = exp(s·m_i) / (exp(b0 + w·g) + Σ_j exp(s·m_j))` and
   `P(none) = exp(b0 + w·g) / (same denominator)`, where `m_i` are the candidate margins and `g`
   six group features (best/second margin, log(1+n candidates), any sequence support, fragmented
   assembly). This is conditional on "exactly one of these loci or none". It sums to one **within
   the group only**; a singleton is not normalised to 1 (its null probability is whatever was
   learned); groups for different reactions, and the per-candidate quantity above, are not tied
   together. Groups with several correct candidates are excluded from this layer's training rather
   than forced into one-of-n labels.

The null layer has **no default**. It exists in an artifact only if the training split contained
enough groups where *no supplied candidate was correct* and enough single-positive groups
(`--min-null-groups`, `--min-pos-groups`) and the calibration split could validate it.
Otherwise the artifact says `locus_model.trained = false` with the reason and
`assignLocus` throws `ModelUnavailableException`. An artifact that is missing, unfrozen,
bound to a different contract, or a mechanism-test fixture makes the learned stage fail or be
skipped; no invented weights, boosts or null probabilities are substituted.

## 3. Artifact

JSON, schema `gspa.learned-context-model/1`: `intended_use` (`production` or `mechanism_test`),
contract id/fingerprint/feature names, `kind` (`logistic`: train-only normalisation + weights;
`gbdt`: init, learning rate, shallow trees as node arrays), Platt `calibration`, advisory
`display_threshold` (rule and value, chosen on calibration only), `locus_model`, and `training`
provenance: `frozen_before_test`, hash of the **train + calibration** rows, per-split identity
(counts and hashes of example ids, groups, genomes), enforced disjointness keys, weighting,
hyper-parameters and CV scores, `synthetic`. No timestamps: identical inputs give byte-identical
artifacts. Test metrics are in a separate `test_report_<kind>.json`, which also carries the
whole-input-file hash and the artifact hash it evaluated.

Both model kinds share this contract; the trees are trained by the trainer itself (Newton-boosted
depth-2/3 regression trees on the weighted logistic loss, numpy only), so no extra dependency is
needed. The comparison logistic vs trees is reported on the test split after freezing; the
"primary" suggestion in `training_summary.json` uses the train-only CV score.

## 4. Training protocol (enforced in code)

```
/tmp/gspa-v3-testenv/bin/python scripts/train_context_model.py train \
    --table features.jsonl --contract feature_contract_v1.json --out-dir out/ \
    --model both --weighting group
```

Input rows (JSONL with a `features` object, or TSV with one column per feature):
`example_id`, `genome_id`, `group` (family/species cluster), `split` ∈ `train|calibration|test`,
`label` (0/1, candidate correct), optional `locus_group_id` (rows of one reaction-locus group),
`role`. Rejected: any missing provenance field, duplicate ids, missing/unknown/non-finite feature,
a group or genome appearing in two splits (default keys `group,genome_id`; `--disjoint-keys group`
relaxes only the genome check, and the choice is recorded), a missing split, a single-class
calibration split. Nothing is ever called out-of-fold on trust.

* `train()` receives train and calibration rows only. Normalisation, weights, hyper-parameters
  (grouped k-fold CV **inside train**), Platt calibration and the display threshold cannot see
  test rows; test rows are read by `evaluate_frozen()` after the artifact is serialised and hashed.
  Tests confirm that flipping or perturbing the test split leaves the artifact byte-identical,
  and that calibration labels move only calibration/threshold, not weights or trees.
* `--weighting group|genome|genome_group|none` weights each example by the inverse size of its
  group/genome (mean 1), so replicated or large families cannot dominate the fit.
* The locus layer's margins on training groups are out-of-fold (grouped CV); its single
  regularisation constant is chosen on the calibration split.

## 5. What the tests do and do not show

Mechanism tests (hand-set weights or planted synthetic data) verify that the code does what the
design requires: repeated hits add nothing; 100 unrelated candidates change nothing; unfavourable
neighbourhood context lowers the learned probability and keeps lowering it (no floor); a weak
candidate is scored and returned; a missing candidate gets non-zero null probability and a
singleton is not forced to 1; complex-subunit and catalytic candidates stay distinct; paralogues
in the same operon are visible as competitors; JVM inference reproduces the Python reference to
1e-12 for logistic, trees and the locus layer; the trainer refuses leaky or unprovenanced input.
The committed fixture (`gspa-core/src/test/resources/learned/mechanism_fixture_*.json`) comes
from planted synthetic data, is marked `mechanism_test`/`synthetic`, and the production loader
refuses it. **None of this measures accuracy on real genomes**; any such claim needs the data
below and the frozen-protocol report.

Regenerate the fixture after a contract change: run `FeatureContractSpec` (it writes
`gspa-core/build/learned-fixtures/feature_contract_v1.json`), then
`train_context_model.py make-mechanism-fixture --contract <that file> --out-dir
gspa-core/src/test/resources/learned`; copy `feature_contract_v1.json` from the same build
directory next to the fixtures, then re-run `FeatureExporterSpec` and copy
`build/learned-fixtures/export_fixture*.jsonl` into the same resources directory (they embed the
contract fingerprint and the model hash).

## What real training needs

Per candidate row (gene × function × role), exported by the production engine rather than a
parallel formula, with the feature columns above plus provenance:

1. **Labels from independent experimental truth**, disclosed per source: exact-EC attribution
   (Price with accession correction), time-split / homology-controlled GO, reference genomes with
   curated functions, and reaction-masking experiments for locus groups. Synthetic or model-derived
   truth is flagged and kept out of headline numbers.
2. **Negative and null groups**: for locus groups, groups where the masked reaction's true locus
   was *not* among the supplied candidates (needed to learn the null outcome), plus groups with
   several correct loci.
3. **Declared splits** by family/species cluster **and** genome (family-heldout and genus-heldout
   variants), with calibration separate from test, plus the per-row `genome_id`/`group`.
4. **Context inputs** as `ExternalContext`: reaction-graph neighbourhood support/conflict,
   cohort recurrence with incomplete-assembly counts (clade-aware, duplicated strains controlled),
   conditional FBA per condition with the objective-forced flag, predictor completion statuses,
   assembly completeness, and explicit roles where structured evidence exists.
5. **Matched baselines** on the same rows (sequence-only vector, per-stage context) so gains are
   attributable, with genome-level uncertainty, correction/harm counts and calibration plots.
