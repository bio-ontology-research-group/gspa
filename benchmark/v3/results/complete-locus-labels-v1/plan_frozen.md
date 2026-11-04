# Frozen plan: honest complete-locus labels for the learned stage (row 4), 2026-11-02

Status: frozen before any completeness count or score is read. The sha256 of this file as frozen is
recorded in `benchmark/v3/results/complete-locus-labels-v1/PLAN_SHA256` before the first count.
Deviations are appended, dated, under "Deviations"; the frozen sections are never edited.

Scope: the MISSING element of acceptance row 4 ("honest complete locus labels"; the null outcome
untrained) in `tasks/v3-acceptance-audit-2026-11-01.md`. Established public data only (UniProtKB/
Swiss-Prot, GOA, Price et al. 2018 reannotations). AllTheBacteria is not used.

## 1. Definitions

**Competition.** For one genome, one function identity (function type, function id, role) and one
candidate source, the set of every protein that the source proposes for that function at any score.
The source is fixed per genome set (section 2). Which proteins belong to a competition is read from
the claim/export files only; no score is read to decide membership or eligibility.

**Label.** A protein-function statement with experimental support:

- positive: (a) a Swiss-Prot GO annotation with an experimental GO evidence code (EXP, IDA, IPI,
  IMP, IGI, IEP, HTP, HDA, HMP, HGI, HEP), propagated to is_a/part_of ancestors with the pinned
  `go-basic.obo` exactly as the production table builder does; (b) a Swiss-Prot four-level EC with
  ECO:0000269 on its `DE EC=` line; (c) for the 22 Price genomes, a Price-148 complete EC
  (Price et al. 2018 reannotation, counted as experimental by instruction).
- negative: an experimental statement denying exactly that function for that protein:
  (a) a GOA `NOT` row with an experimental evidence code (from the audited June 2026 GOA inventory,
  `gaf-negative-audit-v1`), exact term, no propagation; (b) a Swiss-Prot `CC` sentence carrying
  ECO:0000269 that denies the activity, mapped to an EC or GO MF id; (c) the reviewed negatives
  already accepted in the role-truth registry (`benchmark/v3/role_truth_cases.jsonl`, status accepted).
  Every negative that enters a complete competition is read by the curator (this agent) against its
  sentence/row and, where retrievable, its PubMed record; it is accepted only when it denies the
  candidate function itself (not one substrate, not a different process involvement). Dispositions,
  accepted or rejected, are recorded.
- Unconfirmed is never negative. A protein without a label for that function makes the competition
  incomplete; it is never counted as a negative.
- A positive and a negative for the same protein/function is a conflict; the competition is excluded
  and reported.

**Provenance (required for every label in a complete competition).** Accession, source release
(Swiss-Prot 2025_03 flat file sha256 `a6907c4d…`; GOA snapshot 2026-06-17 for NOT rows and for the
PMIDs of GO positives; Price et al. 2018 for Price-148), evidence code, and a PMID. GO positive
PMIDs are taken from GOA rows (same accession, same annotated term, same evidence code); because the
local GOA file is no longer on disk, the GOA annotations for exactly the accessions in candidate-
complete competitions are downloaded once on ws from the EBI QuickGO REST service into an empty
directory, pinned by sha256 and download date. A label without a PMID is not usable; its competition
is incomplete. Price-148 labels carry PMID 29769716 (Price et al., Nature 2018).

**Complete competition.** Every member has a usable label. Classes: `single_positive` (exactly one
positive), `all_negative` (null true), `multi_positive` (isoenzymes; complete but excluded by the
trainer's one-of-n-or-null layer). Size class: `singleton` (one member) or `contested` (two or more).

## 2. Genomes, candidate universes and splits

| set | genomes | candidate universe (score-blind) | split |
|---|---|---|---|
| cohort | 18 label-rich genomes | family-excluded DIAMOND production feature export (`family-excluded-training-v1`, 272,366 candidates; the universe the current default learned models score) | `benchmark/v3/results/training-cohort/splits.tsv`: 11 train / 3 calibration / 4 test, disjoint on genome, group, genus, family |
| model | E. coli K-12, B. subtilis 168, M. tuberculosis H37Rv | role-fenced Swiss-Prot DIAMOND claims (`rolefence-claims`) | none exists; reported as a separate held-out, development-examined stratum, never counted as calibration/test |
| Price | 22 Price genomes | release-candidate claims of TopEC, MAPred, EnzBERT (one competition set per predictor) | genus-blocked split `price-calibration-context-v1/row1/split.json` (sha256 `95020999…`): 7 calibration / 15 test |

Dense Price predictors (CLEAN, GraphEC, DeepProZyme) propose every EC on most proteins; their
universes cannot be complete and are not processed. Only the cohort competitions are scorable by the
current learned models (the model and Price universes carry a predictor contract the models do not
support); those two sets contribute to the label set and its counts only.

## 3. Steps

1. Build universes; build positive and negative label tables with provenance; compute completeness
   (score-blind). Record the registry (`complete_competitions.jsonl`) and its sha256 in the results
   folder (`REGISTRY_SHA256`) before any score is read.
2. Cohort: write `labels.tsv` (existing cohort labels plus accepted negatives) and
   `locus_groups.tsv` (every complete cohort competition with 0 or 1 positives), rebuild the table
   with the unchanged production builder (`prepare_training_data.py build-table`, family masking and
   split audit on), and record which declared groups the builder certifies.
3. Fit with the unchanged trainer (`scripts/train_context_model.py`), both kinds (logistic, gbdt),
   the family-excluded recipe (group weighting, five grouped folds, seed 42, max-F1 threshold) and the
   default null gates (≥20 null and ≥20 single-positive training groups, ≥5 and ≥5 calibration
   groups). Artifacts frozen (sha256) before the test evaluation.
4. Evaluate on test with the trainer's `evaluate`. Independently of whether the null layer trains,
   report on every complete calibration/test competition the frozen family-excluded models'
   (`models/default/`) behaviour: top-1 correct in contested single-positive groups, and whether any
   member reaches the frozen calibration threshold in all-negative groups.

## 4. Pass/fail criteria

| id | criterion |
|---|---|
| C0 | Plan sha256 recorded before the first count; registry sha256 recorded before the first score read |
| C1 | Every label in every complete competition has accession, release, evidence code and PMID (labels lacking any are excluded first; the number excluded is reported) |
| C2 | No negative without an accepted reviewed source (validator over the registry; 0 violations) |
| C3 | ≥ 20 complete competitions usable by the null layer (0 or 1 positive), all 43 genomes |
| C4 | ≥ 5 of them in calibration or test genomes of the existing splits (cohort or Price) |
| C5 | Strict reading, reported separately: C3 and C4 restricted to contested competitions (≥ 2 members) |
| C6 | Every declared cohort group certified eligible by the production builder, except losses to family masking, which are reported |
| C7 | Null layer trained by the unchanged trainer at its default gates (no gate, split or threshold change) |
| C8 | If C7 passes: test locus NLL below the uniform baseline log(n+1) averaged over test groups. If C7 fails: not applicable, and the descriptive evaluation of step 4 is reported |

If a target is missed, the report states exactly how many competitions exist per set, split, class
and size, what blocks each incomplete one (counts of unlabelled members), and what data would be
needed.

## 5. Outputs

`docs/benchmarks/v3-complete-locus-labels.md`; `benchmark/v3/results/complete-locus-labels-v1/`
(plan copy, PLAN_SHA256, REGISTRY_SHA256, registry, labels, locus groups, builder manifest, model
artifacts and evaluation, reviews, SHA256SUMS); code in `benchmark/v3/complete_locus/`; one progress
line in `tasks/v3-supervision-plan-2026-10-28.md`. Raw work in
`/mnt/data2/leechuck/gspa-v3-data/complete-locus-labels-v1/`. Compute on ws (small); no Shaheen.

## Deviations
