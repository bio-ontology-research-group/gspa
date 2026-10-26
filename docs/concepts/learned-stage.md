# Learned integration stage and feature export

Status: **mechanism, API and audit trail only. There is still no trained production model in
the repository, so the stage cannot run on real data and nothing here is evidence of improved
annotation accuracy.** Wiring into the shared engine/CLI/workflow is separate work (see
[Wiring](#wiring-for-the-workflow-owner)). Background, feature contract and training protocol:
[learned-context.md](learned-context.md).

| Piece | Location |
| --- | --- |
| Stage | `gspa.integration.learned.LearnedIntegrationStage`, `LearnedStageSettings`, `LearnedStageReport`, `LearnedStageException` |
| Feature export | `gspa.integration.learned.FeatureExporter` |
| Structured context (incl. context-only candidates) | `ExternalContext` |
| Tests | `LearnedIntegrationStageSpec`, `FeatureExporterSpec`, `ExtractorHardeningSpec`, `LearnedContextModelSpec`; Python `ExportedRowsAgreeWithTrainer` |

## 1. The stage

```groovy
def model  = LearnedContextModel.load(new File('learned_model.json'))   // production artifacts only
def stage  = new LearnedIntegrationStage(new LearnedStageSettings(model: model))
LearnedStageReport report = stage.run(rawClaims, state, integratedAnnotationSet, externalContextOrNull)
```

Inputs: the **original raw `EvidenceClaim`s**, the `IntegrationState`, the `IntegratedAnnotationSet`
the earlier stages produced, the frozen model and an optional `ExternalContext`. The set and the
state are updated in place on success and left untouched on failure.

### Contract

* **Replaces, never adds.** The learned calibrated probability becomes the annotation score and
  `ClaimProvenance.finalProbability`; the state's posterior log-odds are re-derived from it. The
  old posterior, its likelihood and its prior contributions are not inputs and are not added to
  (no extra boost). They are kept for audit in `KeyOutcome.priorProbability/priorLogOdds/priorContributions`
  and in annotation metadata `learned_replaced_posterior`. `ClaimProvenance.priorContributions` of
  rescored keys is replaced by the `learned.probability.<ROLE>` / `learned.margin.<ROLE>` entries.
* **Rescoring uses original raw evidence only.** Features are built by the same extractor as the
  export, from raw claims (strongest score per source and per correlation group), state layout,
  operons and gaps and the optional context. Posteriors, likelihoods, floors and previously calibrated
  probabilities are never read, so nothing the earlier stages folded in is counted twice. Tested:
  two runs with wildly different prior posteriors/likelihoods/floors give identical learned scores.
* **All retained candidates are rescored**, including ones with no raw claim (scored with no
  sequence evidence; they do not keep their old posterior). **Weak new candidates** named only by
  structured context (`ExternalContext.candidateKeys()`: any key in `neighborhood`, `cohort`, `fba`,
  `roles` or `candidates`) are added and scored, not thresholded (`addContextCandidates`, default on).
  Claimed functions the earlier stages did not retain are scored and reported but only added with
  `addUnretainedClaims` (use it when the learned stage replaces evidence combination).
* **Rise, fall, retract.** The learned probability moves freely in both directions. With
  `retractBelow` set, a candidate whose learned probability is below it is removed from the
  annotations, provenance and state (recorded as `retracted` with the score that caused it). Unset
  (default), nothing is ever removed for being low.
* **No later floor.** Pinned posterior floors of every key the stage decided are deleted
  (`report.pinnedFloorsCleared`), so a subsequent `IterativeRefiner.clip` cannot lift a learned
  probability. Stages that re-run evidence combination after this one would still overwrite the
  score; the learned stage belongs last.
* **Atomic.** Claims, state, context, features and scores are all validated/computed before the
  first mutation. Any invalid claim, impossible context value (NaN/infinite, cohort subset larger
  than its superset, percentage out of range, unknown predictor status, malformed key), non-finite
  feature or margin, or incoherent set (annotation without provenance) raises `LearnedStageException`
  naming the offending candidate. If applying still fails, the set, annotations, processed ids, state
  posteriors and pinned floors are restored.
* **State/annotation identity.** Identity is the function key `protein|TYPE|function`; the stage
  never pairs by list position. After a run, provenance keys == annotation keys == state posterior
  keys; two isozymes with the same EC on different proteins stay separate. Annotation objects are
  copied, never mutated in place.
* **No model, no run.** `LearnedStageSettings.model == null` → `ModelUnavailableException`; a
  non-production model → `ModelUnavailableException`. Nothing is defaulted.
* **Mechanism fixtures only through the test API.** The public loaders (`load`, `loadJson`) accept
  production artifacts only (and refuse `synthetic` ones). The fixture path
  `LearnedContextModel.parseArtifact(json, contract, true)` is `protected` and `checkModelUse` is
  overridable; both are used only from `LearnedTestSupport` in the test sources, which the CLI
  cannot see. (Groovy does not enforce `protected`; the barrier is that no main-source code calls
  it, and the stage itself rejects non-production models.)

### Roles and the legacy projection

Each `(protein, function, role)` is scored separately. The legacy annotation has no role, so one
role's probability is **projected** onto it:

| Function | Legacy score is the probability of |
| --- | --- |
| EC | `CATALYTIC` only. An EC with no catalytic candidate gets **no** legacy score (`removed_no_projection` if it was retained; reported otherwise). Accessory, regulator, transporter or complex-subunit candidates can never raise it. |
| GO complex term (`state.complexClassification`) | `COMPLEX_SUBUNIT` |
| anything else | `UNSPECIFIED`, or, failing that, the single candidate if there is exactly one |

A non-EC function with several roles and none of these has no legacy projection. There is **no
maximum or union over roles** anywhere, and the legacy score is never described as a union
probability. The probabilities of *all* roles are retained: `report.candidates` (every role, with
`projected` flagged), `KeyOutcome.roleProbabilities`, provenance `learned.probability.<ROLE>`
entries and annotation metadata `learned_role_scores`.

### Marginal probability versus null

The annotation score is the candidate's **marginal** calibrated correctness (independent of other
candidates). If the artifact has a trained locus layer and `computeLocusAssignments` is on, the
conditional locus ranking with the explicit null outcome is computed per `(function, role)` group
and reported in `CandidateRecord.locusProbability` / `locusNullProbability`. It is never written
into the annotation score or provenance probability. Without a trained layer the fields are null and
`report.locusUnavailableReason` says why; nothing is invented.

### Audit

`LearnedStageReport` carries: model id, **SHA-256 of the artifact text**, kind, intended use,
training provenance; feature contract id and fingerprint; per candidate the exact-bits SHA-256 of
the feature vector (`featureDigest`; the full vector too with `retainFeatureVectors`), extractor
notes (what was ignored or pooled), the original supporting raw claims, margin, probability, the
advisory display threshold flag and the locus fields; per function key the prior vs learned
probability and the action (`rescored | added | retracted | removed_no_projection | scored_only`).
Annotation metadata repeats the model hash, contract fingerprint, role, margin, digest and replaced
posterior. The display threshold is advisory and never gates retention.

## 2. Feature export

```groovy
def exporter = new FeatureExporter()                 // no artifact needed
exporter.exportToDirectory(rawClaims, state, ctxOrNull, outDir)         // rows + contract + manifest
exporter.export(rawClaims, state, ctx, writer, retainedSetOrNull)       // streaming to any Writer
new FeatureExporter(model).export(...)                                  // rows also carry the model's score
```

`candidate_features.jsonl` has one row per `(protein, function, role)`:

```json
{"schema":"gspa.context-feature-row/1","contract":{"id":"gspa.context-features/1","fingerprint":"…"},
 "genome_id":"…","protein_id":"…","function_type":"EC","function_id":"EC:1.1.1.1",
 "function_key":"…","role":"CATALYTIC","candidate_key":"protein|EC|EC:1.1.1.1|CATALYTIC",
 "features":{"<every contract column>":0.1234567890123456, …},
 "available_contexts":{"sequence_evidence":true,"predictor_status":false,"operon":false,"layout":false,
   "gap":false,"neighborhood":true,"cohort":false,"assembly":false,"assembly_completeness":false,
   "fba":false,"fba_suppressed":false},
 "notes":["…"], "model":{"id":"…","sha256":"…","kind":"logistic","margin":…,"probability":…,
   "display_threshold":…,"above_display_threshold":false}}
```

* **Same extractor as the stage** (`ContextFeatureExtractor.prepare`, tested by digest equality),
  same candidate universe when you pass the retained set. `model` is present only when a frozen
  model is given.
* **Full precision**: values are written with `Double.toString` (shortest text that round-trips the
  exact double); tests compare bit patterns after parsing. Non-finite values fail the export.
* **No labels.** The exporter reads no truth set and writes no `label`, `split`, `group`,
  `example_id` or `locus_group_id`; functions with no claim get no row, so absence is never
  auto-labelled negative. A dataset builder joins independently curated truth and declared
  family/genome splits by `candidate_key` and supplies those fields. The Python trainer refuses rows
  without them (tested).
* **Streaming**: claims are indexed once; one function's candidates are held at a time. 30,000 rows
  stream through a counting writer in the tests.
* Written alongside: `feature_contract.json` (byte-identical to `FeatureContract.CURRENT.toJson()`,
  the file `train_context_model.py --contract` reads) and `export_manifest.json` (row count, rows
  SHA-256, contract id/fingerprint, model identity if any, `labels_included: false`).
* **Schema agreement with the trainer** (`scripts/train_context_model.py`): `features` is keyed
  exactly by the contract columns, `role` and `genome_id` are the names the trainer reads, floats
  are bit-exact. `ExportedRowsAgreeWithTrainer` loads committed JVM exports with `read_table` and
  checks that the Python reference reproduces the JVM margins to 1e-12.

## 3. Foundation review: defects found and fixed

| Suspected issue | Finding |
| --- | --- |
| Explicit `roleSet` uses the first claim only | Not a defect: type and id are identical for every claim of a function key and roles come from the key. The role path now takes them from the key explicitly; a test shows claim order changes nothing. |
| Source `other` collapses unknown tools and score scales | **Real.** Pipeline predictors `esm2-deepgoplus` and `eggnog` fell into `other` although columns exist (`esm2-go`, `eggnog-mapper`); aliases added. The remaining pooled `other` column cannot be split without a contract change; the pooled tool names are now in the candidate notes and raw scores outside [0,1] (other scales) are ignored with a note instead of entering the maximum. Percent identity outside 0..100 is unavailable. |
| Missingness flags that cannot be true | **Real.** An empty predictor-status map set `status_available=1` ("all clean"); FBA with zero tested conditions set `fba.available=1`; cohort counts with a subset larger than its superset gave fractions > 1; NaN/infinite context values reached the model. All are now unavailable or rejected up front. |
| GO complex classification values | **Real.** Ambiguous (`a`) terms were indistinguishable from homomeric (`h`); any value counted as a complex. Added column `fn.complex_ambiguous` (the contract fingerprint changed; the fixtures and the committed contract were regenerated) and any value outside `h/n/a` is an error. Conflicting GO aspect labels across claims fall back to the ontology/unknown instead of the first label; GAF letters (F/P/C) are accepted. |
| Scopes | Context maps are keyed by function key and apply to every role (documented). Competitors are same function on local proteins regardless of role. The stage now has an explicit type scope (`types`). |
| Sparse/dense candidate caps and O(candidates × claims) | **Real cost issue.** `extract()` rebuilt the whole extraction per call and `extractAll` materialised every vector. Claims are now indexed once (`prepare`), candidates streamed one function at a time, per-candidate work is bounded by the ±10 kb window; a 20,000-candidate test is linear. There are no count-dependent caps in the extractor or stage. |
| Model malformed values | **Real.** Trees with zero nodes or no leaf, missing/empty `trees`, fractional or < −1 feature indices, non-positive calibration slope, learning rate or locus scale, and artifacts without `model_id` produced NPEs/IndexOutOfBounds or silently odd models; they are now `ModelUnavailableException`. A margin that overflows to infinity fails at scoring time. |

## Wiring for the workflow owner

```groovy
// final stage, after existing integration/enforcement; raw claims are the ORIGINAL ones
def model = LearnedContextModel.load(new File(modelPath))                 // throws ModelUnavailableException if absent
def stage = new LearnedIntegrationStage(new LearnedStageSettings(model: model /*, retractBelow: …, addUnretainedClaims: … */))
def report = stage.run(claims, state, integratedAnnotationSet, externalContext)

// training-data acquisition (no artifact):
new FeatureExporter().exportToDirectory(claims, state, externalContext, outDir, integratedAnnotationSet)
```

Needed from the workflow: the original claim list (not re-derived from the set), an
`ExternalContext` populated from real analyses (neighbourhood, cohort, FBA, predictor statuses,
completeness, roles) or null, and a decision where the stage sits relative to enforcement passes
(enforcement after it would act on learned scores; the stage deletes pinned floors of the keys it
decides). Persist `report` (or at least model hash + per-key outcomes) with the run.

## Pending limitations

* No trained model exists; no accuracy statement is made or supported. The stage's behaviour was
  tested only with hand-set and planted-synthetic mechanism models.
* `seq.src.other` still pools unnamed tools with incomparable score scales (needs a contract
  revision to split by tool or to rank-normalise per source).
* Roles other than the derived EC → catalytic and GO complex → complex-subunit come only from
  `ExternalContext.roles`; no predictor produces them yet.
* The legacy output is role-free; consumers needing role-specific scores must read the report,
  provenance entries or annotation metadata.
* The locus/null layer needs real training groups (see learned-context.md); until then it is
  reported as unavailable.
* Context producers for neighbourhood support/conflict, cohort recurrence and conditional FBA are
  not part of this change.
