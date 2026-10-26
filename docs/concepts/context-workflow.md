# Learned, cohort and metabolic context in the production workflow

Status: **wired, not validated.** The stages below run through the shared `IntegrationWorkflow` (so
`gspa integrate`, batch runs and `AnnotationPipeline.integrate` agree), with a ledger entry for each. No trained
learned-context model exists in the repository: the learned stage refuses to run without a production artifact
and nothing here supports an accuracy statement. Every test uses explicitly synthetic data. Background:
[learned-stage.md](learned-stage.md), [learned-context.md](learned-context.md), [cohort-context.md](cohort-context.md),
[metabolic experiments](../scenarios/metabolic-experiments.md), [community cross-feeding](../scenarios/community-metabolism.md).

## Stage order and ledger

`<out>.stages.json` lists, in this order: `reference_preparation`, `cohort_context`, `evidence_integration`,
`dark_matter`, `reaction_context`, `cross_genome`, `symbolic_enforcement`, `metabolic_refinement`, `agent_completion`,
`metabolic_experiments`, `community_metabolism`, `learned_integration`, `feature_export`. Each is `executed`,
`skipped` (not requested), `failed` (with the reason) or `not_run` (an earlier stage failed). A requested stage
whose input is missing **fails the run**; it is never replaced by a quiet default.

The ledger declares `scoring_mode`: `legacy_heuristic` (the evidence combiner and hand-set priors produced the scores) or
`learned` (a frozen model replaced them). It also names the two baselines, which must not be confused:
`stage0_raw_max_score` (strongest raw predictor score per function; the cohort engine's baseline) and
`legacy_prior_baseline` (the posterior after `evidence_integration`; the ledger's `final_changes` are measured from it).

## Learned stage

| Option (`IntegrationWorkflowConfig` key) | Meaning |
| --- | --- |
| `--learned-model` (`learnedModelFile`) | Frozen production artifact. Naming it requests the stage. |
| `--learned-stage` (`learnedStage`) | Request the stage explicitly; **fails without a model**. |
| `--learned-retract-below` | Retract candidates below a probability; default: never. |
| `--learned-add-unretained-claims` | Also add claimed functions the earlier stages did not retain. |

The stage runs **last**, over the **original raw claims**, so nothing the earlier stages folded in is counted
twice; its calibrated probabilities **replace** the heuristic scores (the old posterior is audit-only). A missing
file, a mechanism-test artifact (`intended_use = mechanism_test` or synthetic) or a contract mismatch raises before
any stage runs and leaves a ledger with `reference_preparation: failed`. Weak candidates are kept and all roles are
scored; the legacy annotation holds one role's projection (see [learned-stage.md](learned-stage.md)). GO terms the
structural taxon constraints forbid for the organism are never added or rescored: the enforcement stage's decision
stands (`forbiddenGoTerms` of the stage settings, taken from the SAT checker). Pairwise disjointness among terms the
learned stage rescores is **not** re-enforced afterwards.

Artifacts: `<out>.learned.json` (model id and SHA-256, contract fingerprint, per-candidate role scores, margins,
feature digests, per-key outcomes) and, in the ledger, the model hash and counts.

## Feature export (for real training data)

`--feature-export-dir` (`featureExportDir`) writes `candidate_features.jsonl`, `feature_contract.json` and
`export_manifest.json` through the **same extractor** the learned stage scores with, with or without a model; a batch
writes `<dir>/<genome_id>/`. Rows carry **no label, split, group or example id**; the trainer refuses them until
independently curated truth and declared family/genome splits are joined by `candidate_key`. Stale exports of an
earlier run are removed at the start of a run. The contract is `gspa.context-features/2`.

## Cohort context (two-pass batch)

```bash
gspa integrate --batch-manifest batch.tsv --cohort \
  --reaction-graph seed_reactions.tsv --reaction-ec-aliases seed_ec_aliases.tsv \
  --cohort-metadata cohort_metadata.json [--cohort-orthogroups og.tsv] [--cohort-manifest refs.tsv]
```

1. **Pass 1.** The stage-0 raw baseline of **every** genome is built from raw claims only (strongest raw score per
   protein and EC; non-EC claims, partial ECs and out-of-range scores are counted and ignored). Its SHA-256 is the
   snapshot id, so genome and claim order change nothing.
2. **One joint inference** (`CohortInference`) over all genomes with their layouts, orthogroup/subfamily maps,
   clade/strain metadata and optional curated interaction edges, role evidence and reaction directions. Candidates are
   direct **EC to reaction** expansions (no EC to GO step); an EC that maps to several reactions is **one** candidate whose
   representative reaction (highest weighted support, then id) carries the features while the others are listed as
   alternatives with their own support and counted in `cohort.w.n_reaction_alternatives`; their evidence is never summed.
3. **Pass 2.** Each genome's workflow consumes only its **held-out** records (the target and its strain cluster or
   clade are removed from its peers, `--cohort-holdout`), as typed `ExternalContext.Cohort.weighted` features.

Protein ids are the FAA sequence ids that layouts, claims and orthogroups share (the same id may occur in several
genomes; candidates are namespaced by genome in memory); a mismatch fails loudly rather than joining zero rows.

**Metadata** (`--cohort-metadata`, strict JSON; see `CohortMetadata`): per genome `clade`, `strainCluster`,
`assemblyStatus` (`COMPLETE | DRAFT | UNKNOWN`), `completenessFraction`, `familySearchComplete`/`familySearchCoverage`,
`contigs` with `length` and `closed`, `plausibleReactions`, and a `tree` or `treeFile`. Unknown keys, an unknown
status, a complete assembly with a fragment contig, or `familySearchComplete` with coverage below one are errors.
Absent facts stay absent: a genome without metadata joins one **unresolved block** with one clade's weight, never an
independent peer; adjusted intervals and permutation p are then withheld. With no metadata at all every genome shares
the block, so nobody has an independent peer. **A genome with no peers gets no cohort entry** (`cohort.w.available = 0`,
`no_peers` in the report): the learned model sees "unavailable", not "recurrence negative". Absence is `TRUE_ABSENT` only
for a complete assembly with a certified all-protein family search; otherwise it is unknown. The complementary
missing-function share is an UNANNOTATED label, never evidence of loss.

A **single run** can name reference genomes in `--cohort-manifest` (`genome_id, claims_path[, layout_path[, orthogroups_path]]`);
the same genome in a batch with the same peers has identical features (tested). `<out>.cohort.json` reports the
snapshot, held-out genomes, peers, metadata gaps, orthogroup rows read, counts and every reaction ambiguity.

The weighted features keep availability: a nullable field is "unavailable", never zero. They are the `cohort.w.*`
columns of the contract (effective N, Jeffreys and raw support, clade-adjusted interval, matched background rate and
enrichment, complementarity, permutation null, copy number, the candidate's own window, reaction/family ambiguity and
raw role evidence). The engine's own re-listing of the direct EC call as catalytic role evidence is **not** passed on
(it is the same score the sequence features already carry); only supplied role evidence enters.

## Typed explicit context

`--external-context` (`externalContextFile`) is a strict JSON rendering of `ExternalContext` (neighborhood, fba, roles,
predictorStatus, assemblyCompleteness, candidates). A key that both the file and a producer (cohort, experiments)
supply for the same candidate is an error; neither silently overrides the other.

## Metabolic experiments and community analysis

See the scenario pages for the models. Production wiring (`--experiment-manifest`, `--community-manifest`):

* Both sidecars run as bounded subprocesses: a hard solver-call budget passed on the command line and re-checked in the
  result, a wall-clock limit (on expiry the process **and its children** are killed), exit status, schema, **model
  SHA-256** and "no probability, no promotion, no independent proof" checked on every record. The sidecar writes to a
  **private fresh directory**; only a result that passed every check is moved to `<out>.experiments/result.json`
  (`<out>.community/` likewise), and an earlier file that cannot be removed stops the run. Failures leave no result.
* `--experiment-mode manifest` runs a user-supplied manifest as written; the ledger says `context_source: user_manifest`
  and that this is **not** an evidence-guided workflow. Only hypotheses that name `function_key` can become features.
* `--experiment-mode evidence_guided` treats the manifest as a **template** (medium, biomass, explicit maintenance are
  its assumptions, copied and recorded). A candidate is an EC claim of this run, **retained or not** (weak claims the
  integration dropped are included). It becomes a **reaction-level** hypothesis only through a validated EC to
  model-reaction mapping (the refiner's `reactions.tsv`, the model's own EC annotation, or a reaction-graph id that is a
  model reaction id). The gene is not required: a protein missing from the model's gene rules, the typical missing
  annotation, is tested as a reaction with `gene_in_model: false` and never as a gene knockout. A gene label without an EC
  mapping makes the candidate unavailable. Above `--experiment-max-hypotheses` candidates are ranked from measurements
  (more distinct evidence types, then probability closest to 0.5, exact ties by key; no fitted weights, not a
  probability) and the ranking and cut-off are in the ledger. Evidence: what the run measured (raw claim, context-positive
  window, a supporting cohort peer), claim metadata naming an **exact model exchange id** (`substrate_exchange` with
  `substrate_role`, `product_exchange`; others are listed as unmapped), plus explicit `--experiment-evidence` records.
  Nothing is inferred from metabolite names. With nothing admissible the ledger says `no_admissible_experiment` and no
  solver starts.
* Results map to `ExternalContext.FbaConditional` from the **per-reaction checks** of each condition: share with
  **required** flux, share with **optional** flux, and the reaction knockout's growth relative to the same condition (only
  when one reaction was assessed). A whole-gene knockout is compared with the baseline medium only: it is kept in the ledger
  (`own_knockouts_not_used_as_features`) and is never a flux class, a deletion effect of a different medium, or proof that each
  reaction of a multi-reaction gene is required; a candidate with only a gene knockout has no conditional feature. A product
  objective, or a reaction the refiner added to make growth feasible (`flux-certified repair`), cannot confirm the candidate:
  those conditions are flagged and excluded, and if nothing else remains the candidate gets `objectiveRequiresCandidate =
  true`. Skipped, infeasible and inconclusive runs count in neither direction.
* Community results stay **organism-specific**. Cross-feeding that restores growth explains a missing route without
  assigning an enzyme or gene and is not used as evidence for any candidate; no community feature is defined, so none is
  fabricated and nothing is read as negative evidence. The raw conditional records are exported verbatim. A batch runs
  the community once and every genome ledger points at the same result.
* `qwen` as `--experiment-strategy` is forwarded only as `--strategy qwen --qwen-model <agentModel> --api-base <agentApiBase>`
  (model `qwen3.8-27b`, endpoint `http://unimatrix01.kaust.edu.sa:8000/v1`); the sidecar reads its own credentials,
  none is read, passed or logged here, and the agent can only choose a validated menu id. The tests never contact it.

## Routes

| Route | Learned | Cohort | Features | Experiments / community |
| --- | --- | --- | --- | --- |
| `gspa integrate` (single, batch) | yes | single run with references, or two-pass batch | yes | yes (community once per batch) |
| `AnnotationPipeline.integrate` / `gspa annotate --integrate` | yes (`integration.workflow` block) | single run with references | yes | yes |
| Nextflow `INTEGRATE` | `--learned_model` | **not run** (batch-level job) | `--feature_export` | **not run** |
| Web service | not wired | not wired | not wired | not wired |

## Remaining gaps

* No trained production model; the root owner trains it from exported rows and benchmarks it. Until then the learned
  stage cannot run on real data and `scoring_mode` stays `legacy_heuristic`.
* The Stage 2B refiner still updates posteriors before the learned stage; a learned run replaces them from raw
  evidence, so the refiner reaches the model only through the conditional experiment features.
* Nextflow has no batch cohort process and stages no experiment files; the web service wiring is separate work.
* No neighbourhood producer fills `ExternalContext.Neighborhood`; the cohort engine's window features
  (`cohort.w.local_*`) are the neighbourhood information the model receives, and `ctx.nbr.*` stays unavailable unless an
  explicit context file supplies it.
* Per-genome plausible reaction sets exist only if the metadata file lists them.
* Transport and BGC evidence for experiments must be supplied; the planner uses cohort evidence only to order its menu.
* Community features for candidates are undefined; the community result is exported, not consumed.
