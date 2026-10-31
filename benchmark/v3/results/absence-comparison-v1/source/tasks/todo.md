# Todo: GSPA improvements

## Audit follow-up — 2026-10-26 (completion reopened)

Robert selected all eleven audit recommendations for implementation and requested
the agent-guided metabolic experiment extension below. The follow-up audit contradicted the blanket completion claim. All twelve original
requirements remain the acceptance scope; see [the current acceptance record](v3-acceptance.md).
The parenthetical claims below are historical implementation reports, not accepted
validation. Audit baseline: commit `7714377`; reviewed release: `498560f`.

Latest progress: [explicit absence and supply comparisons](../docs/benchmarks/v3-absence-comparison.md)
now compare four reaction-level scenarios through the production workflow, with
retained-candidate binding, bounded supply, transport interventions, unknown states
and complete-condition budgets. Nine CLI runs preserve predictions and features;
147 independent GLPK comparisons agree. Gates pass 1,251 JVM cases (two skips),
613 benchmark cases (three skips, ten live Nextflow cases) and 89 metabolic Python
cases. The full 18-genome pathway replay also passes all 72 CLI calls, preserving
56,613 raw candidates and automatic/explicit gap parity across 484,202 homology
candidates. Biological absence/gene-attribution validation, automatic cohort-wide
scenario construction and whole-item acceptance remain open.

Previous progress: [three context modes](../docs/benchmarks/v3-context-modes.md) now
run through the workflow, CLI and Nextflow, with mode-aware model abstention.
The gates pass 1,245 JVM tests (two skips) and 610 benchmark tests (three skips,
ten live Nextflow cases). The complete 18-genome search acquisition verifies.
Four matched export arms each preserve 484,202 candidates and functional TSV
bytes; default mode also reproduces every previous numerical feature. Six models
preserve all 4,262 supervision records; 25,572 JVM/Python scoring comparisons agree.
All paired genome-level Brier intervals cross zero, and null learning remains
unsupported. All 1,281 archived payloads are independently rehashed. Reference independence,
biological validation and whole-item acceptance remain open.

Previous progress: [reaction-signature availability](../docs/benchmarks/v3-signature-availability.md)
now keeps absent or insufficient expected signatures unknown rather than negative.
Eight added JVM cases and ten CLI runs verify the repair; the full gates pass
1,229 JVM cases and 600 benchmark cases, including nine live Nextflow checks.
Matched old/new eighteen-genome exports preserve all 484,202 identities and
functional TSV bytes, correcting 4,780 candidates' context features. Both models
are refitted with all 4,262 supervision records/splits unchanged. All 118 whole-genome
duplicate/shuffle control outputs verify; duplicate probabilities remain exactly
invariant for both models, while calibration transfer still fails. Null learning and biological acceptance
remain open. No whole item is accepted.

Earlier [eight-genome production replay and reference audit](../docs/benchmarks/v3-cohort-pathway-replay.md)
passes all 32 CLI calls, preserving 42,795 raw candidates and matching automatic
and explicit gap features/suggestions across 264,687 homology candidates.
The benchmark gate passes 600 tests, twelve optional skips. Full-sequence checks
confirm 395 target proteins present identically in named references, including
93 test proteins; unseen-sequence performance is not established. Repeated
searches also differ in two weak-hit candidates despite matching input/reference
hashes, with unchanged pathway/gap outputs and an unresolved cause. All findings
and outputs are archived. No whole item is accepted.

Earlier [full-cohort pathway acquisition](../docs/benchmarks/v3-cohort-pathway-acquisition.md)
targets all 18 fixed proteomes (49,463 proteins). All eight bacterial searches
pass the independent audit, with 147,818 reaction rows and 15,384 pathway rows. Real execution
exposed incomplete archaeal references and a directory-marker requirement;
failures are retained and a separate checksum-verified reference acquisition is
running before dependent retry array 10131. Eighteen added checks pass in the 591-test benchmark
gate. Search completion, full-cohort production gap replay and matched context evaluation
remain open; current scheduler and receipt state take precedence over this checkpoint.

Earlier [independent RLC anchor mapping](../docs/benchmarks/v3-rlc-mapping.md)
now honors standalone EC-to-GO input and retains alternatives without counting
correlated labels as independent diversity. Ten added JVM cases and 29 CLI
replays pass. On the real training input, standalone mapping recovers 820 feature
rows and 41 suggestions, exactly matching the earlier pathway-backed route.
Full gates pass 1,221 JVM tests, 573 benchmark tests and nine live Nextflow checks.
These are engineering results; RLC calibration and biological validation remain open.

Earlier [observed pathway-gap integration](../docs/benchmarks/v3-pathway-gap-bridge.md)
now retains both source tables and automatically derives the 50 audited gap
hypotheses through shared integration. Ten real CLI replays verify annotation,
observation-file and explicit-gap feature parity, preserving raw claims and
baseline predictions. Complete/zero/missing/failed cases stay distinct. Full
gates pass 1,211 JVM tests, 573 benchmark tests, and nine live Nextflow checks.
Calibrated context and biological validation remain open; the subsequent RLC
mapping repair is recorded above.

Earlier [local reaction context without peers](../docs/benchmarks/v3-local-reaction-context.md)
now reaches the shared learned/export path with recurrence explicitly unavailable.
Eight real CLI runs on 29,079 unchanged claims and 50 observed gap hypotheses
verify 1,182 restored local contexts, identical sequence evidence, and exact
feature compatibility when independent peers exist. Eleven added JVM and three
Python regressions pass; full gates are 1,172 JVM tests, 572 benchmark tests, and
nine live Nextflow checks. Current models abstain on the new local contexts.
Automatic gap wiring and calibrated/held-out context validation remain open.

Earlier the [predictor-role handoff](../docs/benchmarks/v3-claim-role-bridge.md)
now preserves all 829 observed complex-subunit candidates through learned
candidate construction. Fourteen regressions pass; the old production CLI fails
the new oracle and the repaired engine passes seven CLI replays, the 1,161-test
JVM gate (two skips), and nine live Nextflow checks. Raw claims and functional
TSVs remain byte-identical; both models still abstain on the unsupported source.
Independent role truth and trained role/locus/null inference remain open.

Earlier [genuine gapsmith output and adapter repair](../docs/benchmarks/v3-gapsmith-observed.md)
now preserves all 13,147 observed reaction rows and 3,469 candidates through
production annotation/integration. Input, configuration and score-type defects
are repaired; weak hits and subunit evidence are retained. The 1,147-test JVM
gate, nine live Nextflow checks and seven CLI replays pass. Both learned models
abstain on this new source. A complete-table gap audit yields 50 candidate
records across 46 reactions. Context wiring and biological validation remain open.

Earlier [reaction-local missing evidence](../docs/benchmarks/v3-rlc-missing-evidence.md)
now contributes no artificial probability-0.5 support or penalty. Three new
regressions and ten before/after CLI runs verify the fix; 1,113 JVM tests (two
optional skips) and nine live Nextflow checks pass. The [gap-table converter](../docs/benchmarks/v3-gap-table-contract.md) also
passes 17 CLI regressions after repairing reaction-level eligibility; the broader
569-test benchmark gate passes (12 optional skips). The full Price collection
audit also verifies all 851 bound files across 132 parts. Optimized dense-prior
execution and its independent comparison pass: all 9,790,188 identities and
complete prediction bytes match. No whole item is accepted.

Earlier [heterogeneous-strain evidence identities](../docs/benchmarks/v3-cohort-evidence-identity.md)
repair the reproduced within-cluster duplication defect. Twenty-seven added JVM
cases and the archived production replay pass. Both models have been refitted
with all supervision preserved; all 118 updated production control outputs verify.
Both models now have exactly invariant numeric probabilities across all 18 original
genomes under the real exact-copy intervention. The
1,082-test JVM and 551-test benchmark gates pass (six/three optional skips).
All 132 full Price matrix parts are published and count-verified; a separate
collection hash audit has now passed. No whole item is accepted.

Earlier [InterProScan statistic contracts](../docs/benchmarks/v3-interpro-observed.md)
now distinguish E-values from profile scores and reject malformed statistics.
An actual installed-tool run and seven historical-table parser replays preserve
annotation identity and reproduce/fix negative scores. JVM and live Nextflow
checks pass; the full actual-tool annotation/direct-integration replay also passes. Method-specific
calibration and the other predictor acceptance requirements remain open.

Earlier [supported complex and multifunction assignments](../docs/benchmarks/v3-masked-assignment-experiments.md)
now construct explicit partner AND clauses and evidence-supported joint activities,
with alternative mappings, role restrictions and bounded enumeration. Sixteen new
JVM cases and ten actual CLI runs validate construction and hidden-assignment
invariance. Independent checks agree on 200 new and 54 earlier-fixture phenotype/FVA
comparisons. Two cumulative-knockout solver bugs and an order-dependent independent
oracle defect were reproduced, repaired and regression-tested. These are synthetic
mechanism checks; biological recovery and the other item-12 gates remain open.

Earlier [essentiality indexing and protein identity](../docs/benchmarks/v3-essentiality-index.md)
repair a real pipe-ID rescue crash and remove repeated full-candidate scans.
All 240 synthetic prior states match before/after; the indexed build passes
1,010 JVM tests and 525 benchmark tests. Complete EnzBERT capacity and independent
comparison retain all 349,400 candidates with byte-identical predictions;
dense-prior checks remain running/queued. Earlier [real HMMER execution](../docs/benchmarks/v3-hmmer-observed.md)
exposed and now verifies repairs to availability, malformed statistics and CLI
thread forwarding, including the separate AntiFam adapter. The [AntiFam quality repair](../docs/benchmarks/v3-antifam-quality.md) now preserves
typed raw observations and region metadata through integration and projection.
Current learned models abstain on flagged proteins; a calibrated effect on
functional scores remains open. Annotation and direct integration agree on the observed
synthetic fixture; other adapters' genuine-output coverage remains open.
Verified-batch timeouts now terminate the inference process group, with a
failing-before/passing-after regression archived in `benchmark/v3/results/batch-timeout-v1/`.
The [100-candidate scoring discontinuity](../docs/benchmarks/v3-candidate-count.md)
is fixed and verified through the production CLI; JVM and live Nextflow checks pass.
The [lossless full-matrix export](../docs/benchmarks/v3-price-streaming-export.md)
now independently checks every claim and reproduces the compact artifact's
float32/obsolete-EC statistics without changing raw exports. The complete 132-matrix
run is ongoing; two repaired real cases and complete TopEC/EnzBERT production
capacity checks pass. The 9.79-million-candidate DeepProZyme baseline also passes
exact retention and post-run hash verification; dense priors/cohort measurements
remain pending. [Model input coverage](../docs/benchmarks/v3-model-input-support.md) and the
[GOA source audit](../docs/benchmarks/v3-gaf-negative-audit.md) remain earlier verified
progress. All twelve whole-item acceptance gates remain open.

- [ ] **1. Benchmark the actual application.** Make benchmarks invoke the shared
  production inference engine; generate public claims from versioned experiment
  records. Compare matched sequence evidence with each context stage and the full
  system on held-out genomes/families, with independently selected thresholds,
  calibration, correction/harm counts, and genome-level uncertainty.
  *(Prior implementation report, pending acceptance: `benchmark/v3/run_v3_benchmarks.py`, `benchmark/v3/evaluate_training_cohort_ablations.py`, `benchmark/v3/run_price_context_ablation.py`; 73 runner tests + 429 v3 benchmark tests passing; reports in `docs/benchmarks/v3-*.md`.)*
- [ ] **2. Unify execution routes and prediction state.** Align CLI annotation,
  integration, web, batch, and Nextflow; preserve contributions from combined
  suggesters and synchronize scores, annotations, provenance, and downstream state.
  Require matched-input parity and explicit executed/skipped/failed stage records.
  *(Prior implementation report, pending acceptance: `gspa.integration.workflow.IntegrationWorkflow` shared across `IntegrateCommand`, `AnnotationPipeline.integrate`, web portal, and Nextflow `INTEGRATE`, with schema-2 `<genome>_integrated.tsv.stages.json` stage ledgers and synchronized `Prediction` identity.)*
- [ ] **3. Repair predictor contracts and failure reporting.** Correct eggNOG score
  parsing, HMMER database wiring, genomic-tool dispatch, InterPro EC output, and
  operon GO-aspect handling. Use typed contracts and structured statuses; verify
  genuine output fixtures, complete adapter paths, and injected failures.
  *(Prior implementation report, pending acceptance: `gspa-predictors` and `benchmark/neural/` sidecars updated with calibrated bitscore/e-value parsing, typed adapter statuses `OK`/`SKIPPED`/`FAILED`/`PARTIAL`, genuine output fixtures, and biological failure-injection tests.)*
- [ ] **4. Learn evidence weights and final confidence.** Replace arbitrary boosts,
  candidate-count discontinuities, and irreversible confidence floors with learned,
  calibrated, revisable inference. Retain a no-candidate outcome and account for
  correlated evidence. Compare interpretable and tree-based models with existing
  rules using out-of-fold inputs and family/genus-held-out evaluation.
  *(Prior implementation report, pending acceptance: `LearnedIntegrationStage`, `LearnedContextModel`, `scripts/train_context_model.py`, contract `gspa.context-features/2` with 114 features and explicit `<null>` no-candidate outcome; trained on out-of-fold family-disjoint splits of 18-genome `swissprot_label_rich_v1` (`4,262` rows); held-out test `learned_logistic` AUROC `0.7934`, improving GO $F_{\max}$ from `0.2881` to `0.3221` (`+0.0340`, precision `0.2733 -> 0.3529`); `docs/benchmarks/v3-training-cohort.md`.)*
- [ ] **5. Connect neighborhoods to reaction and interaction networks.** Introduce
  shared protein–reaction candidates with explicit context features, alternative
  mappings, reaction direction, organism-specific feasibility, and distinct
  catalytic/interaction roles. Compare pathway-only, universal-network, and
  organism-specific context on held-out enzyme loci.
  *(Prior implementation report, pending acceptance: `EcCandidateExpander`, `ReactionLocalContextSuggester`, and `ContextFeatureExtractor` bridge reaction-graph distance, directionality, alternative mappings, and roles; evaluated in the 10-arm `swissprot_label_rich_v1` ablation and 42-run Price-148 benchmark; `docs/benchmarks/v3-training-cohort.md`, `docs/benchmarks/v3-price-development.md`.)*
- [ ] **6. Infer recurrence jointly across genomes.** Integrate cohort-wide
  ortholog/subfamily context inference into batch workflows; correct singleton
  normalization, use matched-window backgrounds, account for ancestry and assembly
  completeness, and combine complementary missing-function patterns. Test clade
  holdouts, duplicated strains, shuffled neighborhoods, and fragmented assemblies.
  *(Version-3 replay: 118 production genome/arm outputs across logistic/tree models. Exact-copy probabilities are invariant within 1e-12 over this 18-genome cohort after raw-count removal and permutation repair. A broader heterogeneous-strain-cluster failure is now repaired and regression-tested through complete observed-input identities; the refitted models’ 118 production control outputs verify exact probability invariance across all 18 original genomes. Complementary-gap inference and independent biological acceptance remain open. See `docs/benchmarks/v3-cohort-contract.md`.)*
- [ ] **7. Preserve weak candidates and resolve roles within neighborhoods.**
  Separate candidate retention from display thresholds; distinguish enzymes,
  accessory proteins, transporters, regulators, paralogs, and required complex
  subunits. Evaluate candidate recall and locus ranking with real predictor scores
  and difficult same-operon alternatives.
  *(Prior implementation report, pending acceptance: lossless `--claims-out` and weak-candidate retention down to `1e-9`, `CandidateRole` typing (`CATALYTIC`, `ACCESSORY`, `TRANSPORTER`, `REGULATOR`, `COMPLEX_SUBUNIT`), and same-operon / nearby-locus decoy evaluation on Price-148 across 115,542 proteins; `docs/benchmarks/v3-price-development.md`.)*
- [ ] **8. Model absence, missing sequence, and environmental dependence.** Treat
  missing annotation, incomplete assembly, genuine pathway loss, and external
  supply as competing explanations. Expose medium, biomass, transport, and model
  assumptions. Validate gene attribution separately from growth/knockout prediction
  on complete genomes, auxotrophs, and controlled fragmented assemblies.
  *(Prior implementation report, pending acceptance: explicit stage-ledger FBA medium/biomass/Gram assumptions, competing absence states (`true_absent`, `unobserved_absent`, `pathway_without_family`), 61-case 3-protocol *E. coli* auxotroph validation with independent COBRApy/GLPK cross-verification, and 50% assembly fragmentation control; `docs/benchmarks/v3-metabolic-development.md`, `docs/benchmarks/v3-training-cohort.md`.)*
- [ ] **9. Implement community cross-feeding.** Replace the disconnected heuristic
  scaffold with separate organism models, explicit transport, extracellular
  exchange, and environmental constraints for a defined community. Validate with
  measured coculture/monoculture phenotypes, exchanged metabolites, and transporter
  knockouts; measure reductions in false enzyme assignments separately.
  *(Prior implementation report, pending acceptance: coupled multi-organism SBML FBA in `scripts/community_metabolism.py` with matched no-transfer controls, transport knockouts, and independent COBRApy oracle verification across all 1,035 Wintermute 2010 auxotroph pairs × 6 environmental regimes = 6,210 coupled LPs; `docs/benchmarks/v3-community-development.md`.)*
- [ ] **10. Report prediction changes and uncertainty.** Show stage-wise changes,
  baseline/final scores, evidence, alternatives, and module status; use ontology
  ancestry and EC identifiers for enzyme classification. Reconcile browser counts
  with exports and distinguish internal audit metrics from experimental accuracy.
  *(Prior implementation report, pending acceptance: `scripts/visualize_results.py` and `service/` stage-wise change ledger, GO ancestry and EC hierarchy classification, learned role/null-hypothesis (`<null>`) panels, and headless Chromium browser reconciliation tests.)*
- [ ] **11. Evaluate the agent against equally informed simple methods.** Reserve
  agent use for unresolved cases and useful additional evidence gathering. Compare
  with feature rankings and learned classifiers using identical evidence on unseen
  families/genomes; report correct and harmful decisions, abstentions, and cost.
  *(Prior implementation report, pending acceptance: `benchmark/v3/evaluate_stage3_agent.py` comparing `qwen3.8-27b` (`http://unimatrix01.kaust.edu.sa:8000/v1`, guarded and unguarded) against 6 deterministic and learned comparators on 40 held-out test anomalies of `swissprot_label_rich_v1`; `learned_gbdt` (`+23` net gain, `0.0001 s`) and `learned_logistic` (`+22` net gain, `0.0001 s`) outperform `qwen3.8-27b` (`+15` unguarded / `+7` guarded, `100.0%` precision, `0` harmful, `131.6 s`); `docs/benchmarks/v3-stage3-agent-evaluation.md`.)*
- [ ] **12. Let the agent design and run informative metabolic experiments.**
  Implement the proposal below, building on items 2, 4–8, and 11; use item 9 when
  proposing community exchange experiments.
  *(Prior implementation report, pending acceptance: `scripts/metabolic_experiments.py`, `benchmark/v3/run_planner_benchmark.py`, and `benchmark/v3/audit_condition_dependent_production.py` evaluating `fixed_paired`, `random_item`, `random_paired`, `evidence_item`, `evidence_paired`, and `qwen` (`qwen3.8-27b`) under equal solver budgets $B \in \{8, 16\}$ across 44 genes × 24 substrates: at $B=8$, strict Qwen-attributable discovery is `85.4%` vs `39.0%` for `evidence_paired` (McNemar $p=0.0000$); at $B=16$, `80.5%` strict / `87.8%` hybrid vs `75.6%` for `evidence_paired`, with `75.9%` (`22/29`) binary agreement on visited measured pairs; `docs/benchmarks/v3-planner-evaluation.md`.)*

### 12. Agent-guided metabolic experiments

A single medium and growth objective can leave alternative catabolic pathways and
secondary metabolism unused, so their absence from an FBA solution is weak evidence
against a candidate function. The agent should use sequence/domain evidence,
neighborhoods, cross-genome recurrence, transport predictions, biosynthetic gene
clusters, and existing gaps to shortlist competing functional explanations and the
conditions most likely to distinguish them. It should then design and execute
bounded computational experiments: change carbon/nitrogen sources, oxygen or
nutrient availability, remove supplements, test knockouts, and, where supported,
explore product formation alongside growth or maintenance requirements. Compare
alternative feasible flux solutions and candidate assignments, retain negative or
inconclusive results, and record every assumption. For secondary metabolism, the
relevant reactions must be represented; feasible production is a conditional
hypothesis, not evidence that a pathway is expressed. Choosing an objective that
requires a candidate function must not itself count as independent confirmation.

The expected improvement is recovery of condition-specific functions that a single
growth simulation misses, with fewer unnecessary simulations because prior
evidence guides the search. Benchmark on held-out organisms with measured growth
across media, knockout phenotypes, and characterized secondary-metabolite
production, masking target gene assignments from all model and evidence inputs.
Compare single-medium FBA, a fixed condition panel, random condition selection, a
simple evidence-guided planner, and the agent under the same solver-call budget.
Measure correct gene/function recovery at fixed precision, harmful assignments,
agreement with withheld phenotypes, abstention, and computational cost. Freeze the
planner and confidence rules before testing; count cases where no function can be
recovered rather than evaluating only successful repairs.

Methodological reference: [COBRA Toolbox v3.0: alternative optimal flux solutions
and flux variability analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC6635304/).
Existing starting points: `scripts/stoichiometric_locus_refiner.py` and
`service/metabolic_conditions.py`; current condition probes do not implement this
adaptive experiment-selection loop.

## Historical checklist: GSPA v1.5.0 release

Companion checklist for `tasks/plan.md`. Strike tasks as they complete.

The previous Phase 11 (FM-based operon) work is parked on the
`parking/phase11-glm` branch (NO-GO verdict, 2026-05-05).

## Phase A — Triage and yank

- [x] **A1** `./gradlew clean test` exits 0 (verified 2026-05-06: 335 tests, 0 failures post-yank)
- [x] **A2** Park Phase 11 to `parking/phase11-glm`; remove gLM files; remove `--operon-caller` switch
- [x] **A3** Fix CLI version drift via `VersionProvider` (build.gradle.kts → 1.5.0-SNAPSHOT, picocli reads from generated `version.properties`)
- [x] **A4** Add GPL-3.0-or-later `LICENSE` at repo root; update `README.md` license section

## Phase B — gspa-nf integrator parity

- [x] **B1** New `gspa-nf/modules/integrate.nf` with `BUILD_CLAIMS` + `INTEGRATE` processes
- [x] **B2** Wire into `gspa-nf/main.nf` workflow (guarded by `params.run_integrate`)
- [x] **B3** Add params to `gspa-nf/nextflow.config`
- [ ] **B4** (deferred) Add GO OWL / ec2go / pathways placeholders to `gspa-nf/databases.config` — params already added to `nextflow.config`; databases.config update can ride with the first real cluster run
- [x] **B5** (logical smoke) New jar runs `gspa integrate` on bench9 mgenitalium claims → 4,228 integrated annotations. INTEGRATE process body validated. Full Nextflow Docker smoke deferred (lint pass via `nextflow inspect` already confirmed wiring; dev host has no Docker daemon)
- [x] **B6** Document new flag in `gspa-nf/README.md`

## Phase C — Phase 10 retune (parallel to B)

- [x] **C1** Real gapseq output **already on cluster** (4 of 10 genomes have `*_real.jsonl`; remaining 6 fall back to synthetic 400-gap `*_gaps.jsonl` per the documented zero-byte Reactions.tbl mitigation)
- [~] **C2** Running on unimatrix01 SLURM job 4069 — `phase10_retune.sh` integrates each genome at qBase ∈ {0.50, 0.70, 0.75} plus C1 baseline = 40 integrations
- [~] **C3** Same job scores each (config, tag) with `benchmark_pgap_v2.py` (200 bootstrap, micro + CAFA)
- [ ] **C4** Decide Phase 10 default (default-on iff Δ ≥ 0 with no genome regressing > 0.01)
- [ ] **C5** Update `benchmark/RESULTS.md` Phase 10 section with the verdict

## Phase D — mdF comparison

- [~] **D1** Installing `mDeepFRI` from PyPI into `/data/hohndor/envs/mdf-venv` on unimatrix01 (first attempt stuck on metadata for 13 min, killed; retrying with `--no-cache-dir`)
- [x] **D2** 13-genome PGAP panel: bench9 (9 genomes) + bench10 (10 PGAP genomes) minus duplicates = 13 unique. FAA files confirmed under `/data/hohndor/gspa/proteomes/{tag}.faa` (bench9) and `/data/hohndor/gspa/proteomes/bench10/{tag}.faa` (bench10)
- [ ] **D3** Run mdF on the 13 panel genomes (gated on D1)
- [ ] **D4** Adapter `benchmark/parse_mdf_predictions.py` (with `--self-test`)
- [ ] **D5** Compute F-max micro + CAFA with bootstrap CIs
- [ ] **D6** (Optional) Ensemble: GSPA + mdF substituting for DeepFRI
- [ ] **D7** Add "v1.5.0 — comparison with metagenomic-deepFRI" section to `RESULTS.md`
- [ ] **D8** Cite Bezshapkin et al. 2026 in `README.md`

## Phase E — CHANGELOG and CI

- [x] **E1** `CHANGELOG.md` compiled (Keep-a-Changelog, [Unreleased] entry + v1.0.0 → v1.4.1 history)
- [x] **E2** `.github/workflows/test.yml` (JDK 21, Gradle cache, version smoke)
- [ ] **E3** (Optional, deferred) `.github/workflows/release.yml`
- [ ] **E4** (Optional, deferred) `.github/workflows/docker.yml`

## Phase F — Merge and tag

- [ ] **F1** Pre-merge check: clean test green, shadowJar builds, smoke test passes, RESULTS.md updated
- [ ] **F2** Bump `1.5.0-SNAPSHOT` → `1.5.0`
- [ ] **F3** Open PR `phase11-crossgenome → main`
- [ ] **F4** Tag `v1.5.0` on `main` (annotated)
- [ ] **F5** Branch hygiene (delete `phase11-crossgenome`; preserve `parking/phase11-glm`)
