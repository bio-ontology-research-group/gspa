# Todo: GSPA improvements

## Audit follow-up — 2026-10-26

Robert selected all eleven audit recommendations for implementation and requested
the agent-guided metabolic experiment extension below. All twelve items are
pending; this checklist records scope and validation, not completed implementation.
Audit baseline: commit `7714377`. Preserve independent evaluation of annotation
accuracy; greater model feasibility or annotation completeness alone is not proof
of a correct protein function.

**Production wiring status (branch `agent/v3-context-workflow`).** Implementation of items 2, 4–9 and 12 now has real
production routes in the shared workflow (see `docs/concepts/context-workflow.md`): learned stage (last, over raw
claims, replaces heuristic scores, fails without a production artifact), feature export for real training rows, a
two-pass cohort batch, typed weighted cohort features, bounded metabolic-experiment and community execution with
conditional features. Nothing is checked off: there is no trained model, no empirical benchmark and no accuracy claim
(all tests are synthetic). Open before any release: root trains the model from exported rows and runs the benchmarks;
Nextflow cohort/experiment processes and web-service wiring are not done; transport/BGC evidence for the planner must be
supplied; no community feature contract exists.

- [ ] **1. Benchmark the actual application.** Make benchmarks invoke the shared
  production inference engine; generate public claims from versioned experiment
  records. Compare matched sequence evidence with each context stage and the full
  system on held-out genomes/families, with independently selected thresholds,
  calibration, correction/harm counts, and genome-level uncertainty.
- [ ] **2. Unify execution routes and prediction state.** Align CLI annotation,
  integration, web, batch, and Nextflow; preserve contributions from combined
  suggesters and synchronize scores, annotations, provenance, and downstream state.
  Require matched-input parity and explicit executed/skipped/failed stage records.
- [ ] **3. Repair predictor contracts and failure reporting.** Correct eggNOG score
  parsing, HMMER database wiring, genomic-tool dispatch, InterPro EC output, and
  operon GO-aspect handling. Use typed contracts and structured statuses; verify
  genuine output fixtures, complete adapter paths, and injected failures.
- [ ] **4. Learn evidence weights and final confidence.** Replace arbitrary boosts,
  candidate-count discontinuities, and irreversible confidence floors with learned,
  calibrated, revisable inference. Retain a no-candidate outcome and account for
  correlated evidence. Compare interpretable and tree-based models with existing
  rules using out-of-fold inputs and family/genus-held-out evaluation.
- [ ] **5. Connect neighborhoods to reaction and interaction networks.** Introduce
  shared protein–reaction candidates with explicit context features, alternative
  mappings, reaction direction, organism-specific feasibility, and distinct
  catalytic/interaction roles. Compare pathway-only, universal-network, and
  organism-specific context on held-out enzyme loci.
- [ ] **6. Infer recurrence jointly across genomes.** Integrate cohort-wide
  ortholog/subfamily context inference into batch workflows; correct singleton
  normalization, use matched-window backgrounds, account for ancestry and assembly
  completeness, and combine complementary missing-function patterns. Test clade
  holdouts, duplicated strains, shuffled neighborhoods, and fragmented assemblies.
- [ ] **7. Preserve weak candidates and resolve roles within neighborhoods.**
  Separate candidate retention from display thresholds; distinguish enzymes,
  accessory proteins, transporters, regulators, paralogs, and required complex
  subunits. Evaluate candidate recall and locus ranking with real predictor scores
  and difficult same-operon alternatives.
- [ ] **8. Model absence, missing sequence, and environmental dependence.** Treat
  missing annotation, incomplete assembly, genuine pathway loss, and external
  supply as competing explanations. Expose medium, biomass, transport, and model
  assumptions. Validate gene attribution separately from growth/knockout prediction
  on complete genomes, auxotrophs, and controlled fragmented assemblies.
- [ ] **9. Implement community cross-feeding.** Replace the disconnected heuristic
  scaffold with separate organism models, explicit transport, extracellular
  exchange, and environmental constraints for a defined community. Validate with
  measured coculture/monoculture phenotypes, exchanged metabolites, and transporter
  knockouts; measure reductions in false enzyme assignments separately.
- [ ] **10. Report prediction changes and uncertainty.** Show stage-wise changes,
  baseline/final scores, evidence, alternatives, and module status; use ontology
  ancestry and EC identifiers for enzyme classification. Reconcile browser counts
  with exports and distinguish internal audit metrics from experimental accuracy.
- [ ] **11. Evaluate the agent against equally informed simple methods.** Reserve
  agent use for unresolved cases and useful additional evidence gathering. Compare
  with feature rankings and learned classifiers using identical evidence on unseen
  families/genomes; report correct and harmful decisions, abstentions, and cost.
- [ ] **12. Let the agent design and run informative metabolic experiments.**
  Implement the proposal below, building on items 2, 4–8, and 11; use item 9 when
  proposing community exchange experiments.

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
