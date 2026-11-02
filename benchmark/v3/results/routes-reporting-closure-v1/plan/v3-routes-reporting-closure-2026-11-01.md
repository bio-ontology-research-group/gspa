# Pre-registered plan: closing acceptance points 2 (routes and state) and 10 (reporting)

Date: 2026-11-01. Author: Fable agent for Robert Hoehndorf. Frozen before any check is run; the
SHA-256 of this file is recorded in the closure report and in the evidence folder.

Engine under test: `gspa-cli/build/libs/gspa-3.0.0.jar`, SHA-256
`8391d73fc1f6826492918a6116d40de4ba8147607ab406a895c29c32e067e80d` (release candidate 3, staged as
`/project/k10018/gspa-v3/release-candidate-3/`). Code changes made for this closure are limited to the
web service (`service/`) and are behind options whose defaults leave today's behaviour unchanged; if any
JVM source changes, the JVM gate is rerun and the JAR rebuilt and re-hashed before the checks.

The bars are copied from `tasks/v3-acceptance.md` and are not changed after seeing results.

## Point 2: routes and state

Bar: "Single CLI, annotation, web, batch and Nextflow share scores, provenance and candidate identity,
including joint cohort execution. Requested/missing/no-work/failed states tested."
Empirical: "Matched-input comparisons through each actual route, including multi-genome context and
mixed suggesters."

### Checks that count as meeting the bar

R1. Web route wiring. `service/integration_options.workflow_config` can pass a learned model
(`learnedModelFile`), feature export (`featureExportDir`) and the metabolic-experiment stage
(`experimentManifestFile`, `experimentSidecar`, `experimentPythonBin`, mode/strategy) into the same
typed workflow JSON that `gspa integrate --workflow-config` reads. All three are off unless the deployment
sets them (environment variables); with nothing set the JSON written for a job is byte-identical to the
JSON written before this change (regression test on a fixed job description). Keys must be in
`integration-workflow-config.keys.json`. Pass = tests pass and the unchanged-default test passes.
Cohort: a web job is one uploaded genome; joint cohort context needs supplied orthogroups for the
uploaded proteins against a cohort, which the service does not compute. If a multi-genome web job
is not feasible, that is documented as batch-only (CLI `--batch-manifest`, Nextflow) with the reason;
this does not by itself make point 2 "met" or "not met": the verdict depends on R2 to R4.

R2. Single-genome matched comparison on one real genome (GCA_000011125.1, 1,700 proteins, from the
training cohort, taxon NCBITaxon_272557). Mixed suggesters: DIAMOND against the cohort's reviewed
reference (`training-cohort/baseline/reference.dmnd`) and DeepGO-PlusPlus-Light (sidecar, CPU) run
once inside `gspa annotate --proteins`, with `--claims-out`. Routes:
  a. `gspa annotate --proteins ... --integration-config cfg.json` (annotation route, in-process);
  b. `gspa integrate --claims <exported claims> --workflow-config cfg.json` (single CLI route);
  c. the web portal worker (`service.jobs.worker` running `service.analysis.run_full` on the same
     protein file, with analyses other than homology/names disabled, assets at
     `/home/leechuck/dgw_test/dgpp_assets`, the same reference DIAMOND database, the same
     ontology/EC/KEGG references, the learned model enabled through the new option);
  d. Nextflow: the single-genome case is not run through `main.nf` because its predictors run in
     containers that are not available here; the Nextflow route is tested in R3 on supplied claims,
     which is the route the production pipeline uses after prediction.
Pass = (i) the exported claims of (a) and the `claims.raw.jsonl` of (c) are identical as sets of
(protein, type, function, source, raw score, calibrated score); (ii) integrated TSVs of (a), (b), (c)
are byte-identical after removing nothing (same columns, same order, same probabilities); if the
web route's claim selection (GO/EC only, no operon claims) makes the inputs differ, the comparison is
made on the claims the web route actually integrates, re-run through (b), and the difference is
reported; (iii) the stage ledgers agree on stage names, statuses, reasons and change lists (timing
fields excluded); (iv) learned reports agree on every candidate key, probability, score status and
supporting sources; (v) provenance files agree record-for-record. Fail = any difference not explained
by a documented, intended route difference (for example the web route's claim scope).

R3. Multi-genome cohort comparison on the 18 real training-cohort genomes (the same claims, layouts,
cohort metadata, orthogroups, reaction references, ontology, logistic model and `lite=false`/six-prior
settings as the archived 2026-10-26 run, `benchmark/v3/results/nextflow-cohort-v1/`): CLI
`gspa integrate --workflow-config --batch-manifest` versus Nextflow `gspa-nf/cohort.nf` with the rc3
JAR. Pass = all per-genome integrated TSVs byte-identical, all cohort reports identical, all learned
candidate records and model hash identical, stage names/statuses identical with cohort inference
executed for every genome. Where the run is placed: on Shaheen through sbatch if Java 21 and the
Nextflow launcher can be made to run there; otherwise on this workstation with reduced threads and
`nice`, stated plainly (the archived run took 3 min 16 s). If neither is possible within the
session, the check is reported as not run, and point 2 is not met.

R4. States, each through the single CLI route, the annotation route and the web route (and the
Nextflow cohort route for cohort inputs), on the real small genome of R2:
  - requested: learned stage requested with a model file → `learned_integration` executed,
    `scoring_mode` learned or mixed, report written;
  - missing: learned stage requested with a model path that does not exist → the run fails with an
    explicit message naming the missing file, no integrated TSV is published, Nextflow publishes
    nothing;
  - no-work: a requested stage that finds nothing to do (reaction-local context without gaps, or
    dark matter with no candidates) → status `executed` or `skipped` with a reason that says so and
    an empty change list, never a silent success;
  - failed: metabolic refinement requested with a Python interpreter path that does not exist →
    stage `failed` with the tool error as reason, non-zero exit, ledger retained (web: job status
    `failed`, run log downloadable; Nextflow: task fails, work directory retained, nothing
    published).
Pass = every route shows the same state with the same reason text for the same input. Fail = a route
reports success where another fails, or hides a failed/missing state as zero changes.

### Verdict rule for point 2

Met only if R1, R2, R3 and R4 all pass. If the web route cannot run the full predictor stage here,
R2(c) may be run from the exported claims of (a) through the real portal integration code path
(`build_claims` → `workflow_config` → `gspa integrate`), and that substitution is stated; the point
is then reported as met for the integration path and not met for the predictor path of the web
route. Any cohort gap for the web route is reported as a scope limit, not hidden.

## Point 10: reporting

Bar: "Stage-wise baseline/final changes, evidence, alternatives and uncertainty; ontology/EC enzyme
classes; browser/export reconciliation."
Empirical: "Real production output rendered with correct no-work, failed, missing, learned/null/role
and partial-data states."

### Outputs (real production outputs unless marked controlled)

| Id | State(s) | Output | Origin |
|---|---|---|---|
| T1 | enforcement removals, taxon applied with lineage, priors executed | `price-taxon-rc2`, a TopEC cell, arm `stage1_taxon` | Shaheen Price run 2026-10-30/31 |
| T2 | stage 0: priors skipped, no organism taxon asserted, no-work stages | same cell, arm `stage0` | same |
| T3 | failed stage | `price-taxon-rc2` cell 101 (GCF_002797435.1), arm `stage1_taxon` | same |
| T4 | Stage 2B repairs (metabolic refinement executed with changes), FBA outputs | `price-stage2b-rc1`, the CLEAN cell with the most changed scores, arm `stage2b` | Shaheen Price run 2026-10-30 |
| T5 | learned (mixed scoring), cohort context, taxon not in hierarchy (strain taxon without lineage) | one genome of the R3 run (CLI batch output) | new run in this closure |
| T6 | null / abstained (model abstained on every candidate; probabilities null) | `benchmark/v3/results/gapsmith-observed-v2/final-cli-replay/logistic.tsv` | stored 2026-10-27 |
| T7 | role (candidate roles CATALYTIC / COMPLEX_SUBUNIT / UNSPECIFIED in the learned report) | `benchmark/v3/results/claim-role-bridge-v1/after-cli/logistic.tsv` | stored 2026-10-27 |
| T8 | partial data: protein-only input without coordinates; analyses not applicable / insufficient evidence | the web-route job of R2(c) (`analyses.json`, no GFF) | new run in this closure |
| T9 | missing: stage ledger absent; learned report absent while the ledger says the stage ran | controlled: copies of T5 with the ledger file, then the learned report, removed; labelled controlled | derived |

If a listed state does not occur in the named real output when inspected, the report says so and
uses the controlled input only for that state.

### Checks per output

C1. Packaged CLI: `java -jar gspa-3.0.0.jar visualize` on the output succeeds and writes the
browser HTML.
C2. Export reconciliation: browser `DATA.ledger.reconciliation.counts_match` is true with zero
probability mismatches and no net change without a stage action; the number of candidate rows in
the browser equals the rows in the integrated TSV; all input proteins are present (denominator).
C3. Stage-wise ledger: browser stage list equals the raw ledger (name, status, reason) and every
structured stage note is rendered (no `[object Object]`); failed stage shows its reason; skipped
stages show their reason; no-work stages show an empty change list, not a missing one.
C4. Enforcement removals (T1): number of `removed` changes in `symbolic_enforcement` equals the
rows of the enforcement-actions file for that cell (if present) and equals the number of candidates
above threshold in T2 that are absent or demoted in T1; the browser and portal show that count;
`reference_preparation` note shows `taxon_consistency_state` applied and the organism taxon.
C5. Stage 2B repairs (T4): the `metabolic_refinement` change list equals the per-candidate
difference between the stage 0 TSV and the stage 2B TSV of the same cell; the browser's metabolism
panel reports the repaired reactions; the portal `fluxes` and `model` artifacts are served.
C6. Learned/null/role (T5, T6, T7): every learned candidate identity, probability (including null),
score status, role and supporting sources in the browser equals the learned report; the portal
serves the exact report bytes; abstention is labelled as such, never as zero.
C7. Taxon not applied (T5): the `reference_preparation` note `taxon_consistency_state` reads
"taxon not in hierarchy: constraints not applied" and that text is visible in the rendered ledger.
C8. Enzyme classes: browser enzyme-class counts (complete EC, partial EC, GO catalytic) equal a
recomputation from the TSV and the EC-to-GO table.
C9. Missing (T9): portal ledger summary says unavailable (not zero) for the absent ledger; the
absent learned report under learned/mixed scoring is "unavailable", and the browser renders without
page errors.
C10. FastAPI: `/api/jobs/{id}` ledger summary agrees with the browser; `/api/jobs/{id}/learned`,
`/analyses`, `/fluxes`, `/model` return the stored bytes when present and 404 when absent.
C11. Chromium: results browser and portal page load with no page errors; the ledger tab text contains
the expected state words for each output (failed + reason, skipped, "model abstained",
"constraints not applied", removal count, repaired reactions).
C12. Partial data (T8): the browser says "Not measured" / not applicable for absent analyses and
coordinates, never zero; the portal lists the analyses with their recorded states.

### Verdict rule for point 10

Met only if C1 to C12 pass for every output where the state exists, with real outputs for every
state except those explicitly recorded as not occurring in any real output (then controlled, labelled).

## Failure handling

Any check that fails is reported as failed with the observed values; no bar is loosened and no output
is swapped for an easier one after inspection. Harness defects are fixed and rerun, with the failed
attempt kept in the evidence folder. Evidence: `benchmark/v3/results/routes-reporting-closure-v1/`
with `SHA256SUMS`; report: `docs/benchmarks/v3-routes-reporting-closure.md`.
