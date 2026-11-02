# Routes and reporting: closure checks for acceptance points 2 and 10

Date: 2026-11-01. Pre-registered plan: `tasks/v3-routes-reporting-closure-2026-11-01.md`
(SHA-256 `02050c93ba53b360a38f7a6827e64fc6d97befd47e4dd9f9e6c91cdc00705a9b`, frozen before any check ran).
Evidence: `benchmark/v3/results/routes-reporting-closure-v1/` (`SHA256SUMS`; files over 5 MB are listed with their hashes in `LARGE_FILES_REMOVED.tsv` and kept under `~/tmp/routes-closure` on ws). Code:
`benchmark/v3/routes_closure/` and `benchmark/v3/shaheen/routes-cohort18.sbatch`. Every number below is
read from files in the evidence folder; Shaheen job 15733547 holds the bulky cohort outputs
(`/scratch/project/k10018/gspa-v3/out/routes-closure-v1/cohort18`).

Engines: release candidate 3 (`8391d73f…`) for the Price outputs and the Shaheen cohort run; a rebuilt JAR
(`08f8b5a1…`, called rc4 below) after one defect fix in the annotation command. An entry-by-entry
comparison of the two JARs (`jvm/jar-diff-rc3-rc4.json`: 31,223 vs 31,230 entries) shows only
`gspa.cli.AnnotateCommand` classes differ (8 changed, 8 removed and 15 added inner classes of that one
command); the integration engine, the Nextflow route and the visualization resource are byte-identical.

## Summary

| Point | Verdict against the frozen plan | Shown | Not shown |
|---|---|---|---|
| 2. Routes and state | **Met for the integration routes; the web route remains single-genome.** R1 passed; R2: annotation, single CLI and web routes agree on the real genome (CLI and web byte-identical; annotation identical to the CLI run on the same raw claims, byte for byte); R3: 18 real genomes, CLI batch and Nextflow byte-identical for every file; R4: requested, missing, failed and no-work states agree across the four routes. | Scores, provenance, candidate identity and stage states shared by all routes on matched inputs; joint cohort on CLI batch and Nextflow; mixed suggesters (reaction-local context, cohort recurrence, learned stage) and four-source claims. | A multi-genome web job (the service has no orthogroup step, so the joint cohort stays batch-only); predictor-level parity through Nextflow `main.nf` (containers not available here); a learned-scored output, because every claim file at hand lacks the predictor contract the contract-3 models require, so the learned stage abstains everywhere (explicitly). |
| 10. Reporting | **Met for every state that exists in a real output; one state controlled.** 11 cases, 192 checks, all pass: enforcement removals, taxon applied and not applied, stage-0 skipped stages, Stage 2B repairs with model and fluxes, failed stage, no-work stage, learned abstention with null probabilities, roles, cohort context, feature export, partial data, missing ledger, missing learned report. | Packaged CLI render, FastAPI summary and artifact endpoints, Chromium results browser and portal page, counts reconciled with the raw files, enzyme classes recomputed independently. | A learned-scored (mixed) output: none exists any more (the 2026-10-27 stored outputs were removed by the 2026-10-31 clean-up and no current model scores these claims); that state keeps only the earlier reporting-parity evidence. The missing-ledger and missing-learned-report states are controlled copies of real outputs, labelled as such. |

## Defects found and fixed on the way

1. **The web worker crashed at the quality step of every default job.** `service/analysis.py` referenced an
   undefined name (`pathway_args`) when building the `gspa evaluate` command, so a job with the pathway
   analysis enabled (the default) stopped with a `NameError` after integration. Fixed by building the
   pathway arguments from the same KEGG references the integration used; a new test parses the worker
   with `ast` and fails on any undefined name.
2. **`gspa annotate --preset X --integrate` ignored the preset's integration defaults.** Only
   `gspa integrate --preset X` applied `WorkflowPresets.integrationDefaults` (priors, stage-2 enforcement
   and its mode, essential profile, organism taxon, `lite`). The annotation route therefore ran the shared
   workflow with the built-in defaults on the same claims: no stage-2 enforcement, no organism taxon and
   full ELK ontology inference (the first run on the real genome spent its time in the `GoReasoner`,
   `r2/routeA-rc3-divergence/run.log`), while the CLI and web routes, which pass the preset to
   `gspa integrate`, ran enforcement in lite mode. `AnnotateCommand.applyIntegrationPreset` now applies the
   same defaults with the same precedence (preset < integration-config file; an explicit `--taxon` or
   `--infer-taxon` keeps its own organism). A new Spock case in `WorkflowRouteParitySpec` compares the
   resolved settings of the two commands. JVM gate after the fix, all test tasks re-executed: 1,584 tests,
   0 failures, 6 optional skips (`jvm/test-results-summary.json`, `jvm/jvm-gate-rerun.log`).
3. **Web route options (plan R1).** `service/integration_options.deployment_options` reads
   `GSPA_LEARNED_MODEL`, `GSPA_FEATURE_EXPORT` and `GSPA_EXPERIMENT_MANIFEST` (with mode, strategy, sidecar
   and Python) and passes `learnedModelFile`, `featureExportDir` and the `experiment*` keys into the same
   typed workflow JSON the CLI reads. With nothing set the JSON is byte-identical to the one written before
   (`test_defaults_are_unchanged_without_deployment_options` compares against a literal). The container
   paths became overridable (`GSPA_APP_ROOT`, `GSPA_FBA_PYTHON`, `GSPA_LIGHT_MODEL`) so the same worker
   runs outside the image; a tool that is not installed is recorded as "not installed" in the manifest
   instead of crashing the provenance step. The deployment options are recorded in each job's
   `manifest.json`. Service tests: 49, all pass (`service/service-tests.log`).

The joint cohort cannot be offered to a web upload: cohort context needs orthogroups of the uploaded
proteins against the cohort, which the service does not compute (`gspa-nf`/CLI take them as supplied
evidence). It stays a CLI-batch and Nextflow feature; the web route now carries every single-genome option.

## Point 2: routes and state

### R2. One real genome through the annotation, CLI and web routes

Genome GCA_000011125.1 (1,700 proteins, training cohort, archaea preset). Predictors: DIAMOND against the
cohort's reviewed reference database (0 claims: that database carries no function mapping, so there is no
second predictor source here) and DeepGO-PlusPlus-Light (57,570 GO claims). Learned model requested through
the new deployment option. Wall time 2.6 to 3.9 min per route, up to 5.2 GB.

| Comparison | Claims | Integrated table | Stage ledger | Learned report | Provenance | Enforcement actions |
|---|---|---|---|---|---|---|
| web worker vs CLI `integrate` on the web job's selected claims | 57,570 = 57,570 | byte-identical (`fe0f44f9…`, 52,873 rows) | identical | identical | identical | identical (5,917 actions) |
| annotation route vs CLI `integrate` on the annotation route's raw claims | 57,570 = 57,570 | byte-identical (`d1f4d058…`) | identical | identical | identical | identical |
| annotation route vs web worker (as run) | identical as sets | 52,873 of 52,873 rows differ in the `go_aspect` column only; every posterior identical | same stage states and change contents; record order differs | identical except `feature_digest` in 417 of 57,570 candidates | identical records, different order | identical as a multiset |

The one cross-route difference is the web route's intended edit: it fills a GO claim's missing aspect from
the reference before integration (`build_claims`), the annotation route integrates the claims as the
predictor produced them (no aspect). Scores, candidate identities, stage states, roles and actions do not
depend on it. Record order differs because the annotation route integrates claims from memory and the CLI
reads them from a file. Files: `r2/comparison.json`, `r2/comparison-rawclaims.json`.

Mixed sources on the single-genome routes: the portal's integration path (claim selection, workflow
config, `gspa integrate --preset`) and the plain CLI flag route were also run on real four-source Price
claims (GCF_000236665.1: CLEAN 5,245, DeepProZyme 2,762, EnzBERT 106, MAPred 69 claims): integrated table,
learned report, provenance and 18 enforcement actions byte-identical; the ledgers differ only in elapsed
seconds (`r2e/`).

### R3. Eighteen real genomes: CLI batch versus Nextflow (Shaheen job 15733547)

Same inputs as the archived 2026-10-26 run (claims, layouts, cohort metadata and orthogroups, reaction
references, ontology, KEGG tables, reasoner cache; `r3/inputs.sha256`), six priors, clade hold-out,
universal context, full ontology inference, feature export on. Deviation from the plan: the archived
run's logistic model is a contract-2 artifact that rc3 rejects, so the contract-3 model
`predictor-support-v1/models/default/model_logistic.json` (`c9df93d6…`) was used on both sides.

| Route | Wall time | Peak memory | Output |
|---|---|---|---|
| `gspa integrate --workflow-config --batch-manifest` (16 threads) | 4 min 13 s | 28.5 GB | 18 genomes, 484,202 rows |
| `gspa-nf/cohort.nf` (Nextflow 25.10.4, one task, 16 cpus) | 4 min 33 s task | 15.4 GB | 18 genomes, 484,202 rows |

For every genome the integrated table, cohort report, learned report (model hash, counts and every
candidate record), stage ledger (names, statuses, reasons and change lists; elapsed seconds excluded),
quality summary, rescue table, provenance and the feature-export files are identical
(`r3/comparison.json`, 18 of 18 `identical: true`). Cohort context executed for every genome; the learned
stage abstained for every genome (reason recorded; below). No suggestion table was written by either route
(both requested it; the stage found nothing to write), which the record states rather than counting as agreement.

### R4. Requested, missing, failed and no-work states

Single CLI route, annotation route and web worker on the real genome of R2; Nextflow and the CLI batch on a
two-genome real sub-cohort (GCA_000011125.1, GCA_000011585.1) of the training cohort.
(`r4/agreement-combined.json`; the first Nextflow attempt in `r4/results.json` failed on a harness error,
cohort metadata describing 18 genomes for a 2-genome cohort, identically on both routes; it was rerun with
the metadata and orthogroups subset, `r4/nextflow-rerun/`.)

| State | Input | CLI | Annotation | Web | Nextflow | Agree |
|---|---|---|---|---|---|---|
| requested | contract-3 model | executed run; `learned_integration` skipped: "no candidates supported by the model input-coverage manifest; baseline scores retained"; table written | same | same, job complete | same; published outputs byte-identical to the CLI batch | yes |
| missing | model path that does not exist | `reference_preparation` failed: "learnedModelFile does not exist: …"; downstream not run; exit 1; no table | same | same; job failed, log downloadable | Nextflow refuses to stage the missing file and never starts the engine; nothing published; the CLI batch on the same manifest writes the same failed ledger | yes (Nextflow fails earlier, documented) |
| failed | contract-2 model file | `reference_preparation` failed: "ModelUnavailableException: artifact feature contract gspa.context-features/2/… does not match …/3/…"; exit 1; ledger retained, no table | same | same; job failed | task fails; two failed ledgers retained in the work directory; nothing published; CLI batch identical | yes |
| no-work | reaction-local context requested, no metabolic gaps | `reaction_context` skipped: "no metabolic gaps to evaluate"; 0 changes | same | not exposed on the web route (reaction context is not a web option) | not run (no overlay key) | yes (CLI, annotation) |

### Why no learned-scored output exists

Every claim file available to this closure (the 18-genome homology claims, the Price pooled claims, the
DeepGO-PlusPlus-Light claims the annotator exports) lacks the structured predictor contract that the
contract-3 models require; the learned stage records `candidate_predictor_contract_unavailable` for every
candidate and keeps the baseline scores. This is the explicit abstention state, exercised on all routes,
but it means neither point can show a learned-scored (mixed) output today. The stored 2026-10-27 mixed
outputs used in the earlier reporting-parity check were removed in the 2026-10-31 scratch clean-up
(dangling links in `~/.cache/gspa-v3-data/model-input-support-browser-v2-final/`).

## Point 10: reporting

Eleven cases (192 checks) rendered through the packaged `gspa visualize`, the FastAPI job-summary and artifact
endpoints (`/api/jobs/{id}`, `/learned`, `/analyses`, `/fluxes`, `/model`, `/ledger`, `/log`) and Chromium
(results browser with every stage note expanded, and the portal page served through the real endpoints).
Checks C1 to C12 of the plan; harness `benchmark/v3/routes_closure/replay_reporting.py`; results in
`p10/batch1`, `p10/batch2`, `p10/batch3`, `p10/batch3-failed-web`. The first harness attempt failed on
two harness defects (browser protein key, dangling stored input) and is kept in `p10/batch1-failed-attempt`;
`p10/batch1-rc3` is the same batch rendered with rc3 (identical results, visualization unchanged).

| Case | Output (origin) | States shown | Key reconciled numbers | Checks |
|---|---|---|---|---|
| T1 | Price cell 0 (TopEC, GCF_000236665.1), arm with taxon fix, Shaheen job 15721368 | enforcement removals, taxon applied with lineage, learned not requested | 37 removals in the ledger = 37 rows lost between stage 0 and the fixed arm = 37 shown; 3,444 rows reconciled; 3,493 input proteins; 37 complete-EC enzyme changes recomputed | 21 pass |
| T2 | same cell, stage 0 arm | "taxon not in hierarchy: constraints not applied", enforcement and priors skipped, no net change | 3,481 rows; state text visible in the rendered ledger | 17 pass |
| T4 | Price Stage 2B cell 0 (TopEC), Shaheen job 15712377 | Stage 2B metabolic refinement executed with repairs; model and fluxes served | 6 ledger changes = 6 posterior differences between the stage-0 and Stage 2B tables; 123 repaired reactions in the note; `metabolic_model.xml` and `reactions.tsv` served | 19 pass |
| T3 | real failed web job (R4, contract-2 model) | failed stage with its reason, downstream stages not run, job failed | portal status failed with the worker's message; run log served with the exception; browser, download, learned (404) and ledger (409) refused; failed ledger retained without a table | 9 pass |
| T10 | real CLI run (R4 no-work) | requested stage with nothing to do | `reaction_context` skipped "no metabolic gaps to evaluate"; 52,873 rows reconciled | 18 pass |

A failed run writes no prediction table, so the results-browser harness cannot render it (`p10/batch3`,
case `T3-failed-stage`, recorded as a harness limit); the failed state is checked where production shows
it, on the portal (T3).
| T5 | new 18-genome cohort run, GCA_000007005.1 | cohort context executed, learned requested and abstained (26,328 null probabilities), feature export executed, taxon applied | 26,328 rows and candidates; roles CATALYTIC / COMPLEX_SUBUNIT / UNSPECIFIED; 2,994 input proteins | 21 pass |
| T6 | stored gapsmith replay (protein-only) | model abstained on every candidate; no coordinates | 3,469 candidates all null; 2,661 input proteins, 802 annotated; "Not measured" for sequence length and contigs | 18 pass |
| T7 | stored claim-role-bridge replay | candidate roles in the learned report | roles CATALYTIC / COMPLEX_SUBUNIT / UNSPECIFIED on 3,469 candidates | 18 pass |
| T8 | real web job (rc4, GCA_000011125.1 proteins) | protein-only input; analyses skipped by the user; symbolic enforcement executed (relax, 5,917 actions); learned abstained | 52,873 rows, 1,700 proteins, 847 annotated; 57,570 candidates null; "Not measured" | 21 pass |
| T9a | controlled: T1 without its ledger | stage ledger missing | browser and API say unavailable, "not zero" | 13 pass |
| T9b | controlled: T5 without its learned report | learned report missing | browser and API: not applicable / unavailable, no candidates, no page error | 17 pass |

In every case the browser's reconciliation reports matching candidate counts, zero probability
mismatches and no net change without a stage action; the stage list equals the raw ledger (names,
statuses, reasons, structured notes); learned candidate identities, probabilities (including null), score
status, roles and supporting sources equal the report; the portal serves the exact report bytes; and no
page error occurs. Enzyme classes of the final changes recomputed from the EC pattern and an independent
is_a closure over `go-basic.obo` agree with the browser.

No real output with a failed stage exists among the 264 Price cells (a grep hit on cell 101's 6 GB ledger
was inside a note, not a stage status; `r4`/`p10/batch3-failed-web` therefore use the real failed runs of R4).

## Limits

* The web route's predictor stage ran here with DeepGO-PlusPlus-Light and a DIAMOND database without
  function mapping; PSORTb, AMRFinderPlus and antiSMASH were not available and those analyses were
  disabled at submission (recorded as "skipped" and "not installed"). The portal worker ran under Python
  3.10 with a four-line `hashlib.file_digest` shim in the harness (the image uses a newer Python).
* Nextflow `main.nf` (predictors in containers) was not run; Nextflow parity is shown for the cohort and
  integration routes on supplied claims, which is where the production pipeline enters the shared workflow.
* No learned-scored output could be produced (predictor contracts absent from every available claim file);
  the learned state is shown as explicit abstention only.
* The Stage 2B reporting case is a TopEC cell (3,481 rows), not the CLEAN cell with the most changed top
  calls named in the plan: CLEAN cells are dense (3 GB table, 9 GB ledger) and cannot be rendered here.
* These are mechanism and reporting checks on real outputs, not biological validation. The twelve-item
  acceptance record remains authoritative for the whole-item verdicts.
