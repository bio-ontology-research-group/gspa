# Learned-scored output on every route, and a joint cohort job for the web service

Date: 2026-11-01. Pre-registered plan: `tasks/v3-routes-learned-cohort-2026-11-01.md`
(SHA-256 `96abb5042b80b799a1edc29d5c2e5de375516a19c037ce45645b4f33a874807c`, frozen before any check ran;
deviations are appended to the plan, dated). Evidence: `benchmark/v3/results/routes-learned-cohort-v1/`
(`SHA256SUMS`; files over 5 MB are listed in `LARGE_FILES_REMOVED.tsv` and kept under `~/tmp/routes-learned` on ws
and under `/scratch/project/k10018/gspa-v3/out/routes-learned-v1` on Shaheen). Code: `benchmark/v3/routes_closure/`
(`predict_contract_claims.py`, `run_learned_routes.py`, `build_learned_cohort_inputs.py`, `run_cohort_routes.sh`,
`web_cohort_job.py`, `compare_cohort_outputs.py`, `replay_reporting.py` extended with check C6b),
`benchmark/v3/shaheen/routes-learned-cohort3.sbatch`, and the service (`service/cohort.py`, `service/portal.py`,
`service/jobs.py`, `service/test_cohort.py`). Every number below is read from files in the evidence folder.

This closes the two gaps the previous closure left for acceptance points 2 and 10
([routes and reporting closure](v3-routes-reporting-closure.md)): no learned-scored output existed, and the web
service could not run a joint cohort.

Engine: the release-candidate-4 JAR of the previous closure, SHA-256 `08f8b5a1…`, pinned as a copy whose hash
is checked at every launch (`jvm/engine-jar.sha256`). No JVM source was changed. The JVM gate was still run:
`./gradlew build` then `./gradlew cleanTest test`, 1,584 tests, 0 failures, 6 optional skips
(`jvm/test-results-summary.json`). The JAR that `./gradlew build` rebuilt from the unchanged sources has a different
file hash: the comparison of its 31,230 entries with the pinned copy (`jvm/jar-diff-pinned-vs-rebuilt.json`) shows
one differing entry, the nested archive `downloads/gspa-resources.pyz`, which is repacked at build time; every class
and resource entry has the same CRC and size. Learned model: contract-3 logistic model
`benchmark/v3/results/predictor-support-v1/models/default/model_logistic.json` (`c9df93d6…`).

## Summary

| Point | Verdict against the frozen plan | Shown | Not shown |
|---|---|---|---|
| 2. Routes and state | **A1, A2, A3, A6, B1, B2 pass; A4 and A5 record the explicit abstention reason.** Three real Price genomes, genuine DIAMOND claims with the training predictor's contract (A1); the learned model rescores 25,838 to 28,201 candidates per genome on the single CLI route (A2), the CLI batch route is byte-identical to it (A3); CLI batch, Nextflow `cohort.nf` and the new web cohort job are byte-identical on every output file with cohort context executed and 27,870 to 30,385 candidates rescored per genome (A6, on ws and on Shaheen). The web cohort job is the CLI batch route itself, off by default, with the single-genome worker byte-identical (B1, B2). | Learned (mixed) scoring with scores, provenance, candidate identity and stage states shared by the single CLI, CLI batch, Nextflow and web cohort routes on matched inputs; joint cohort on all three multi-genome routes including the web service; requested state with a real model that scores. | The annotation route and the web single-genome upload route cannot produce a contract-bearing claim: their predictors (the in-process DIAMOND adapter, DeepGO-PlusPlus-Light) emit no `predictor_contract`, and the annotation route has no option to take supplied claims, so on those two routes the learned stage abstains explicitly on every candidate (A4, A5: 146,160 of 146,160, reason `candidate_predictor_contract_unavailable`). Orthogroups and cohort metadata remain supplied evidence on all three cohort routes (no route computes them). Predictor containers of Nextflow `main.nf` were not run. |
| 10. Reporting | **A7 passes.** Four real learned-scored outputs (three single-CLI runs and one genome of the web cohort job) rendered through the packaged `gspa visualize`, the FastAPI job-summary and artifact endpoints and Chromium: 103 checks, all pass, including the new C6b (every learned probability in the browser equals the report, the scored count equals the report's candidates minus unsupported, and the portal reports the learned stage as available with the report's own counts). | The mixed-scoring state on real outputs of this plan: scored and unsupported candidates side by side with null probabilities for the unsupported ones, the rescored count, cohort context, taxon applied (single CLI with lineage) and taxon not applied (cohort routes without lineage), reconciled candidate counts, enzyme classes recomputed independently, the portal serving the exact report bytes. | A learned-scored output with a `learned` (all-supported) scoring mode: none of these inputs yields one (the model abstains on 10,169 to 11,503 candidates per genome for the profile reasons tabulated below), so only `mixed` is shown; the states of the previous closure keep their earlier evidence. |

## What the contract requires, and how the claims were regenerated (A1)

The contract-3 models abstain with `candidate_predictor_contract_unavailable` unless every sequence claim of a
candidate carries `metadata.predictor_contract` (`gspa.predictor-contract/1`: tool, version, binary hash, reference
digest, parameter digest, adapter hash, score model). The production models were trained on one source only,
`swissprot_diamond_baseline`: DIAMOND 2.1.11 (binary `d6bea1f3…`) against the training cohort's genus-fenced
experimental Swiss-Prot reference (`reference.faa` `1d8c4ebb…`, labels `db8a8ac9…`), threads 4, 10 hits, e-value
1e-5, sensitive mode, adapter `benchmark/v3/homology_baseline.py` (`f84f9413…`). All four identity components
reproduce on ws today (Shaheen's DIAMOND is 2.1.16, a different binary), so the predictor ran on ws.

`predict_contract_claims.py` ran the unchanged adapter on three Price genomes and attached the contract built by
the training replay's own function; the contract identity (`ae823453…`) equals the one in the model's input-support
profiles for every genome (`a1/report.json`).

| Genome | Proteins | Proteins with hits | Claims | Contract equals the model's |
|---|---|---|---|---|
| GCF_000236665.1 | 3,498 | 2,225 | 38,173 | yes |
| GCF_000015245.1 | 3,662 | 2,461 | 41,203 | yes |
| GCA_003443575.1 | 3,814 | 2,369 | 41,448 | yes |

## Single-genome routes (A2 to A5)

Single CLI route: `gspa integrate --workflow-config --preset bacteria` with the claims, the genome layout,
the organism taxon with its lineage, the learned model, provenance and feature export (`a2/<genome>/workflow.json`).
Every run executed the learned stage with `scoring_mode = mixed`; no candidate was unsupported for a missing contract.

| Genome | Candidates | Rescored | Scored only | Unsupported | Learned probability range | Taxon |
|---|---|---|---|---|---|---|
| GCF_000236665.1 | 38,173 | 25,838 | 2,166 | 10,169 | 0.0006 to 0.780 | applied |
| GCF_000015245.1 | 41,203 | 27,571 | 2,342 | 11,290 | 0.0011 to 0.877 | applied |
| GCA_003443575.1 | 41,448 | 28,201 | 2,384 | 10,863 | 0.0008 to 0.870 | applied |

Unsupported reasons (single CLI route): `profile_lacks_both_classes_in_train_and_calibration` 9,560 / 10,565 / 10,246
(mostly GO molecular-function candidates: the training calibration split has no positive of that profile),
`unseen_source_function_role_or_availability_profile` 529 / 627 / 578, `sequence_score_outside_training_calibration_range`
80 / 98 / 39. `candidate_predictor_contract_unavailable`: 0 on every genome. Scored candidates per genome: 25,882 /
27,501 / 27,526 GO and 2,122 / 2,412 / 3,059 EC.

CLI batch route without cohort inference (A3, the three genomes in one `--batch-manifest`): for every genome the
integrated table, provenance, enforcement actions, stage ledger (elapsed seconds and the route's own feature-export
directory excluded), learned report (whole document, byte-identical) and the three feature-export files are
identical to the single CLI run (`a2/results-A2-A5.json`, `A3.all_identical_to_A2 = true`).

Annotation route (A4) and web single-genome upload route (A5), genome GCF_000236665.1, the routes' own predictors
(DIAMOND against the portal's reviewed reference, which carries no function mapping, and DeepGO-PlusPlus-Light):
both complete (annotation route exit 0; web job `complete`), both export 146,160 raw claims from
`deepgo-plusplus-light`, none with a predictor contract, and on both the learned stage is recorded as requested and
skipped with the reason "no candidates supported by the model input-coverage manifest; baseline scores retained",
146,160 of 146,160 candidates `candidate_predictor_contract_unavailable`, `scoring_mode = legacy_heuristic`. These
two routes cannot carry a contract-bearing claim today: the in-process DIAMOND adapter and the neural sidecar emit no
contract, and `gspa annotate` has no option to take supplied claims. That is the real reason, recorded, not a
mechanism failure. (The first web attempt used a Python without FastAPI and never submitted; kept as
`a5/failed-attempt-no-fastapi`.)

## Three cohort routes on matched inputs (A6)

Inputs: the A1 claims, the Price layouts, the three-genome subset of the Price cohort metadata and orthogroups
(with the layout bundle's `unplacedProteins` merged in, see deviations), reaction graph and EC aliases, ontology,
EC2GO, the learned model; settings are the Nextflow cohort module's defaults (lite, bacteria profile, four priors,
clade hold-out, universal context), feature export on, no taxon lineage (the module has no such column), so the
taxon state is "taxon not in hierarchy: constraints not applied" on all three routes (the strain taxids are not in
the bundled hierarchy). Routes: `gspa integrate --workflow-config --batch-manifest --batch-threads`; Nextflow 25.10.4
`gspa-nf/cohort.nf`; the web cohort job submitted through the real `POST /api/jobs` endpoint (kind `cohort`, a ZIP
bundle) and run by the real worker. Same pinned JAR on every route.

| Run | Route | Wall time | Peak memory | Output |
|---|---|---|---|---|
| ws | CLI batch (8 threads) | 3 min 7 s | 7.1 GB | 3 genomes, 120,824 rows |
| ws | Nextflow `cohort.nf` (one task, 8 cpus) | 3 min 17 s | 7.4 GB | 3 genomes, 120,824 rows |
| ws | web cohort job (worker, 8 threads, three browsers rendered) | 3 min 52 s | 7.7 GB | 3 genomes, 120,824 rows |
| Shaheen 15734950 | CLI batch (16 threads) | 2 min 50 s | 11.4 GB | 3 genomes |
| Shaheen 15734950 | Nextflow `cohort.nf` (one task, 16 cpus) | 2 min 48 s | 12.0 GB | 3 genomes |
| Shaheen 15734950 | web cohort job (worker, 16 threads) | 28 min 9 s | 10.7 GB | 3 genomes |

For every genome and every route the integrated table (38,173 / 41,203 / 41,448 rows), the stage ledger, the learned
report (model hash, counts, every candidate record; byte-identical), the cohort report, the provenance and the
three feature-export files are identical on ws (`a6-ws/comparison.json`, `identical: true` for all three genomes;
no route wrote suggestions or enforcement actions, which is recorded, not counted). The cohort stage executed for
every genome; the learned stage rescored 27,870 / 29,700 / 30,385 candidates with 10,303 / 11,503 / 11,063 unsupported
(reasons as above; 0 for a missing contract). The Shaheen run of the same three routes (job 15734950, one node,
same pinned JAR, Java 21.0.12, Nextflow 25.10.4, the web job submitted through the real endpoint) gives the same
verdict: all three routes exit 0 and every compared artifact is identical for all three genomes
(`a6-shaheen/comparison.json`). Across machines, the seven CLI-batch files that were transferred back before the
ssh session stalled (all five artifacts of GCA_003443575.1 and the table and cohort report of GCF_000015245.1) are
identical to the ws run (`a6-shaheen/a6-cross-machine.json`); the remaining files stay on Shaheen under
`/scratch/project/k10018/gspa-v3/out/routes-learned-v1/cohort3`, uncompared. The web job's wall time on Shaheen
(28 min against 3 min for the other two routes on the same node, and 3 min 52 s for the same job on ws) was not
diagnosed; its integration outputs are identical, so the time went into the per-genome browser rendering or the
Python worker on that node, which the retained run log can show.

The first ws attempt failed identically on all three routes at the cohort pass: the Price cohort metadata, built
from an older layout, did not declare protein AXT23839.1 of GCA_003443575.1 as unplaced, and the engine refuses a
cohort whose claims name a protein absent from the layout (`a6-ws/failed-attempt-metadata/`). This is an input
preparation defect of the harness, corrected by merging the layout bundle's `unplacedProteins`; the three routes
agreed on the failed state (one more state agreement, not a planned check).

## The web cohort job (B1, B2)

`service/cohort.py` adds the input kind `cohort`: one ZIP bundle with `cohort.tsv` (exactly the manifest columns of
`gspa-nf/modules/cohort.nf`), the claims and layout files it names, `cohort_metadata.json` and
`cohort_orthogroups.tsv`. The worker writes the typed workflow JSON with the Nextflow cohort module's settings and
the deployment's cohort resources (`GSPA_COHORT_RESOURCES`, validated against the module's allowed keys) and runs
`gspa integrate --workflow-config … --batch-manifest … --batch-threads …`, that is the CLI batch route; no cohort
inference exists in Python. Per genome it renders the results browser with the packaged `gspa visualize` and serves
the table, ledger, learned report, provenance, cohort report and browser at
`/api/jobs/{id}/genomes/{genome}/{artifact}`; `/api/jobs/{id}` carries one ledger summary per genome (no job-level
ledger is invented). On the real ws job every per-genome artifact was served byte-equal to the stored file, an
unknown genome and the single-job browser artifact are refused (404), download and log are served
(`a6-ws/web/api-surface.json`).

Default off: without `GSPA_COHORT_RESOURCES` the kind is refused with 400 and no folder is created. Single-genome
behaviour is unchanged: the literal-baseline configuration test still passes, `jobs.run` dispatches the three
single kinds to `run_full` without importing the cohort module, and `service/analysis.py` and
`service/integration_options.py` have the SHA-256 recorded before the change (`plan/service-before.sha256`,
checked by a test). Service tests: 66, all pass (`service/service-tests.log`); the 16 new ones pin the column list
and resource keys to the Nextflow module by parsing it, the composed JSON to the typed schema and the module's values,
the command to the batch route, the bundle validation (absolute and parent paths, unknown columns, duplicate or
unsafe ids, missing files, metadata mismatch, bad orthogroup header, unsafe archive member, non-ZIP) and the render
layout. The portal page offers no cohort upload control; the kind is an API feature of deployments that set the
resources.

One defect found and fixed on the way: the first web cohort job integrated all three genomes identically to the CLI
batch and then failed at the browser step, because the batch route names a genome's companion files after its table
(`integrated.tsv.provenance.json`) while the renderer reads the single-job layout (`gspa_out/provenance.json`).
`render_run_dir` now links every batch artifact under the single-job name; a regression test covers it; the job was
rerun (`a6-ws/web/job.id.failed-attempt-render` names the failed job, kept).

## Reporting of the learned-scored outputs (A7)

Harness `replay_reporting.py` (previous plan's checks C1 to C3, C6, C7, C8, C10, C11, plus C6b of this plan),
rendered with the pinned JAR, the FastAPI portal and headless Chromium.

| Case | Output | States shown | Reconciled numbers | Checks |
|---|---|---|---|---|
| single CLI, GCF_000236665.1 | A2 | mixed scoring, learned executed, taxon applied, enforcement, feature export | 32,264 rows; 28,004 candidates with a learned probability equal to the report, 10,169 null; portal: learned available, counts equal the report, "25,838 rescored" shown in the browser | 26 pass |
| single CLI, GCF_000015245.1 | A2 | same | 34,677 rows; 29,913 scored, 11,290 null | 26 pass |
| single CLI, GCA_003443575.1 | A2 | same | 34,896 rows; 30,585 scored, 10,863 null | 26 pass |
| web cohort job, GCF_000236665.1 | A6 (ws) | cohort context executed, mixed scoring, taxon not in hierarchy, feature export | 38,173 rows; 27,870 scored, 10,303 null; "Cohort context" and "constraints not applied" visible in the ledger tab | 25 pass |

The first render attempt failed one check per case on the harness's own expectation text ("rescored 25838" instead
of the browser's "rescored 25,838"); kept as `a7/batch-a2-failed-attempt-expectation-format`, the expectation was
corrected and the batch rerun (no rendering defect).

## Limits

* Learned scoring is shown as `mixed`, never as `learned`: the model abstains on about a quarter of the candidates
  for profile reasons that are properties of the training data, not of the routes.
* The annotation route and the web single-genome upload route abstain for a real reason (no contract-bearing
  predictor on those routes); closing that would need the in-process adapters to emit contracts and the models to be
  trained on them.
* Orthogroups and cohort metadata are supplied evidence on all three cohort routes; none computes them.
* These are route-parity and reporting checks on real outputs, not biological validation; the twelve-item
  acceptance record remains authoritative for the whole-item verdicts.
