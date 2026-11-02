# Pre-registered plan: learned-scored output on every route, and a joint cohort job for the web service

Date: 2026-11-01. Author: Fable agent for Robert Hoehndorf. Frozen before any check is run; the SHA-256
of this file is recorded in the evidence folder (`benchmark/v3/results/routes-learned-cohort-v1/plan/`)
and in the report (`docs/benchmarks/v3-routes-learned-cohort.md`). Nothing above the "Deviations"
heading is changed after the hash is taken; anything that has to change is appended there, dated.

This plan closes the two gaps left by `tasks/v3-routes-reporting-closure-2026-11-01.md` for acceptance
points 2 (routes and state) and 10 (reporting) in `tasks/v3-acceptance.md`:

A. no real learned-scored output existed, because every claim file at hand lacked the predictor
   contract that the contract-3 learned models require;
B. the web service could not run a joint multi-genome cohort.

The bars are the ones in `tasks/v3-acceptance.md` (items 2 and 10) and are not changed after seeing results.

## What the contract requires (found in the code, not assumed)

* `gspa.integration.learned.ModelInputSupport.unsupportedReason` returns
  `candidate_predictor_contract_unavailable` unless every sequence claim of a candidate carries
  `metadata.predictor_contract` (`gspa.integration.learned.PredictorContract`: schema
  `gspa.predictor-contract/1`, keys `tool`, `tool_version`, `tool_sha256`, `reference_sha256`,
  `parameters_sha256`, `adapter_sha256`, `score_model`, all 64-hex digests), read in
  `ContextFeatureExtractor.build` from `claim.metadata.predictor_contract`.
* A candidate is then supported only when its (sources, availability flags, context mode, contract
  identities) profile is one the model saw in train and calibration with both classes, and its
  sequence scores lie inside the profile's observed range (`unseen_predictor_contract`,
  `unseen_source_function_role_or_availability_profile`, `profile_lacks_both_classes_in_train_and_calibration`,
  `sequence_score_outside_training_calibration_range` otherwise).
* The production contract-3 models (`benchmark/v3/results/predictor-support-v1/models/default/model_logistic.json`
  and the family-excluded refit) were trained on one source only, `swissprot_diamond_baseline`, with one
  contract: DIAMOND `diamond version 2.1.11`, tool SHA-256 `d6bea1f3…`, reference digest `2ef431fe…`
  (reference.faa `1d8c4ebb…` + reference_labels.tsv `db8a8ac9…` of the training cohort's genus-fenced
  Swiss-Prot reference), parameters digest `f8a5d09a…` (threads 4, 10 hits, e-value 1e-5, sensitive,
  outfmt qseqid sseqid pident qcovhsp evalue bitscore, transfer mode `homology`), adapter SHA-256
  `f84f9413…` (= `benchmark/v3/homology_baseline.py`), score model
  `best_reference_label_identity_times_capped_query_coverage_round6/v1`. Supported profiles need a genome
  layout (`ctx.layout.available` = 1), no predictor-status block, GO candidates of aspect BP/CC with role
  unspecified or EC candidates with role catalytic, homology score in [0.034, 0.941] (per profile), and for
  EC also some cohort-context profiles.
* All four identity components reproduce on the workstation ws today: `/usr/bin/diamond` is 2.1.11 with
  the same SHA-256, `benchmark/v3/homology_baseline.py` has the adapter hash, the reference files are in
  `~/.cache/gspa-v3-data/training-cohort/baseline/`, and the parameter digests recompute to the recorded
  values. Shaheen's DIAMOND is 2.1.16 (different binary hash), so the predictor runs on ws; that is a
  tool-identity requirement, not a convenience. Integration runs are placed as stated per check.

## Part A: genuine contract-bearing claims and learned-scored output on every route

### Inputs

Three real Price benchmark genomes, the smallest of the 22 (`~/.cache/gspa-v3-data/price-layout-v4/`,
`price-cohort/`): GCF_000236665.1 (3,498 proteins, NCBI taxid 640081), GCF_000015245.1 (3,662, taxid
326297), GCA_003443575.1 (3,814, taxid 221822). They are not training-cohort genomes; their genera are
excluded from the reference by construction (`homology_baseline.py`). Protein FASTA per genome from
`price-layout-v4/all_proteins.faa` (namespace `<assembly>__` stripped), layout from
`price-layout-v4/<genome>.layout.tsv`, cohort metadata and orthogroups as the subset of
`price-cohort/cohort_metadata.json` and `cohort_orthogroups.tsv` for these three genomes, lineages from
`benchmark/v3/taxonfix/price22_taxon_lineages.tsv` where a route can take them.

Engine: `gspa-cli/build/libs/gspa-3.0.0.jar`, SHA-256 `08f8b5a132f4f1e5bf458f45b7604b1aa3b94388c96f170682bd2c27f2559c5a`
(release candidate 4 of the previous closure; no JVM source is newer than it). No JVM source is changed
by this plan; if one has to be, the JVM gate is rerun and the JAR rebuilt and re-hashed before the checks.
Model: `benchmark/v3/results/predictor-support-v1/models/default/model_logistic.json` (contract 3), hash recorded.
References: `~/tmp/routes-closure/r2/resources/{go/go.owl,ec2go/ec2go,kegg_pathways.tsv,kegg_modules.tsv,kegg_module_defs.json,reasoner-cache}`,
reaction graph `~/Public/software/gapsmith/data/{seed_reactions.tsv,seed_Enzyme_Class_Reactions_Aliases_unique.tsv}`
(the same files the previous closure used; on Shaheen their mirrors under `/project/k10018/gspa-v3/`), all hashed.

### A1. Regenerate the predictor claims with the contract (ws)

New script `benchmark/v3/routes_closure/predict_contract_claims.py`: per genome, runs
`homology_baseline.run_diamond` (the unchanged adapter; `makedb` on the training reference.faa, then
`blastp` with exactly the training parameters), `claims_from_hits(..., 'homology')`, `write_claims_jsonl`
(source `swissprot_diamond_baseline`), and attaches `metadata.predictor_contract` built by
`replay_predictor_contract.contract(...)` (the same function the training replay used). It records every
input and tool hash and the hits file hash.
Pass = for all three genomes: the contract's `identity_sha256` equals the digest of the contract stored in
the model's input-support profiles; the claims file is non-empty; reference and adapter hashes equal the
recorded ones. Fail = any hash differs (then no learned scoring is attempted with that file).

### A2. Single CLI route with the learned model (ws), one run per genome

`gspa integrate --workflow-config cfg.json --preset bacteria` with `claimsFile` (A1), `genomeLayoutFile`,
`organismTaxon NCBITaxon_<taxid>`, `taxonLineage` (lineage table), `learnedModelFile`, provenance on,
`featureExportDir` on, ontology/EC/KEGG references, no cohort, no operons, no external context.
Pass = `learned_integration` executed; ledger `scoring_mode` is `mixed` or `learned`; the learned report
`counts.rescored + counts.added > 0` on at least one genome; zero candidates with
`unsupported_reason = candidate_predictor_contract_unavailable`; every unsupported reason is tabulated.
Fail = the learned stage abstains on every candidate; then the tabulated reasons are reported as the real
reason and the point stays not met.

### A3. CLI batch route without cohort (ws): the three genomes in one `--batch-manifest`

Same settings as A2 (manifest columns taxon, genome_layout_path, taxon_lineage; `--preset bacteria`).
Pass = for each genome the integrated TSV is byte-identical to A2, the learned report identical (model hash,
counts, every candidate record), the stage ledger identical except `elapsed_seconds`, provenance identical,
feature export files identical. Fail = any other difference.

### A4. Annotation route (ws), one genome (GCF_000236665.1)

`gspa annotate --proteins <faa> --preset bacteria --integration-config cfg.json --claims-out ...` with the
same learned model, DIAMOND against the portal's Swiss-Prot database and DeepGO-PlusPlus-Light (the
predictors the route has). The annotation route has no option to take supplied claims and its DIAMOND
adapter emits no predictor contract, so the expected outcome is explicit abstention.
Pass (as documentation of a real reason) = the run completes; the ledger records the learned stage as
requested with the explicit abstention reason; 100 % of candidates carry `unsupported_reason =
candidate_predictor_contract_unavailable` (or another single explicit reason, reported verbatim); the
exported raw claims contain no `predictor_contract`. This check cannot make the route "learned-scored"; it
is reported as not shown, with this reason. Fail = a silent success (learned stage marked executed with
0 unsupported while no contract is present) or a crash.

### A5. Web single-genome upload route (ws), the same genome

`service.analysis.run_full` through the real worker on the protein FASTA (same environment as the previous
closure: `GSPA_LEARNED_MODEL` set, analyses other than homology/names disabled).
Pass / fail = as A4 (same predictors, same expected reason); the job must finish `complete` with the
abstention visible in `integrated.tsv.learned.json` and the ledger.

### A6. Cohort routes on the three genomes: CLI batch, Nextflow `cohort.nf`, web cohort job (Shaheen, one node, one job)

Matched inputs: A1 claims, layouts, the three-genome subset of the Price cohort metadata and orthogroups,
reaction graph and EC aliases, ontology, EC2GO, the learned model, feature export on. Settings are the
Nextflow cohort module's defaults (`lite` true, essential profile bacteria, priors
`essentiality,coherence,gap_filling,genomic_context`, hold-out CLADE, context UNIVERSAL), no taxon lineage
(the cohort module has no such column; all three routes omit it, and the resulting taxon state is recorded).
Routes: (i) `gspa integrate --workflow-config --batch-manifest --batch-threads N`; (ii) `gspa-nf/cohort.nf`
with Nextflow 25.10.4; (iii) the new web cohort job (`service.cohort.run_cohort`) on the uploaded bundle,
run by the real worker code with `GSPA_BIN` the same JAR launcher. Same JAR (hash checked in the job).
Pass = for each genome all three routes give byte-identical integrated TSVs, identical learned reports
(model hash, counts, every candidate record), identical stage ledgers (except `elapsed_seconds`),
identical cohort reports, provenance, suggestions (or identically absent) and feature-export files; the
cohort stage executed for every genome; and the learned report shows `counts.rescored + counts.added > 0`
on at least one genome (if not, the unsupported reasons are reported and the learned-scored cohort output
is "not shown"). Fail = any difference not explained by a documented intended route difference.
Also recorded, informational (not pass/fail): the CLI batch of (i) rerun on ws and compared with the
Shaheen copy, to show the engine is machine-independent.

### A7. Reporting of the learned-scored outputs (ws)

Through `benchmark/v3/routes_closure/replay_reporting.py` (extended, not loosened): the A2 output of each
genome and one genome of the web cohort job (A6 iii) are rendered through the packaged `gspa visualize`,
the FastAPI job-summary and artifact endpoints and Chromium. Checks C1, C2, C3, C6, C8, C10, C11 of the
previous plan apply, with `expect.scoring_mode = mixed` (or `learned`), `expect.null_candidates` = the
report's `counts.unsupported`, and a new check C6b: the number of candidates with a non-null probability in
the browser equals the report's `counts.candidates - counts.unsupported`, each probability equal to the
report's, and the portal summary shows the mixed-scoring label and the rescored count.
Pass = every check passes for every rendered output. Fail = any check fails (reported with observed values).

## Part B: a multi-genome cohort job in the web service

Design (service only; no JVM change): a new input kind `cohort` accepts one ZIP bundle containing a
manifest `cohort.tsv` with exactly the columns of `gspa-nf/modules/cohort.nf` `cohortColumns()`
(`genome_id, claims_path, out_path, taxon, essential_profile, operons_path, gaps_path, max_iter,
enforce_stage2, enforce_mode, enable_priors, genome_layout_path, orthogroups_path, cluster_consensus_path,
fba_medium, fba_gram`; `out_path` ignored, the service owns outputs), the per-genome claims and layout files,
`cohort_metadata.json` and `cohort_orthogroups.tsv`. The job runs **the CLI batch route itself**:
`gspa integrate --workflow-config <json> --batch-manifest <tsv> --batch-threads <n>` with a workflow JSON
composed of the same keys and defaults the Nextflow cohort module composes (`lite`, `essentialProfile`,
`enabledPriors`, `cohortConfig`, `cohortInference`, `provenanceFile`, `suggestionsOut`, optional
`featureExportDir`, plus the deployment's cohort resources with the allowed-key list of `cohortResources()`).
No cohort inference is implemented in Python. Per genome the service renders the results browser with the
packaged `gspa visualize` and serves ledger summaries through the existing `ledger_summary.summarize`.

Default off: the kind is accepted only when the deployment sets `GSPA_COHORT_RESOURCES` (a JSON file of
cohort resources). Without it, `POST /api/jobs` with `kind=cohort` returns 400 and creates nothing.

B1. Tests (`service/test_cohort.py`, run with the other service tests). Pass = all pass:
  - the column list equals `cohortColumns()` parsed from `gspa-nf/modules/cohort.nf`, and the allowed
    resource keys equal `cohortResources()`'s list;
  - the composed workflow JSON uses only keys of `integration-workflow-config.keys.json` and equals the
    Nextflow module's settings for the same inputs (same key set and values, paths aside);
  - the argv is `integrate --workflow-config … --batch-manifest … --batch-threads …` (no private flags);
  - bundle validation rejects absolute paths, `..`, unknown columns, duplicate or unsafe genome ids,
    missing claims/layout/metadata/orthogroups, and does not leave a job folder behind;
  - default off: with the variable unset the kind is rejected with 400 and no folder is created;
  - single-genome behaviour unchanged: the existing literal test `test_defaults_are_unchanged_without_deployment_options`
    still passes; `jobs.run` dispatches `gff3`/`genome`/`proteins` to `run_full` exactly as before (mocked)
    and never imports the cohort module for them; and the SHA-256 of `service/analysis.py` after the change
    equals its SHA-256 before this plan (recorded in the evidence folder: a proof that the single-genome
    worker code is byte-identical).
B2. The full service test suite passes (`python -m unittest discover -s service`), and the JVM gate
    (`./gradlew build`) passes unchanged (no JVM source touched; the run is still made and recorded).
B3. A6 above is the empirical check of the web cohort job against the CLI batch and Nextflow.

## Verdict rules

Point 2 (routes and state): "Met" only if A1, A2, A3, A6 and B1–B2 pass and A4/A5 record the explicit
abstention reason; the web route's single-genome upload path then remains documented as unable to produce
contract-bearing claims (its predictors emit none), which is stated in the row text. If A2 or A6 shows no
learned scoring, the point stays "Not met (mechanism verified)" with the tabulated reasons.
Point 10 (reporting): "Met" only if A7 passes on a learned-scored output from a real run of this plan and
every state of the previous closure keeps its evidence. Otherwise "Not met (mechanism verified)" with the
missing state named.

## Rules carried over

No label is read by any run; unconfirmed labels are never negatives. Harness defects are fixed and rerun
with the failed attempt kept. Heavy runs go to Shaheen via Slurm (account k10018, `--exclude=nid00022,nid00497,nid01033`,
REPO_ROOT passed explicitly); nothing is written to `/tmp`; the jobs and files named in the task are not
touched. Evidence: `benchmark/v3/results/routes-learned-cohort-v1/` with `SHA256SUMS` (files over 5 MB
listed with hashes and kept under `~/tmp/routes-learned` on ws); report:
`docs/benchmarks/v3-routes-learned-cohort.md`; one progress line appended to
`tasks/v3-supervision-plan-2026-10-28.md`. `tasks/v3-acceptance.md` and `README.md` are not edited;
row text is proposed in the final report. Nothing is committed.

## Deviations and incidents (appended after the checks ran; the plan above is unchanged)

2026-11-01, during execution:

1. A5 (web single-genome upload route): the first attempt drove the portal with a Python that has no FastAPI and
   never submitted a job (harness defect; `a5/failed-attempt-no-fastapi`). The harness gained a `--web-python`
   option and the check was rerun with the portal's test environment. No route behaviour involved.
2. A6 (cohort routes): the first attempt failed identically on all three routes at the cohort pass, because the
   Price cohort metadata (built from an older layout) does not declare protein AXT23839.1 of GCA_003443575.1 as
   unplaced while the layout bundle's metadata does. The input builder now merges the layout bundle's
   `unplacedProteins`/`contigs`/`assemblyStatus` into the cohort metadata (`--layout-metadata`); the failed attempt
   is kept (`a6-ws/failed-attempt-metadata`). The three routes agreed on the failed state.
3. A6 (web cohort job): the first job integrated all three genomes identically to the CLI batch and then failed at
   the browser step (renderer expects `gspa_out/provenance.json`, the batch route writes
   `integrated.tsv.provenance.json`). This is a service defect in the new code, fixed in `service/cohort.py`
   (`render_run_dir`) with a regression test; the web route was rerun on ws, and the corrected code was staged to
   Shaheen before its job started. The failed job id is kept (`a6-ws/web/job.id.failed-attempt-render`).
4. Comparisons: the stage ledger's `feature_export` stage carries the route's own export directory at the stage
   level (`directory`); the comparison scripts exclude that one field (as they already excluded `elapsed_seconds`).
   Every other ledger field is compared.
5. A7: the first render batch failed one check per case on the harness's own expectation text (thousands
   separator in "rescored 25,838"); the expectation was corrected and the batch rerun; the failed attempt is kept
   (`a7/batch-a2-failed-attempt-expectation-format`). No rendering defect.
6. A6 informational cross-machine comparison and the Shaheen run: recorded in the report; the ws run of the three
   routes (run first, while the Shaheen job was staged) is the one the web-route reporting case (A7) uses.
