# GSPA v3.0.0 implementation and release evidence

Started 2026-10-26 from `7714377`; development branch `development/v3.0.0`.
The authoritative scope is all twelve items in `tasks/todo.md`. The release is
authorized after implementation, biological validation, and full benchmark reporting.
No item is complete merely because its scaffolding or synthetic test passes.

## Acceptance and evidence

| Item | Required evidence before completion | Status |
| --- | --- | --- |
| 1 | Production-engine benchmark runner; frozen splits, input hashes, matched baselines and all-case reports | In progress: production runner and 73 tests including real CLI smoke; full datasets pending |
| 2 | Equivalent single/batch/web/Nextflow execution; stage status; coherent gene-keyed predictions | In progress: synchronized prediction state, identical single/batch stages, per-stage change ledger; annotate/web/NF parity pending |
| 3 | Correct tool contracts and end-to-end biological parser/dispatch/failure tests | In progress: adapters repaired with typed statuses, strict output contracts and biological dispatch tests; external-tool fixtures still need validation |
| 4 | Trained context model and held-out calibration/precision evidence; revisable inference and null assignments | In progress: typed feature extraction, logistic/tree training and JVM inference with explicit null layer; synthetic agreement tests pass; real training and production wiring pending |
| 5 | Reaction/interaction context wired into production candidates and attribution; biological network ablations | Pending |
| 6 | Cohort recurrence with evolutionary/missing-data controls; no singleton certainty bug; independent-clade evaluation | Pending |
| 7 | Weak-candidate retention; role-specific enzyme/complex attribution; same-operon decoy evaluation | Pending |
| 8 | Explicit media/biomass/assembly assumptions and alternative absence hypotheses; phenotype validation | Pending |
| 9 | Executed community models with transport/exchanges; established coculture/auxotroph validation | In progress: coupled SBML community LP and 26 mechanistic tests pass; empirical validation and workflow wiring pending |
| 10 | Auditable change ledger and truthful UI/statistics reconciled to predictions | In progress: integration ledger implemented; web presentation pending |
| 11 | Qwen on unimatrix evaluated against equally informed simple models on held-out cases | Pending |
| 12 | Evidence-guided experimental planning and executed metabolic perturbations; equal-budget controls, independent phenotype evaluation | In progress: SBML experiment engine, fixed/random/evidence/Qwen policies, 40 mechanistic tests pass; workflow and learned attribution pending |
| Release | Complete matrix above; source/tests/docs/artifacts consistent; v3.0.0 tag and GitHub release with full results | Pending |

## Benchmark scope

Use established, traceable datasets; **AllTheBacteria is excluded**. Required
tracks cover exact EC attribution (Price with accession correction), prospective
and homology-controlled GO, reference genomes across domains, operon and complex
role discrimination, taxonomic/anatomical constraints, reaction masking, metabolic
growth/knockouts, BacDive phenotypes, incomplete genomes/communities, secondary
metabolism, and agent/condition-planning comparisons. Each track must disclose
experimental versus synthetic truth and missing coverage. Dataset preparation,
model training, calibration and final test labels are separate. Benchmark
execution must use production inference, not parallel scoring formulas.

## Current work and process handles

- Root: core prediction-state repair and integration; then context/learning and
  execution parity. Local Gradle tests are permitted; large jobs use compute nodes.
- Claude Code Sonnet 5.5 worktrees `/tmp/gspa-v3-predictors`,
  `/tmp/gspa-v3-benchmarks`, and `/tmp/gspa-v3-experiments`: predictor contracts,
  production benchmark harness, and executed metabolic experiment planning.
  Prompts and streaming logs are under `/tmp/gspa-v3-agents/`.
- agy attempts ended without edits because headless auto mode was unavailable.
  Claude Code with explicit file/test permissions is running successfully.
- Unimatrix SSH verified; `/data/hohndor/gspa` exists; no hohndor Slurm jobs at
  initial inspection. Qwen endpoint is configured in the existing agent sidecar.

## Validation log

- Audit baseline: 488 JVM tests passed, 2 skipped; seven targeted Python tests
  passed. These establish the starting suite only, not v3 correctness.
- Existing untracked research data and paper build files are retained unchanged.
- Full JVM suite passed after prediction identity synchronization and recurrence
  mass-preservation fixes. Biological regressions cover competing enzyme loci,
  evidence retraction, negative recurrence, and retained alternative candidates.
  These are implementation checks, not held-out benchmark results.

## Development evidence (not release benchmark results)

- `b15c3df`: condition search continues after a non-growing/infeasible reference;
  explicit maintenance-only product experiments support non-growing cells;
  redundant isozyme deletions remain executable. Forty-seven mechanistic tests
  pass. The phenotype runner has nine tests, including real LPs and truth-blind
  simulation checks. These mechanisms are not an annotation-accuracy result.
- `e63e735`: learned feature/training/inference foundations. Full JVM tests and
  31 Python training tests pass. No real-data model has been fitted or enabled.
- Empirical auxotroph benchmark: all 61 published cases ran at `d4e9895` in
  Slurm job 9727, with 25/44 binary observations agreeing; 17 partial-growth
  cases are separate. This is unregulated FBA using repository/model defaults,
  not the paper's full regulatory configuration. See
  `docs/benchmarks/v3-metabolic-development.md` and its linked all-case report.
- Price development runs use full nonzero target candidates. Initial completed
  DeepProZyme, TopEC and pool4 runs show identical EC top-1 before/after current
  pathway priors (49/148, 5/148, 60/148 respectively). Remaining runs and raw-rank
  checks are pending. These are reused-data, limited-context regressions.
- Panel13 job 9726 failed its input audit before inference: 420 exact duplicate
  records remain in the historical merged claims (distinct homolog hits are now
  allowed). A separate, hashed deduplication preparation is required; originals
  are retained. Job 9716 on node008 was cancelled after no output; job 9726
  was the separate node002 retry.

- `822ab69`: gene identity, single/batch stage parity, ledger, adapter foundation,
  production benchmark harness, metabolic experiment mechanics. The built jar
  still reports 2.2.2; its source revision and SHA-256 distinguish this development
  artifact from released software. Version changes wait for release readiness.
- Unimatrix job 9714 ran the panel13 input audit under Slurm and terminated before
  inference: the validator incorrectly treated distinct DIAMOND targets from one
  predictor as duplicate records. Full claim identity fixes that check. The
  production combiner already collapses correlation groups; rescue support counts
  are also being corrected to count predictor channels, not repeated hits.
- Remote run directory: `/data/hohndor/gspa-v3/development-822ab69/`. This reused
  panel is a development regression, explicitly not a new blind evaluation.
- Local metabolic test interpreter: `/tmp/gspa-v3-testenv/bin/python` with
  numpy 2.2.6, scipy 1.15.3, python-libsbml 5.21.2. Forty experiment tests and
  twenty-six community tests pass. They verify mechanisms on synthetic SBML
  networks, not accuracy on measured organisms.
- Candidate empirical metabolic sources verified against the primary papers:
  iML1515/KEIO 16-carbon-source growth screens (DOI 10.1038/nbt.3956;
  https://pmc.ncbi.nlm.nih.gov/articles/PMC6521705/; supplementary datasets 1, 10,
  11) and the Harcombe E. coli/S. enterica/M. extorquens consortia
  (DOI 10.1016/j.celrep.2014.03.070;
  https://pmc.ncbi.nlm.nih.gov/articles/PMC4097880/). Raw files still need
  acquisition and the experimental protocol must distinguish measured data from
  model-imposed export assumptions. Model-fitting phenotypes cannot serve as an
  independent function-recovery test without masking and split controls.
