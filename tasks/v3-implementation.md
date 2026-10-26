# GSPA v3.0.0 implementation and release evidence

Started 2026-10-26 from `7714377`; development branch `development/v3.0.0`.
The authoritative scope is all twelve items in `tasks/todo.md`. The release is
authorized after implementation, biological validation, and full benchmark reporting.
No item is complete merely because its scaffolding or synthetic test passes.

## Acceptance and evidence

| Item | Required evidence before completion | Status |
| --- | --- | --- |
| 1 | Production-engine benchmark runner; frozen splits, input hashes, matched baselines and all-case reports | In progress: production runner and 73 tests including real CLI smoke; full datasets pending |
| 2 | Equivalent single/batch/web/Nextflow execution; stage status; coherent gene-keyed predictions | In progress: shared workflow used by integrate/annotate/web/NF, route tests and failure ledgers; full external-tool/Nextflow execution still pending |
| 3 | Correct tool contracts and end-to-end biological parser/dispatch/failure tests | In progress: adapters repaired with typed statuses, strict output contracts and biological dispatch tests; external-tool fixtures still need validation |
| 4 | Trained context model and held-out calibration/precision evidence; revisable inference and null assignments | In progress: learned stage and feature export wired into the shared workflow, contract /2, logistic/tree training and JVM inference; real-data fitting and held-out evaluation pending |
| 5 | Reaction/interaction context wired into production candidates and attribution; biological network ablations | In progress: direct EC/reaction and interaction context bridge implemented; real network ablations pending |
| 6 | Cohort recurrence with evolutionary/missing-data controls; no singleton certainty bug; independent-clade evaluation | In progress: joint two-pass batch inference and weighted held-out features wired; empirical clade evaluation pending |
| 7 | Weak-candidate retention; role-specific enzyme/complex attribution; same-operon decoy evaluation | In progress: raw weak candidates and distinct roles retained in learned/cohort paths; independent locus benchmark pending |
| 8 | Explicit media/biomass/assembly assumptions and alternative absence hypotheses; phenotype validation | In progress: 61 measured auxotroph cases completed under 3 registered protocols; regulatory closures improve binary agreement only 25/44 to 26/44; more organisms and function masking pending |
| 9 | Executed community models with transport/exchanges; established coculture/auxotroph validation | In progress: coupled SBML community LP and 26 mechanistic tests pass; empirical validation and workflow wiring pending |
| 10 | Auditable change ledger and truthful UI/statistics reconciled to predictions | In progress: browser/portal ledger implemented, EC/GO ancestry classification, counts and downloads tested in Chromium; full empirical UI reconciliation pending |
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

- Root branch `development/v3.0.0`; core workflow source through `3c75a38`; community benchmark through `0501c14`. All twelve
  acceptance gates remain open. Real-data learned-model fitting is still pending.
- Claude Sonnet 5.5 `/tmp/gspa-v3-context-workflow`: phase A (`d2e5ef7`) merged
  as `3c75a38`; phase B (`c924448`) awaits review fixes before merge. Required fixes include matched-medium
  knockout comparisons, reaction-specific EC bindings, weak/unmapped candidate retention,
  evidence-based shortlisting and private fresh subprocess outputs. Active review session
  `24084`, prompt/logs `context-review.txt` / `.jsonl`; original session `7544` completed.
- Community runner and review fixes merged as `4ff862f` / `0501c14`: all-case reporting,
  matched no-transfer controls, strain resampling, explicit failed/unknown coverage,
  model-owned cache. Job 9741 runs all 1,035 pairs under six oxygen/trace conditions
  in `/data/hohndor/gspa-v3/community-development-0501c14`.
- Claude Sonnet 5.5 `/tmp/gspa-v3-ledger-learned` is joining actual learned-model
  role/null artifacts to the UI and preserving stage order and scoring-mode labels.
  Session `72970`, logs `ledger-learned.jsonl`.
- Ontology concurrency agent and two reviews completed; merged as `45c9622`,
  `d9c5ea1`, `7d36adb`. Single-flight exact inferred has_part extraction, immutable
  results and source/revision-aware cache validation have biological concurrency tests.
- Panel13 job 9732 ended TIMEOUT after 30 minutes; stage0 outputs are preserved
  under `/data/hohndor/gspa-v3/panel13-dedup-6e446bb`. Thread dumps identified repeated
  full ontology extraction by concurrent genomes. Separate job 9739 uses the
  repaired `3c75a38` jar, identical deduplicated claims, persistent reasoner cache,
  2-hour allocation and 2,100-second per-arm timeout. It is a development regression.
- Metabolic job 9736 completed all 183 case/protocol combinations. Reports committed
  under `benchmark/v3/results/development-e8bcc07`; no large accuracy gain.
  Job 9738 independently replays all matched control/deletion LPs with COBRApy/GLPK.
- All 132 original full-proteome Price predictor files (23.6 GB) copied from IBEX
  to `/data/hohndor/gspa-v3/price-raw-sources`; all 132 destination hashes and sizes verified in job 9740.
  Matrix schema/FASTA ID probes are queued as job 9742 after the COBRA environment setup.
  Transfer manifest/log: `/tmp/gspa-v3-agents/price-transfer.jsonl` and `.log`.
- All 22 exact assembly versions acquired via NCBI Datasets with protein FASTA/GFF3;
  package MD5 and returned assembly IDs checked. Local sources and acquisition log:
  `/tmp/gspa-v3-price-genomes`. Functional text in GFF is not independent predictor
  evidence; context extraction must use coordinates/identifiers only.
- The Qwen endpoint is configured but the required empirical Qwen-versus-simple-method
  benchmark has not yet run. AllTheBacteria remains excluded.

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
- Price development runs use full nonzero target candidates. All seven variants completed with identical EC top-1
  before/after current pathway priors: CLEAN 60/148, DPZ 49/148, EnzBERT 53/148,
  MAPred 50/148, GraphEC 18/148, TopEC 5/148, pool4 60/148. Raw-ranking checks
  and input/output hashes are verified; reports are committed. These are reused-data, limited-context regressions.
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

- `3c75a38`: complete JVM suite and fat-jar build passed; Python learned trainer
  34 tests passed with feature contract /2. No production-trained weights yet.
- `b37e553` / `aeb602d`: ledger browser/portal tests passed (36 Python tests);
  three opt-in Chromium tests passed, including removed enzymes, filtered downloads,
  non-enzyme changes and unavailable ledgers. No JavaScript page errors.
- `e8bcc07`: 62 metabolic mechanism/runner tests and 187 benchmark tests passed
  (2 opt-in skips). Active regulatory closures are recorded separately from gene
  deletion; missing sidecars or malformed constraints fail as unknown.
- `dae00b3`: independent COBRApy/GLPK replay checker; two biological/hash-integrity
  tests passed locally. This validates implementation agreement, not biological truth.

- `0501c14`: 221 benchmark tests passed (2 opt-in skips), plus 29 community
  mechanism tests. Six-arm real coculture execution is running, not yet a reported result.
- `b2bbcf1`: complete Price source-transfer verification and 22 exact NCBI assembly
  acquisition records committed; see `docs/benchmarks/v3-price-context-data.md`.
- Latest complete JVM run: core 561 (2 skipped), CLI 66, predictors 241;
  0 failures/errors. This establishes implementation checks, not release readiness.
