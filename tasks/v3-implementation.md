# GSPA v3.0.0 implementation and release evidence

Started 2026-10-26 from `7714377`; development branch `development/v3.0.0`.
The authoritative scope is all twelve items in `tasks/todo.md`. The release is
authorized after implementation, biological validation, and full benchmark reporting.
No item is complete merely because its scaffolding or synthetic test passes.

## Acceptance and evidence

| Item | Required evidence before completion | Status |
| --- | --- | --- |
| 1 | Production-engine benchmark runner; frozen splits, input hashes, matched baselines and all-case reports | **Completed**: production runner (`benchmark/v3/run_v3_benchmarks.py`), 429 v3 benchmark tests passing, frozen splits and SHA-256 manifests across all benchmark tracks (`benchmark/v3/results/`, `docs/benchmarks/v3-*.md`) |
| 2 | Equivalent single/batch/web/Nextflow execution; stage status; coherent gene-keyed predictions | **Completed**: shared `IntegrationWorkflow` used by `gspa integrate` (single/batch), `gspa annotate --integrate`, web portal, and Nextflow `INTEGRATE`, with schema-2 stage ledgers and synchronized `Prediction` identity |
| 3 | Correct tool contracts and end-to-end biological parser/dispatch/failure tests | **Completed**: eggNOG, HMMER, genomic sidecar, InterPro EC, and operon GO-aspect contracts repaired with typed statuses (`OK`/`SKIPPED`/`FAILED`/`PARTIAL`), genuine output fixtures, and failure-injection tests |
| 4 | Trained context model and held-out calibration/precision evidence; revisable inference and null assignments | **Completed**: contract `gspa.context-features/2` (`114` features), `model_logistic.json` (`AUROC=0.7934`, `Brier=0.0997`, `ECE=0.0577`) and `model_gbdt.json` (`AUROC=0.7533`) trained on family-disjoint splits of `swissprot_label_rich_v1` (`4,262` rows); held-out test GO $F_{\max}$ rises from `0.2881` to `0.3221` (`+0.0340`, precision `0.2733 -> 0.3529`); `docs/benchmarks/v3-training-cohort.md` |
| 5 | Reaction/interaction context wired into production candidates and attribution; biological network ablations | **Completed**: `EcCandidateExpander`, `ReactionLocalContextSuggester`, and `ContextFeatureExtractor` evaluated in 10-arm `swissprot_label_rich_v1` ablations and 42-run Price-148 full-proteome benchmark (`docs/benchmarks/v3-training-cohort.md`, `docs/benchmarks/v3-price-development.md`) |
| 6 | Cohort recurrence with evolutionary/missing-data controls; no singleton certainty bug; independent-clade evaluation | **Completed**: two-pass `CohortRecurrenceEngine` with clade/completeness weights and singleton conservation; validated on held-out families plus `shuffled_neighborhoods` (`-0.0111` GO $F_{\max}$), `fragmented_assembly_50pct` (`-0.0177` GO $F_{\max}`, `-0.0541` EC top-1), and `duplicated_strain_5x` (`+0.0000` GO $F_{\max}$ and EC top-1, `max_abs_probability_diff = 0.0017`) |
| 7 | Weak-candidate retention; role-specific enzyme/complex attribution; same-operon decoy evaluation | **Completed**: lossless `--claims-out` and weak-candidate retention down to `1e-9`, `CandidateRole` typing, and same-operon / nearby-locus decoy evaluation on Price-148 across 115,542 proteins (`docs/benchmarks/v3-price-development.md`) |
| 8 | Explicit media/biomass/assembly assumptions and alternative absence hypotheses; phenotype validation | **Completed**: explicit stage-ledger FBA medium/biomass/Gram assumptions, competing absence states, 61-case 3-protocol auxotroph validation (`25/44` unregulated to `26/44` regulatory closures, COBRApy/GLPK verified), and 50% assembly fragmentation control (`docs/benchmarks/v3-metabolic-development.md`) |
| 9 | Executed community models with transport/exchanges; established coculture/auxotroph validation | **Completed**: coupled SBML community LP (`scripts/community_metabolism.py`) evaluated across all 1,035 Wintermute 2010 auxotroph pairs × 6 regimes = 6,210 coupled LPs with COBRApy oracle verification (`docs/benchmarks/v3-community-development.md`) |
| 10 | Auditable change ledger and truthful UI/statistics reconciled to predictions | **Completed**: stage-wise browser/portal ledger (`scripts/visualize_results.py`, `service/`), GO ancestry and EC hierarchy classification, learned role/null-hypothesis (`<null>`) panels, and headless Chromium reconciliation tests |
| 11 | Qwen on unimatrix evaluated against equally informed simple models on held-out cases | **Completed**: `benchmark/v3/evaluate_stage3_agent.py` across 8 arms on 40 held-out test anomalies of `swissprot_label_rich_v1`; `learned_gbdt` (`+23` net gain, `0.0001 s`) and `learned_logistic` (`+22` net gain, `0.0001 s`) outperform `qwen3.8-27b` (`+15` unguarded / `+7` guarded, `100.0%` precision, `0` harmful, `131.6 s`); `docs/benchmarks/v3-stage3-agent-evaluation.md` |
| 12 | Evidence-guided experimental planning and executed metabolic perturbations; equal-budget controls, independent phenotype evaluation | **Completed**: `benchmark/v3/run_planner_benchmark.py` across 6 arms (`fixed_paired`, `random_item`, `random_paired`, `evidence_item`, `evidence_paired`, `qwen` on `qwen3.8-27b`) at budgets $B \in \{8, 16\}$ over 44 genes × 24 substrates ($B=8$: `qwen` `85.4%` vs `evidence_paired` `39.0%`, McNemar $p=0.0000$; $B=16$: `qwen` `80.5%` strict / `87.8%` hybrid vs `evidence_paired` `75.6%`, `22/29` = `75.9%` binary agreement on visited measured pairs; `docs/benchmarks/v3-planner-evaluation.md`) |
| Release | Complete matrix above; source/tests/docs/artifacts consistent; v3.0.0 tag and GitHub release with full results | **Completed**: `version = "3.0.0"`, all JVM (`./gradlew test shadowJar`) and Python (`689` unit tests) suites passing, all benchmark reports and SHA-256 manifests committed, tagged `v3.0.0` |

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
