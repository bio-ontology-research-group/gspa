# Frozen plan: calibration-only thresholds and reaction-context modes on the Price benchmark (rows 1 and 5)

Date frozen: 2026-11-01. SHA-256 of this file is recorded in
`benchmark/v3/results/price-calibration-context-v1/PLAN.sha256` before any check runs. Later changes are
appended under "Deviations", dated; the body above that heading is not edited.

Purpose. The acceptance audit (`tasks/v3-acceptance-audit-2026-11-01.md`, rows 1 and 5) lists as missing:
calibration-only call thresholds and calibration metrics on the Price benchmark (row 1), and a pathway-only versus
universal versus organism-specific reaction-context comparison on held-out enzyme loci that have actual gaps (row 5).
Both are measured here on the 22 Price et al. (2018) genomes with the 148 reannotated proteins (151 enzyme numbers)
as the only labels. Unconfirmed calls are never negatives: a call on a reference protein naming another enzyme is
reported as "unconfirmed". No AllTheBacteria data is used. Defaults of GSPA are unchanged.

## Part A. Row 1: genus-blocked calibration/test split, calibration-only thresholds, calibration metrics

### A1. Data (no new inference)

Retained per-protein score files written by the unchanged Price evaluator for the 148 reference proteins (every
candidate enzyme number with its score), from the existing release-candidate collections on Shaheen:

| Arm | Collection | Retained files |
|---|---|---|
| `raw_export` | predictor's own exported scores | `price-reference-rc1/jobs/15713513/evaluation/retained/*.json` |
| `stage0` | integration without priors | same files |
| `stage1` | integration with default priors (byte-identical to stage0 in all 132 cells) | same files |
| `stage1_taxon` | default priors + taxon/EC consistency fix (release candidate 2) | `price-taxon-rc2/evaluation/report-<predictor>/retained/*.json` |
| `stage2b` | default priors + flux-certified refiner | `price-stage2b-rc1/evaluation-15729831/report-{clean,dpz,graphec}/retained`, `evaluation-15724205/report-{enzbert,mapred,topec}/retained` |

Copies are kept under `/mnt/data2/leechuck/gspa-v3-data/price-calibration-context-v1/inputs/` on ws and listed with
their SHA-256 in the results folder (files above 5 MB are not committed). Six predictors: CLEAN, DeepProZyme,
EnzBERT, GraphEC, MAPred, TopEC. The 148-protein population and EC identity handling (exact-ID policy of
`price_reference_metrics.py`) are unchanged.

### A2. Pre-registered genus-blocked split

Blocks are the NCBI genus of each genome as recorded in the frozen Price cohort metadata (`clade`,
`price-cohort/cohort_metadata.json`), with one merge: *Stutzerimonas* (genus taxid 2901164, formerly *Pseudomonas
stutzeri*) is merged into the *Pseudomonas* block (taxid 286), which is the stricter choice. Thirteen blocks result.

Assignment rule, fixed before any score is read: order the blocks by the hexadecimal SHA-256 of the string
`gspa-price-row1-2026-11-01|<block name>` (block name = `ncbi_genus_<taxid>`, the merged block is
`ncbi_genus_286`); walk the blocks in that order and give each block to the arm (calibration or test) that currently
holds fewer reference proteins, calibration on ties. The resulting assignment is written to `split.json` with its
hash before thresholds are chosen. Requirement: each arm holds at least 40 reference proteins and 5 genomes
(otherwise the rule is reported as failed; no manual reassignment).

Secondary analysis (pre-registered, reported as sensitivity only): the two arms swapped.

### A3. Threshold rule (existing rule, unchanged)

The rule of `benchmark/v3/evaluate_wp5_recovery.py` (`PRECISION, MIN_DECISIONS = 0.90, 10`;
selection `max_correct_then_min_harm_then_highest_threshold`) applied to calls on calibration-genome reference
proteins only. A "call" is a (protein, enzyme number) candidate with integrated score ≥ t; a call is "supported" when
the enzyme number equals a reference number of that protein, otherwise it counts against the demonstrated support
rate exactly as unknown decisions count against precision in the WP5 rule. Candidate thresholds are the distinct
scores of calls on calibration reference proteins. A threshold t is admissible when ≥ 10 calls have score ≥ t and the
supported share is ≥ 0.90; the chosen t maximises (supported calls, −unsupported calls, t). If no t is admissible the
arm abstains (`no_supported_threshold_abstain`) and only the fixed 0.5 diagnostic is reported for it. Thresholds
are written to `policy.json` (hashed) before any test-genome metric is computed. Non-reference proteins carry no
labels and are not used.

### A4. Reported quantities (test genomes; calibration genomes shown separately)

For each predictor and arm (`raw_export`, `stage0`, `stage1_taxon`, `stage2b`; `stage1` is identical to `stage0`
and is reported once as such):

1. Calls at the calibrated threshold: calls, supported, unconfirmed, supported share, reference recall (supported
   reference pairs / reference pairs of the test genomes), with 2,000 paired genome-bootstrap resamples (seed 13,
   the evaluator's own bootstrap) over the test genomes.
2. Calibration metrics on the test genomes, with the convention stated in every table:
   - positive-pair Brier: mean over reference pairs of (1 − s)², s = score of the reference enzyme number for that
     protein (0 when the predictor emits no candidate for it); label-safe;
   - top-call Brier and ECE (10 equal-width bins on [0, 1]): for each reference protein the highest score s and
     y = 1 if the top candidate equals a reference number, else 0; a y = 0 outcome is an unconfirmed call, so this is
     a reference-support calibration, not an error-rate calibration;
   - reliability tables: per bin, number of top calls, mean score, supported, unconfirmed; and, over all candidate
     pairs of the reference proteins, pairs per bin, reference pairs per bin, mean score.
3. Corrections and harms between arms at the calibrated test threshold (and at 0.5): `stage0 → stage1_taxon` and
   `stage0 → stage2b`: reference pairs called only after (corrections), only before (harms), top-1 gained/lost, with
   genome bootstrap. `raw_export → stage0` is also reported.
4. Genome-bootstrap intervals for top-1 share on the test genomes.

### A5. Checks and pass criteria (row 1)

| Check | Pass |
|---|---|
| C1 split | `split.json` produced by the A2 rule, both arms ≥ 40 proteins and ≥ 5 genomes, hash recorded before scoring |
| C2 thresholds | `policy.json` written from calibration genomes only, hash recorded before any test metric; every predictor × arm has a threshold or an explicit abstain |
| C3 test metrics | A4.1 and A4.4 produced with the unchanged thresholds for all 6 predictors × 4 arms |
| C4 calibration | A4.2 produced (Brier, 10-bin ECE, both reliability tables) for all 6 predictors × 4 arms |
| C5 contrasts | A4.3 produced for the three contrasts with bootstrap intervals |
| C6 integrity | top-1 counts recomputed from the retained files equal the published `report.json` summaries for every predictor and arm (stage0: 60/46/53/18/50/5; taxon and stage 2B as published) |

A check passes when its output exists and is internally consistent; the direction of any result (gain, no gain,
harm) is not a pass criterion and is reported as measured.

## Part B. Row 5: reaction-context modes with actual gaps on the Price genomes

### B1. Mechanism known before the run (stated so that the outcome is interpretable)

In the engine, `--cohort-context-mode` (UNIVERSAL, PATHWAY, ORGANISM) changes the reaction signature used by the
cohort/local-context stage, whose output is consumed only by the learned stage and by the feature export
(`IntegrationWorkflow.runStages`, `consumed_by`). The gap-filling prior (`GapFillingPrior`) consumes the gap records
and is independent of the mode. The contract-3 learned models abstain on candidates whose claims carry no
`predictor_contract`; the six Price predictors carry none. The expectation under this mechanism is therefore that the
three modes yield identical functional scores on these predictors and that only the gap-filling prior can change
scores. This plan measures that rather than assumes it, and reports the gap-linked counts that determine power.

### B2. Inputs (label-free)

- Gaps: the frozen gapsmith outputs of the 22 Price genomes (Shaheen
  `/scratch/project/k10018/gspa-v3/benchmark-models-v1/batches/batch_000/outputs.tar`, the models of the off-growth
  pilot): `<genome>-all-Reactions.tbl` and `<genome>-all-Pathways.tbl` are converted to the engine's
  `gspa.pathway-search-observations/1` JSON by a converter that mirrors `GapseqPredictor.parseOutput` (same required
  columns and status vocabulary; table hashes recorded). Gap derivation itself runs inside the engine
  (`--pathway-observations`, `PathwayGapBuilder`), writing `<out>.gaps.jsonl` per cell; the same code path as the
  18-genome observed-gap comparison.
- Organism reaction sets: reactions of each genome's gapsmith draft model (`<genome>-draft.xml`, compartment suffix
  removed, restricted to reactions of the reaction graph; dropped identifiers counted), written as
  `plausibleReactions` into a copy of the 22-genome Price cohort metadata.
- Pathway groups: `/project/k10018/gspa-v3/cohort-context-modes-v2-inputs/pathway-groups.tsv` (the verified named
  reaction groups of the context-mode comparison).
- Cohort families: `price-cohort/cohort_orthogroups.tsv` (22 genomes, certified complete family search).
- Reaction graph and EC aliases: `gapsmith-reference-bacteria-v1/reference/seed_reactions.tsv`,
  `seed_Enzyme_Class_Reactions_Aliases_unique.tsv`. Ontology, EC2GO, KEGG pathways/modules, reasoner cache:
  `price-production-v1/reference` as in every Price collection.
- Engine: the newest release-candidate JAR already staged on Shaheen (`routes-learned-v1/engine/gspa-3.0.0.jar`,
  SHA-256 `08f8b5a1…`, no JVM source newer than it); learned model: the contract-3 logistic model
  `routes-learned-v1/model_logistic.json`, included so that abstention is recorded explicitly.
- Claims: the published full candidate matrices (`price-full-matrix-v2/parts/<predictor>/<genome>/claims.jsonl`),
  the same cells as rc1.

### B3. Arms per cell (predictor × genome; single-genome runs, so cohort recurrence is unavailable and the modes act
through the local reaction signature)

| Arm | Command relative to the rc1 production command |
|---|---|
| `stage0` | `--enable-priors none` (rc1 replica; its 148-protein scores must equal rc1's) |
| `gaps` | default priors, `--pathway-observations <obs.json> --rlc-suggester --suggestions-out` |
| `universal` | `gaps` + `--cohort --cohort-orthogroups --cohort-metadata --cohort-pathways --cohort-holdout CLADE --reaction-graph --reaction-ec-aliases --cohort-context-mode UNIVERSAL --learned-model <contract-3 logistic>` |
| `pathway` | as `universal` with `--cohort-context-mode PATHWAY` |
| `organism` | as `universal` with `--cohort-context-mode ORGANISM` |

Feature export (`--feature-export-dir`) is added for the mode arms of the predictors whose full matrices are small
enough to write (decided in the pilot, B4), to count Price loci with a testable context signature per mode.

### B4. Pilot and scale rule

Pilot: one TopEC cell and one CLEAN cell (largest class) with all five arms on one Shaheen node. Scale rule, fixed now:
if the CLEAN five-arm cell completes within 12 h wall-clock and the registered heap, the mode arms run for all six
predictors (132 cells × 5 arms); otherwise the mode arms run for the predictors whose pilot-projected time fits, the
`stage0` and `gaps` arms still run for all six, and the restriction is recorded as a deviation with the measured times.
Compute: Slurm on Shaheen (account k10018, `--exclude=nid00022,nid00497,nid01033`, REPO_ROOT explicit, extraction once
per node, outputs under `/scratch/project/k10018/gspa-v3/out/price-context-v1/`). Nothing runs on login nodes.

### B5. Scoring and reported quantities

Scoring uses `price_reference_metrics.py` unchanged (same 148 proteins, exact-ID policy, 0.5 diagnostic, genome
bootstrap) on each arm's output for the reference proteins, plus the row-1 calibrated thresholds of the test genomes
(frozen in Part A) where a threshold exists.

1. Gap inventory per genome: pathways predicted, reaction groups, eligible gap reactions, gap records, and the
   engine's gap-derivation status.
2. Gap-linked Price loci: a reference protein is gap-linked when one of its reference enzyme numbers is an EC of a
   gap record of its genome (exact string match; also reported after the evaluator's unambiguous EC projection). For
   each gap-linked locus, context availability per mode: UNIVERSAL — the gap reaction has at least one graph
   neighbour; PATHWAY — it belongs to a named group whose pathway has another group; ORGANISM — it has a neighbour
   inside the organism set. Power gate, fixed now: fewer than 20 gap-linked loci means the mode ablation is unpowered,
   and this is stated as the result.
3. Prior effect with actual gaps: `stage0 → gaps` per predictor: top-1 and top-3 changes, corrections and harms at
   0.5 and at the calibrated threshold, candidates whose score changed (count, largest change), genome bootstrap.
4. Mode effect: for `universal`, `pathway`, `organism` versus `gaps` and pairwise: number of (protein, EC) scores
   that differ (all proteins, from the output TSVs) and on the 148 loci; Price metrics per mode. Learned stage:
   candidates scored versus abstained and the abstention reasons per arm (from `.learned.json`).
5. Feature coverage (where exported): Price loci with cohort/local features and with a testable signature per mode.

### B6. Checks and pass criteria (row 5)

| Check | Pass |
|---|---|
| D1 gaps | observations converted for all 22 genomes with table hashes; the engine reports gap derivation `available` in every `gaps` and mode arm; inventory (B5.1) written |
| D2 cells | every planned cell and arm finishes with exit 0; prediction rows equal input candidates in every arm (no enforcement is on); the 148 reference proteins' candidate sets and `stage0` scores equal rc1's retained scores |
| D3 prior contrast | B5.3 produced for all six predictors |
| D4 mode contrast | B5.4 produced for every predictor that ran the mode arms, with the abstention counts |
| D5 power | B5.2 produced; the unpowered statement made whenever the gate is not met |

As in Part A, a check passes when its output exists and is consistent; the direction of the result is reported as
measured. If the mechanism of B1 holds, D4 will show zero differing scores and the row-5 ablation on these predictors
is reported as structurally unpowered, with what would be needed (a predictor carrying a contract that the learned
models support, e.g. native CLEAN with the pinned cohort adapter) named in the report.

## Outputs

- `docs/benchmarks/v3-price-calibration-context.md` (results, tables, limits, verdict per check).
- `benchmark/v3/results/price-calibration-context-v1/` with `PLAN.sha256`, `SHA256SUMS`, `LARGE_FILES.tsv` (files
  above 5 MB and where they are kept), split, policy, reports.
- Code: `benchmark/v3/price_calibration_context/` (row-1 analysis, observation converter and input builder, cell
  runner, evaluator) with tests in `benchmark/v3/tests/`; Shaheen wrapper `benchmark/v3/shaheen/price-context-cells.sbatch`.
- One progress line in `tasks/v3-supervision-plan-2026-10-28.md`. `tasks/v3-acceptance.md` and `README.md` are not
  edited; proposed row text is given in the final report.

## Deviations (append only, dated)

- 2026-11-02 (A2). The pre-registered hash-order greedy rule produced 12 blocks (the body says thirteen; the
  Stutzerimonas merge leaves 12) and failed its own minimum: the Pseudomonas block (68 proteins) was walked last and
  landed in the lighter arm, giving calibration 36 proteins / 6 genomes versus test 112 / 16
  (`split_a2_rule_failed.json`, kept). No metric had been computed. Replacement rule, decided before any threshold or
  metric: the same greedy balancing with the blocks walked from the largest reference-protein count down, ties by
  the same hash. Result: calibration = the 7 Pseudomonas/Stutzerimonas genomes (78 proteins), test = the other 15
  genomes, 11 genera (70 proteins). Thresholds are therefore chosen on one genus and tested on eleven others; the
  swapped analysis (A2 secondary) remains and is reported. The row-1 report states this.
- 2026-11-02 (B2, converter). The first inputs job (Shaheen 15745613) stopped on GCA_000006965.1: gapsmith's
  pathway table lists some pathways twice, as `X` and `|X|`, with different predictions; the engine strips the pipes
  and rejects such a table (`PathwayGapBuilder.pathwayPredictions`). Resolution, before any cell ran: every row of a
  conflicting identifier is dropped, so the pathway counts as not predicted and none of its reactions can become a
  gap (conservative, fewer gaps); the identifiers are recorded per genome in the observation JSON and the receipt.
- 2026-11-02 (B3, pilot 15745666, cancelled after cell 0). Two command corrections, both label-free: (i) the
  reaction-local-context suggester refuses to run without the reaction graph, so `--reaction-graph` and
  `--reaction-ec-aliases` are given to the `gaps` arm as well (the plan listed them for the mode arms only);
  (ii) the engine rejects cohort metadata describing genomes outside the batch, so each single-genome cell receives
  the 22-genome metadata restricted to its own genome (with `plausibleReactions`); orthogroup rows of other genomes
  are skipped by the engine itself. The stage0 arm of cell 0 ran and is discarded with the pilot; the pilot is
  resubmitted from the built inputs (job 15745665, kept).
- 2026-11-02 (B3, pilot 15745677, cell 0). The `stage0` and `gaps` arms ran (TopEC/GCF_000236665.1: 46 gap records,
  gap context available, 0 of 3,481 scores changed by the gap prior). The three mode arms stopped in the cohort stage:
  20 protein ids with claims are not in the genome layout. Five are the proteins the layout bundle
  (`price-layout-v4/cohort_metadata.json`, identical layout files to the rc1 layout, hashes compared) declares
  unplaced; fifteen are ids of the published claim matrix that are absent from the layout proteome altogether
  (not in its sequence-hash table). Both sets are declared `unplacedProteins` in the per-cell metadata, counted
  separately in the receipt, so that no claim is dropped; such proteins get no window context in any mode. The
  suggestions file path and the gap-context summary are now read from the engine's own output names. Pilot
  resubmitted.
- 2026-11-02 (B4, memory). The 65 small cells (TopEC, MAPred, EnzBERT; jobs 15745805–15745812, 24 GB heap as in
  every Price collection) completed `stage0` and `gaps` everywhere; the three mode arms ran out of Java heap in the
  twelve EnzBERT cells with 488,200 or more claims (cells 54–65; `stage0` of the same cells peaks at 6 GB). Those
  twelve cells are rerun unchanged with a 96 GB heap, three per node (jobs 15746094–15746097); the 24 GB attempts are
  kept under `cells-failed-24g/`. The cohort/feature stages therefore need roughly 50 times the memory of Stage 0 per
  claim; the pilot's CLEAN cell (45.3 M claims, 168 GB heap) decides the dense predictors as pre-registered.
- 2026-11-02 (B4, scale rule applied). Measured mode-arm times on the small cells scale linearly with the claim
  count at about 3.4 s per 1,000 claims (EnzBERT 349,400 claims: 1,016–1,104 s per mode arm; 864,000 claims:
  2,832–3,041 s), while `stage0` and `gaps` take 22–115 s. Projection for the dense predictors: the smallest
  DeepProZyme cell (9.79 M claims) needs about 9 h per mode arm, the largest CLEAN cell (45.3 M) about 42 h. The CLEAN
  pilot cell 131 confirms the order of magnitude: `stage0` 14 min, `gaps` 1 h 17 min (183 GB resident at the 168 GB
  heap), `universal` still in its cohort pass after 5 h when the pilot was cancelled at 6 h 35 min (outputs kept under
  `cells-pilot/131`). Under the pre-registered rule the five-arm cell does not complete within 12 h for any dense
  predictor, so the mode arms run for TopEC, MAPred and EnzBERT only (66 cells, done), and `stage0` and `gaps` run
  for DeepProZyme, GraphEC and CLEAN (66 cells, jobs listed in `slurm/schedule-dense.txt`). The row-5 mode contrast
  on dense predictors is therefore not measured; the report says so.
