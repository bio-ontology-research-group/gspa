# Frozen plan: contract-bearing claims on the 22 Price genomes, learned stage on Price (rows 1, 4, 5, 11)

Date frozen: 2026-11-02. The SHA-256 of this file is recorded in
`benchmark/v3/results/price-contract-learned-v1/PLAN.sha256` before any check runs. Later changes are appended
under "Deviations", dated; the body above that heading is not edited.

## Purpose

In rows 1, 5 and 11 of the acceptance record the learned stage abstained on Price because the six Price enzyme
predictors' claims carry no predictor contract. The learned models (contract 3) were trained on one predictor only,
the genus-fenced Swiss-Prot DIAMOND transfer (`swissprot_diamond_baseline`). This plan produces genuine claims from
that predictor, with its contract, for all 22 Price genomes, and measures on them:

- rows 1 and 4: the learned stage on the Price benchmark (scored / abstained counts, calibration, corrections and
  harms, supported precision) under the unchanged Price evaluator and the frozen row-1 genus-blocked split and rule;
- row 5: the universal / pathway-only / organism-specific context modes with the learned stage executing;
- row 11: the production anomaly extractor on GO-bearing outputs, and the truth intersection;
- a mandatory leakage check between the Price-148 proteins and the DIAMOND reference.

Labels: the 148 reannotated Price et al. (2018) proteins (151 enzyme numbers; manifest
`benchmark/manifests/suite6_price149_manifest.tsv`, the replicon row excluded as in every Price analysis), exact-ID EC
policy of `price_reference_metrics.py`, genome bootstrap 2,000 resamples seed 13. A call naming another enzyme on a
reference protein is "unconfirmed", never a negative. No AllTheBacteria data. GSPA defaults unchanged. No JVM
source change; the worktree JAR is not rebuilt.

## Fixed components (identities to be verified at run time)

| Component | Identity |
|---|---|
| DIAMOND (ws `/usr/bin/diamond`) | 2.1.11, SHA-256 `d6bea1f3…` |
| Reference (training baseline) | `~/.cache/gspa-v3-data/training-cohort/baseline/reference.faa` `1d8c4ebb…`, `reference_labels.tsv` `db8a8ac9…`; contract reference digest `2ef431fe…` |
| Adapter | `benchmark/v3/homology_baseline.py` `f84f9413…`; contract builder `benchmark/v3/replay_predictor_contract.py` |
| Parameters | threads 4, 10 hits, e-value 1e-5, sensitive (parameter digest `f8a5d09a…`) |
| Driver | `benchmark/v3/routes_closure/predict_contract_claims.py` (unchanged) |
| Proteomes | `~/.cache/gspa-v3-data/price-layout-v4/all_proteins.faa` `d84a1da6…` (the Price layout proteome) |
| Engine | pinned copy of the release-candidate JAR `08f8b5a1…` (Shaheen `routes-learned-v1/engine/gspa-3.0.0.jar`) |
| Learned models | contract-3 `benchmark/v3/results/predictor-support-v1/models/default/model_logistic.json` (`c9df93d6…`, primary) and `model_gbdt.json` (secondary) |
| Row-1 split | `benchmark/v3/results/price-calibration-context-v1/row1/split.json` (`950209…`; calibration = 7 Pseudomonas/Stutzerimonas genomes, test = 15 genomes in 11 genera), primary and swapped |
| Row-5 inputs | `price-context-v1/inputs` (observations, 22-genome cohort metadata with `plausibleReactions`), `price-cohort/cohort_orthogroups.tsv`, `layout_metadata.json`, pathway groups, reaction graph and EC aliases, Price reference (go.owl, ec2go, KEGG) — the same files as the row-5 cells |

## Step L. Leakage check (mandatory; run before any score is read)

L1. Price-148 sequences (from the layout proteome) are searched against `reference.faa` with the same DIAMOND binary
(`--sensitive`, e-value 1e-5, unlimited targets). A reference entry is *near-identical* to a Price-148 protein when
identity ≥ 95 % with query coverage ≥ 80 % and subject coverage ≥ 80 %; *exact* when the sequences are identical
(SHA-256 of the uppercase residue string). Reported: number of Price-148 proteins with an exact or near-identical
reference entry, the reference accessions, their organisms and NCBI taxids (from the Swiss-Prot flat file
`~/deepgoplus/data/uniprot_sprot.dat.gz`, `a6907c4d…`), and whether that entry carries the protein's Price EC.
L2. Same-organism reference entries: the number of reference entries whose Swiss-Prot `OX` taxid is one of the 22
Price organism taxids, and the genus-fence gaps that let them through (the fence uses the first token of the organism
name).
L3. Sensitivity set: every Price-148 protein with an exact or near-identical reference entry is excluded; all row-1/4
metrics are reported on the full population and on the reduced population. L3 also flags the Swiss-Prot
experimental Price labels (`benchmark/v3/offgrowth/data/swissprot_experimental_ecs_price22.tsv`) whose accession is
itself in the reference (circular truth for row 11).

## Step C. Claims (ws)

`predict_contract_claims.py` on all 22 genomes (one run, `--threads 4`). Pass C: 22/22 genomes with claims, every
contract identity equal to one in the model's input-support profiles, input hashes unchanged at the end.
Per-genome proteins, proteins with hits, claims (GO / EC) tabulated.

## Step I. Integration cells (Shaheen, Slurm)

One cell per genome (22), single-genome `gspa integrate --batch-manifest`, the row-5 command family of
`run_context_cell.py` (taxon from the layout report, essential profile bacteria, Price reference files), arms:

| Arm | Command |
|---|---|
| `stage0` | `--enable-priors none` |
| `stage1` | default priors |
| `learned` | default priors + `--learned-model` logistic |
| `learned_gbdt` | default priors + `--learned-model` gbdt |
| `gaps` | default priors + `--pathway-observations` + `--rlc-suggester` + reaction graph/aliases (row-5 `gaps` arm) |
| `gaps_learned` | `gaps` + `--learned-model` logistic (no cohort) |
| `universal`, `pathway`, `organism` | `gaps` + cohort inputs + `--cohort-context-mode` + `--learned-model` logistic (row-5 mode arms, unchanged) |

Retained per arm: EC scores of the 148 reference proteins (evaluator layout), the stage ledger, and from the learned
report: counts, unsupported reasons, and every EC candidate of a reference protein with its status, probability and
reason. Outputs stay on Shaheen scratch; small files are copied to ws.

Pass I: 22 × 9 arms exit 0; in every learned arm the learned stage executed with
`candidate_predictor_contract_unavailable` = 0; the output rows and claims are recorded (enforcement may remove rows
when the taxon applies; not a failure).

## Step E. Rows 1 and 4 (ws, after all cells)

Arms scored: `raw` (the DIAMOND claim scores themselves, read from the claims files), `stage0`, `stage1`, `learned`,
`learned_gbdt`. Metrics with the unchanged evaluator functions (`row1_calibration.py`: `pairs`, `choose_threshold`,
`evaluate_arm`, `calibration_metrics`, `contrast`):

E1 thresholds: the WP5 rule (≥ 10 calls, ≥ 90 % reference support, unconfirmed counts against) on the calibration
genomes of the frozen split, primary and swapped; `policy.json` written and hashed before any test metric.
E2 test metrics: top-1 (tie-aware, bootstrap), top-3, fixed 0.5 matched/unconfirmed, calibrated calls (calls,
supported, unconfirmed, support share = "supported precision", reference recall, bootstrap); also at the model's own
frozen display threshold (0.3212, chosen on its training calibration split) for the learned arms.
E3 calibration: positive-pair Brier, top-call Brier, 10-bin ECE, reliability tables; for the learned arms also
restricted to learned-scored candidates.
E4 contrasts at 0.5 and at the calibrated thresholds: `stage0 → stage1`, `stage1 → learned`, `stage1 → learned_gbdt`,
`raw → stage1`; corrections (reference pairs called only after), harms (only before), top-1 gained/lost, bootstrap.
E5 learned counts: candidates scored / unsupported (all types, and EC candidates of reference proteins), reasons.
E6 leakage sensitivity: E2–E4 recomputed without the L3 proteins.
Pass E: every element produced for every arm in both splits.

## Step F. Row 5 (ws)

F1: score every arm with the Price metrics (as E2, all 22 genomes). F2: differing scores between `gaps_learned`,
`universal`, `pathway`, `organism` (whole output and the 148 loci), learned counts per mode. F3: power gate restated:
gap-linked loci (unchanged definition: a reference EC is an EC of a gap record of its genome; pre-registered gate 20)
and, of those, the loci where these claims propose the reference EC. Pass F: F1–F3 produced. If fewer than 20
gap-linked loci exist the ablation is reported as unpowered whatever the score differences.

## Step G. Row 11

G1: the production extractor through `run_price_anomalies.run` (unchanged code; arms stage1, learned_logistic,
learned_gbdt; its own genus split, seed 20261101) on the 22 genomes' claims (Shaheen). Cases per arm and type.
G2: truth intersection with `prepare_price_agent_truth.py` (Price-148 complete ECs + Swiss-Prot experimental ECs,
exact sequence) and, because cases key on GO, additionally with GO-mapped truth: a case target GO term is
"positive" when ec2go maps a labelled EC of that candidate protein to that term; no GO negatives are derived
(unconfirmed is never negative). Truth is read only after the case bundles are frozen (hash recorded).
Truth whose Swiss-Prot accession is in the reference (L3) is reported separately.
G3: the six simple methods of `evaluate_stage3_agent.py` run only if at least one case with known truth exists in
both the calibration and the test split; the live agent only if its endpoint answers non-interactively with the
configured key (key never printed). Otherwise G3 is reported as not applicable with the measured reason.
Pass G: G1 and G2 produced; G3 run or its non-applicability reported with counts.

## Criteria for the acceptance rows (stated now, reported as measured)

- Row 1 "full system" element: shown only if the learned arm executes and scores Price-148 EC candidates.
- Row 4 bar elements unchanged; this plan adds Price-side calibration/harm evidence for the learned models; the
  null-outcome (complete locus labels) element is not addressed here.
- Row 5: Met requires a powered ablation (≥ 20 gap-linked loci); otherwise Not met with the measured differences.
- Row 11: correctness/harm counts as evidenced only if some method makes ≥ 1 known test decision under a frozen
  policy.

## Outputs

`docs/benchmarks/v3-price-contract-learned.md`; `benchmark/v3/results/price-contract-learned-v1/` with
`PLAN.sha256`, `SHA256SUMS`, `LARGE_FILES.tsv` (files > 5 MB and where they are kept); code in
`benchmark/v3/price_contract_learned/`; Shaheen wrapper `benchmark/v3/shaheen/price-contract-learned.sbatch`; one
progress line in `tasks/v3-supervision-plan-2026-10-28.md`. `tasks/v3-acceptance.md` and `README.md` are not edited.

## Deviations (append only, dated)
- 2026-11-02 (L, before any score was read). L1 finds no Price-148 protein with an exact or near-identical reference
  entry (best identity 94.1 %, one protein; three at ≥ 80 %), so the L3 sensitivity set is empty and the
  "excluded" population equals the full one. L2 finds 25 reference entries of a Price organism itself: Sinorhizobium
  meliloti 1021 (taxid 266834, genome GCA_000006965.1), named "Rhizobium meliloti" in Swiss-Prot, which the genus
  fence (first token of the organism name; "ensifer" fenced, "rhizobium" not) lets through; 15 of the Swiss-Prot
  experimental Price labels are these entries. Added sensitivity analysis, fixed now: every row-1/4 metric is also
  reported without genome GCA_000006965.1, and for row 11 truth from those 15 labels is reported separately as
  circular.
- 2026-11-02 (I, pilot Shaheen 15753032, GCF_000236665.1, no label read). All nine arms exit 0 and the learned stage
  scores 28,004 of 38,173 candidates in `learned` and `learned_gbdt` (0 contract-unavailable), but in `gaps_learned`
  and the three mode arms it abstains on every candidate with `unseen_source_function_role_or_availability_profile`
  (contract-unavailable 0): adding the gap observations and the reaction-local-context suggester gives every candidate
  a feature-availability profile absent from training. Five arms are added, fixed now and run on all 22 genomes:
  `obs_learned` (gap observations without the suggester, learned), `rlc_learned` (suggester without observations,
  learned) to locate the cause, and `cohort_universal`, `cohort_pathway`, `cohort_organism` (the three modes with
  cohort inputs and the learned model, without observations and suggester) so that the modes can be compared while
  the learned stage scores. The pre-registered mode arms remain the primary row-5 arms; the `cohort_*` arms are
  reported as a secondary analysis whose loci are not restricted to actual gaps. The pilot outputs are kept as
  `cells-pilot/`; the full run starts from a fresh `cells/`.
- 2026-11-02 (G1, execution). The first anomaly job (Shaheen 15753636) failed on every genome before any output:
  Java's temporary directory defaulted to the container's read-only `/tmp`. Rerun unchanged with
  `JAVA_TOOL_OPTIONS=-Djava.io.tmpdir=<job work dir>` (job 15753686); the failed attempt is kept on Shaheen as
  `anomalies-failed-readonly-tmp/`. The case bundle was frozen by `evaluate_stage3_agent.py prepare` inside Slurm
  (job 15754529, via the wrapper's new `MODE=py`), SHA-256 `88566c82…`, before any label was joined.
- 2026-11-02 (G2, recorded after the join). The unchanged EC wrapper `prepare_price_agent_truth.py` reports 0 pairs
  on this bundle because it reads `case["candidates"]`, whereas `prepare` nests them under `case["anomaly"]`; the
  plan's own join (`row11_truth.py`) reads the nested list (1,167 pairs). Both outputs are kept. G3 was not run: no
  known-truth pair exists in either split (the pre-registered condition), although the agent endpoint answered
  (HTTP 200 on `/models`).
- 2026-11-02 (I, recorded). In these single-genome batch commands (the row-5 command family) the strain taxon is not
  in the bundled hierarchy and symbolic enforcement is not enabled, so no taxon constraint is applied in any arm and
  output rows equal claims (1,234,318) in every arm; this matches the row-5 cells.
