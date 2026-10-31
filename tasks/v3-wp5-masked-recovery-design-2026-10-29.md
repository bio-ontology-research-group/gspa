# WP5 design: masked condition-dependent function recovery (items 8 and 12)

Status: **design only, nothing registered or run.** Requirement: handoff
`tasks/v3-resume-2026-10-28.md` §8 item 4. Acceptance rows 8 and 12 of
`tasks/v3-acceptance.md` stay open regardless of the outcome. The design builds on
the existing masked engine (`scripts/masked_assignment_experiments.py`). It does
not reuse the condition-planner benchmark (`benchmark/v3/run_planner_benchmark.py`),
which reads the target GPR when it builds its evidence.

## 1. Question and what can be claimed

Take a gene whose model assignment has been removed from every input. Given the
same target-blind candidate assignments (including null) and the same solver-call
budget, does a planner choose conditions whose simulated knockout phenotypes,
compared afterwards with the measured phenotypes, recover the hidden function at a
fixed precision? Function recovery, phenotype agreement and capacity are scored
separately. A result can show which arm recovers the most functions at that
precision under this protocol. It cannot show that a simulation confirms a
function. Agreement between two solvers is not biological correctness either.

## 2. Candidate sources found in the repository

| Source | In repo / local | Organism, model | Panel is gene-independent? | Role in WP5 |
|---|---|---|---|---|
| Maranas 2026 Keio auxotrophs (`prepare_metabolic_benchmarks.py`; derived files under `~/tmp/evacuated-from-tmp-2026-10-26/gspa-v3-metabolic-data/derived/`: `auxotroph_validation_labels.csv`, `cases_manuscript_regulatory.jsonl`, `masking.json`, `regulatory_closures.json`) | yes | *E. coli* Keio (BW25113 parent), iML1515 | **No**: the measured substrate per locus was chosen by an iML1515 design that used the target GPR | Tier A, calibration split; menu = full 24-substrate grid + glucose (as in the planner benchmark) |
| Baichman-Kass 2026 nutrient auxotrophs (`benchmark/v3/results/nutrient-environments-v1/data/derived-v3/{inputs,observations}.jsonl.gz`) | yes | *E. coli* EcNR1 / iML1515; *B. subtilis* 168 / iYO844 (`register_nutrient_models.py`) | yes (44 environments ± supplement, same for all six deletions) | Tier A, test split |
| Eichhorn 2000 sulfonate transport (`benchmark/v3/results/sulfonate-diagnostic-v1/dataset/`) | yes | *E. coli* EC1250 (MC4100), iML1515 | yes (20 sulfur sources) | Tier A, test split, **transport-role negative control** (single-gene targets tauA, ssuA). Multi-gene genotypes serve only as background-knockout conditions |
| Price 2018 RB-TnSeq (`benchmark/manifests/suite6_price149_manifest.tsv`, six predictor matrices × 22 assemblies) | EC truth and predictions only | 22 non-model bacteria, **no organism models** | yes in the source (all genes in all experiments), **but no condition-level fitness data in the repo** | Tier B, gated (§9) |
| gapseq/gapsmith media | **none in repo** (only the adapter and one observed gapsmith run) | — | — | Needed only for Tier B |

> [!WARNING]
> Price manifest columns `cognate_exchange_rxn`/`cognate_substrate_name` are a
> constant glucose placeholder (`build_all_manifests.py`), not measured conditions.
> `metabolic_stratum`, `seed_reaction_ids` and `in_seed_universal_model` are derived
> from `ground_truth_ecs`. Any of them would leak truth and none may be used.

## 3. Case selection: the "unambiguous" checklist

A (case, condition) pair enters the primary analysis only if every item holds.
Failures go to `ineligible.jsonl` with the reason; they are never dropped
silently. **No item may read a measured outcome or the model's agreement with it.**

1. **Strain.** The measured strain is the model organism or a documented derivative.
   Background deletions of the parent strain are applied as fixed background
   knockouts, or the affected conditions are excluded. Examples: BW25113
   Δ(araD-araB) Δ(rhaD-rhaB) ΔlacZ; MC4100 Δ(argF-lac). The derivation is cited.
2. **Gene identity.** The deleted locus maps to exactly one model gene by locus tag
   (symbol mapping only through the source's own index). The genotype is a
   single-gene deletion. Multi-gene genotypes are admissible only as
   background-knockout conditions for another target. For Tier B transposon data
   there must be no same-strand downstream gene within the operon rule (≤300 bp)
   with the same phenotype; this polar-effect rule is frozen before the join.
3. **Function truth.** T(g) is the set of reactions whose original GPR contains g
   (`masking.json` style, private). Primary: 1–2 internal (non-transport)
   reactions. Stratum *complex subunit*: g sits in an AND clause with ≥1 partner.
   Stratum *transport-only*: every reaction in T(g) crosses a compartment. The
   correct answer for this stratum is null or abstain, and any catalytic call is
   harmful. More than 2 reactions, or mixed transport and internal reactions:
   `ambiguous_truth`, reported only.
4. **Medium.** Every varied component maps to an exact model exchange. Unrepresented
   ingredients are listed. The protocol is the registered one (Keio
   `manuscript_regulatory`; nutrient `printed` and `additional_inorganic` arms,
   both reported; sulfonate frozen M63 binding). No medium is tuned on labels.
5. **Measured phenotype.** The wild type grows in that condition: a measured control,
   or the parent phenotype stated by the source; otherwise `wt_unverified`. The
   label is binary under the source's own thresholds and stable across its
   threshold grid. Partial growth, borderline rows, rows whose raw OD does not
   reproduce, and unknown values stay in a separate sensitivity class.
6. **Assay ↔ model match.** The endpoint is growth or yield, not a product, and the
   condition can be expressed as structured `mods` (carbon/nitrogen/sulfur swap,
   supplement, oxygen, measured background gene knockout).

Splits (frozen; gene-disjoint): calibration = Keio (*E. coli* catabolic loci); test =
nutrient (both species, includes an organism holdout) + sulfonate controls (+ Tier B
when it exists). Preliminary Keio scoping from the private truth file only:

| Keio truth class (44 loci) | Loci | Use |
|---|---:|---|
| internal, ≤2 reactions (includes sdhCDAB-type subunits) | 25 | primary / complex stratum |
| transport-only | 16 | transport negative control |
| mixed or internal >2 reactions | 3 | `ambiguous_truth`, reported |

## 4. Masking: every channel, and mechanical verification

| Channel | Where the target can leak | Rule |
|---|---|---|
| Model | SBML `fbc:geneProduct` id/`fbc:label`/`fbc:name` (e.g. `tpiA`), ncbigene/UniProt annotations, GPR clauses | `mask_model` removes whole AND clauses and every alias of the target. External tools (COBRA checker) only ever receive a re-serialized masked SBML with the geneProduct deleted. Orphaned reactions stay in the network. |
| Candidate construction | sequence DB self-hits; neural predictors trained on the target (CLEAN etc. on *E. coli*/*B. subtilis*; possible Price-derived curated entries in Swiss-Prot); `reactions.tsv`/refiner models; learned scores trained on labelled cohorts | Homology only (DIAMOND vs Swiss-Prot) with the target accession and same-genus entries removed, using a frozen DB hash. Neural predictors and learned scores are excluded from the primary arms; neural predictors run as a labelled leaky sensitivity. EC→reaction mapping comes from reaction-level SBML EC annotations only. Gene rules and gene annotations are never read. |
| Sequence/context inputs | NCBI FASTA headers and GFF `product`/`Note`/`EC_number`; claim `other_annotations`/GO; cohort orthologs carrying curated annotations | Headers reduced to IDs (as `prepare_price_layout.py`); coordinates/strand only; same-species genomes excluded from recurrence. |
| Case files | `cases*.jsonl` `case` block, `masking.json`, labels CSV, Price `ground_truth_ecs`/`seed_reaction_ids`/`metabolic_stratum` | Stored under `private/`; only `evaluate_wp5_recovery.py` opens them. |
| Conditions | reaction-level knockouts of T(g); menus restricted to measured design substrates (Keio) | Condition mods may not name reactions; Keio menu = full grid. |
| Product objectives | `product_exchange`, `product_demand`, `baseline.product_demands`, `bgc` evidence; objectives that require the candidate | **Banned in WP5**: no measured product labels with negatives exist. The runner refuses any public problem containing them. |
| Agent text | gene symbol/locus in the prompt (`run_planner_benchmark.make_prompt` uses `locus`); IDs embedding the locus (`H:{locus}:{substrate}`); organism name; partner/background gene IDs Qwen can decode from memory | Opaque `H<n>`/`C<n>` (already enforced). The chat-transport wrapper replaces every model gene ID with a private alias `P<n>` and omits organism/strain. Only opaque condition IDs come back. |
| Evidence free text | `--experiment-evidence` `source` strings | Private construction report only; public evidence is numeric channels. |

Mechanical verification (`audit_wp5_masking.py`; must pass before any arm runs;
output `masking_audit.json` bound to input hashes):

- **Alias scan.** Every target alias (locus tag, old locus tags, symbol, b-number,
  UniProt/RefSeq/GenBank protein IDs, GFF/UniProt product name) is searched with
  token boundaries over every byte of: public problem, claims and evidence, FASTA
  headers, menus and constraints, the protocol, and **every Qwen request in the
  transcripts**. Zero hits are required.
- **Hidden-truth invariance.** Rebuild the public bundle with (a) the target's GPR
  moved to a decoy reaction set, (b) truth ECs permuted, and (c) all measured
  labels flipped. Public-problem SHA-256, deterministic-arm selections and the
  first Qwen prompt must be byte-identical.
- **Orphan invariance.** Candidate construction from the masked and the unmasked
  model gives identical hypotheses. A builder that reads orphan status fails.
- **DB exclusion.** No hit to the target accession or to the excluded taxa in the
  candidate DB; the DB hash is recorded.
- **Objective/condition check.** Zero product conditions and zero reaction-level
  condition mods.

## 5. Arms and the identical budget

All arms receive the identical public problem (hash checked per run). The cost of
one condition is c = Σ_h (2 + 2·|a_h|) LPs: a control and a target knockout per
hypothesis, plus an FVA min and max per assigned reaction (alternative flux). The
existing engine reserves the full rectangular cost before it selects a condition.
Budget per case: `max_solver_calls = k·c`, `max_experiments = k`, with k ∈ {3, 6}
frozen. Unused calls are reported and never reallocated.

| Arm | Engine strategy | Definition |
|---|---|---|
| `minimal_medium` | `fixed`, panel = [baseline] | Only the source's minimal medium (M9 glucose, MSgg glucose/NH4, sulfate). Same cap. |
| `fixed_panel` | `fixed`, panel = frozen source-table order | Gene-independent; registered in `protocol.json`. |
| `random` | `random`, seeds 0–9 | Seeded shuffle of the full menu. |
| `evidence` | `evidence` | Existing rule: sum of per-channel maximum hypothesis evidence over `evidence_channels`; ties broken by ID. |
| `qwen` | `qwen` | `qwen3.8-27b`, unimatrix01 endpoint via `GSPA_LLM_URL`; key loaded by `ME.load_api_key()` and never printed. Temperature 0, seed 0, thinking off, timeout 600 s. An invalid reply or endpoint failure stops the arm, with no fallback credit. |

Reference rows (not arms): `sequence_only` (0 LPs: top evidence hypothesis) and
`exhaustive` (every menu condition, outside any budget). Each arm cell must equal
the exhaustive cell for the same (case, hypothesis, condition).

## 6. Outcome categories

- **Cell** (hypothesis × condition): `both_grow`, `knockout_no_growth`,
  `knockout_infeasible_with_maintenance`; `no_growing_control`; `inconclusive`
  (solver error, budget). FVA per assigned reaction: `required`, `optional`,
  `no_flux_at_floor`, `unassessed`.
- **Measured** (case, condition): `growth`, `no_growth`, `partial`, `unknown`,
  `unmeasured` (common for the Keio grid), `wt_unverified`.
- **Comparison**: `agree`, `contradict`, `not_evaluable`.
- **Case decision**, under a frozen rule identical for all arms. C is the set of
  non-null hypotheses with zero contradictions. A *supporting* condition is an
  agreeing condition where hypothesis predictions differ.
  - `correct_call`: unique top-evidence member of C with ≥1 supporting condition,
    whose reactions match T(g).
  - `wrong_call` (harmful): the same, but the reactions do not match; this includes
    any catalytic call on a transport-only target.
  - `null_call`: C is empty and null is supported. It is `correct_negative` if T(g)
    has no candidate or the case is transport-only, otherwise `false_null`.
  - `abstain_tied`, `abstain_unsupported`, `no_consistent_hypothesis` (every
    hypothesis contradicted).
  - `inconclusive` (solver/agent failure), `no_admissible_experiment`.

Partial T(g) overlap is reported separately, never as correct.

## 7. Metrics and denominators

| Metric | Unit | Denominator | Notes |
|---|---|---|---|
| **Function recovery at fixed precision** (primary) | case | **all eligible test cases**, including candidate-generation failures | Per-arm threshold on decision score (evidence score of the call) frozen on the calibration split with the `evaluate_stage3_agent.calibrate` rule: precision ≥ 0.90, ≥ 10 decisions, unknowns count against precision. If no threshold qualifies, the arm abstains and scores zero (reported). |
| Harm, abstention, false null | case | same | Reported beside recovery. |
| Candidate recall | case | all eligible | Truth in the candidate set. Identical across arms; separates generation failure from planning. |
| Conditional recovery | case | cases with truth in candidates | Descriptive. |
| Test precision–recovery curve | case | all eligible test | Descriptive only; not the fixed-precision result. |
| **Phenotype agreement** | (case, visited condition) | binary-evaluable visited measured pairs per arm; also the common exhaustive set | Prediction of the called hypothesis vs measurement. Reference rows: unmasked original model and constant growth/no-growth comparators on the same pairs. Partial growth reported separately. |
| **Capacity** (assumption: model-conditional, not recovery) | case / cell | eligible cases with ≥2 hypotheses; growing-control cells | C1: visited set contains a discriminating condition. C2: FVA class distribution of called and true assignments. |
| Cost | run | all runs | LPs used vs cap, wall time, Qwen requests, failures, tokens, seconds. |

Statistics: exact McNemar on paired case outcomes (Qwen vs each comparator; every
arm vs `minimal_medium` and `sequence_only`). Wilson intervals for proportions.
Bootstrap over genes (10,000 resamples, frozen seed). Random arm: per-case mean
over seeds. Per-stratum tables cover catalytic, complex, transport-control and
species.

## 8. Scripts, tests and order

1. `benchmark/v3/prepare_wp5_cases.py`: checklist, splits, menus, `public/` vs
   `private/`, `ineligible.jsonl`, hashes. Test `tests/test_wp5_cases.py`: each
   exclusion reason fires; label-flip invariance; Keio menu independent of the
   measured pair; Price placeholder columns rejected.
2. `benchmark/v3/build_wp5_candidates.py`: filtered-DB DIAMOND claims (cluster) →
   production `MaskedAssignmentBuilder` via `gspa integrate --experiment-mode
   evidence_guided` (panel conditions supplied through `assignment_recovery.conditions`,
   cap `--experiment-max-hypotheses 8`) → frozen public problem. Test
   `tests/test_wp5_candidates.py`: header stripping, DB exclusion, no GPR reads.
3. `benchmark/v3/audit_wp5_masking.py` (§4). Test `tests/test_wp5_masking_audit.py`:
   an alias injected into each channel (FASTA header, claim metadata, prompt,
   condition mod, product condition, SBML `fbc:name`) fails the audit; a GPR-reading
   stub builder fails invariance.
4. `benchmark/v3/run_wp5_recovery.py`: stages `protocol`, `exhaustive`,
   `deterministic`, `qwen --shard i/n`. Has no private-path argument and wraps the
   Qwen transport for aliasing and usage. Test `tests/test_wp5_recovery_runner.py`:
   identical public hash across arms, equal caps, `minimal_medium` visits only the
   baseline, agent failure gives no credit, no private file opened (guarded `open`).
5. `benchmark/v3/evaluate_wp5_recovery.py`: phases `join`, `calibrate`, `evaluate`,
   `report`; refuses a policy whose hash does not match. Test
   `tests/test_wp5_evaluate.py`: decision truth table for every §6 class, every
   denominator, transport-control harm, partial-growth exclusion.
6. `benchmark/v3/check_wp5_cobra.py`: COBRApy/GLPK replay through
   `check_masked_assignment_cobra.py` with one model copy per cell.
7. `benchmark/v3/wp5_recovery.sbatch`, then the report
   `docs/benchmarks/v3-wp5-masked-recovery.md`.

Run order: 1 → 2 → 3, with the audit gate before any LP. Then 4 on a synthetic
fixture, 5 on synthetic labels, and focused tests locally (stdlib/mock only;
pattern `test_wp5_*.py`). Then on the cluster: DIAMOND, audit, exhaustive +
deterministic arms, Qwen shards, COBRA replay. Calibrate, freeze the policy and
only then join test labels.

## 9. Compute budget (cluster only; no local LP, JVM or DIAMOND)

| Job | Estimate | Resources |
|---|---|---|
| Filtered Swiss-Prot DB + DIAMOND on 2 proteomes | <1 CPU-h | node002/003, 8 CPU, 16 GiB |
| Builder through production CLI (JAR 10526, hash-bound) | minutes per case | 2 CPU, 8 GiB |
| Exhaustive sweep (~60 cases × ≤90 conditions × c ≈ 45 LPs) | ~2.4·10⁵ LPs, ≤12 CPU-h | 8 CPU, 16 GiB, 24 h |
| Arms (14 runs × 2 budgets × ≤6 conditions × c) | ≤4.5·10⁵ LPs, ≤24 CPU-h | same |
| Qwen (≤9 requests per case across budgets; ~90–115 s each, as observed earlier) | ~4 h wall at 4 concurrent | no GPU job; do not touch node005/node008 services |
| COBRA replay (all arm-visited cells; 10% seeded sample of the rest if >24 CPU-h) | ≤24 CPU-h | 4 CPU, 16 GiB |

Frozen Price/CLEAN jobs (10297/10298 → 10410 → 10432; 10367 → 10368 → 10406) are
not touched.

**Tier B (Price) gate.** It may start only after all four hold:
- condition-level fitness and experiment-media tables are acquired from the Fitness
  Browser / paper deposit, with URLs, versions and hashes;
- the phenotype thresholds are quoted verbatim from Price 2018 and frozen;
- target-blind organism models exist: a curated exact-strain model where one exists
  (e.g. *S. meliloti* 1021, *S. oneidensis* MR-1), otherwise a gapseq draft built
  from the proteome with the target removed and gap-filled only on a minimal medium
  outside the test panel;
- the engine supports adding candidate reactions that are absent from the model. It
  currently requires assignments to name existing reactions.

The 10432 outputs only describe raw predictor candidates for Price targets. They
are not phenotype labels.

## 10. What counts as a negative result

All of the following are reported plainly. None triggers retuning.

- Qwen's fixed-precision recovery does not exceed `evidence` (difference ≤ 0, or
  McNemar p ≥ 0.05): **no agent benefit**.
- No arm exceeds `minimal_medium` or `sequence_only`: **experiments add no recovery**.
- No calibration threshold reaches 0.90 with ≥10 decisions: that arm **cannot
  recover at fixed precision**.
- Low candidate recall: the **bottleneck is candidate generation**, not planning.
- Harmful calls ≥ correct calls in any arm, or any catalytic call on a transport
  control.
- Called-hypothesis phenotype agreement is at or below the constant comparator.

A failed masking audit is an **invalid run**, not a negative result. It is fixed
and the run repeated, with the failure retained. Small denominators (tens of test
cases) are stated together with every interval.

## 11. Open decisions for the supervisor

1. The reading of "capacity" (C1/C2 above) versus product/secretion capacity, which
   is excluded because no labels exist.
2. Whether Keio may serve as calibration, given its design-selected panel.
3. The hypothesis cap (8) and budgets k ∈ {3, 6}, frozen before any run.

## 12. Supervisor decisions (2026-10-29)

1. Capacity is read as model-conditional discrimination (C1) and FVA class distribution (C2); product or
   secretion capacity is out of scope until measured product labels with negatives exist.
2. Keio serves as the calibration split with the full 24-substrate grid plus glucose as the menu; the
   design-selected origin of its measured pairs is stated in every table that uses it.
3. Hypothesis cap 8 and budgets k ∈ {3, 6} are frozen as written. Price (Tier B) stays gated on the four
   conditions of §9.

## 13. Supervisor rulings on the implementation's open decisions (2026-10-29, second round)

The first implementation pass (Opus sub-agent, report `~/tmp/gspa-v3-agy/wp5-impl/report.md`) stopped before the
evaluator, tests, COBRA checker and job scripts existed. Its open decisions are ruled as follows; the completion
(pi Task N) implements them.

1. Nutrient cases: oxygen uptake fixed at 20 mmol gDW⁻¹ h⁻¹ (aerobic plate assay) — accepted, stated in every
   table; no oxygen sweep in WP5. The gene-independent menu of 44 environments × (none + six supplements) is
   accepted; its cost is reported, not trimmed per case. `printed` is the primary medium arm and
   `additional_inorganic` the sensitivity arm. The supplemented wild type is the unsupplemented wild type.
2. Keio: wild-type growth as stated by the source; the glucose baseline is labelled unmeasured and never a
   label; substrates affected by the BW25113 background deletions are excluded — accepted.
3. Sulfonate: sulfur ceiling 1.0; the background-knockout sets exclude both targets, so tauA and ssuA share a
   menu — accepted. EcNR1 background deletions cited but not applied — accepted (no condition varies biotin).
4. Candidate construction with a JVM budget of one experiment and one solver call (Java rejects zero) is
   accepted only if the construction receipt records that budget and the result validation accepts a cell
   list of at most one; a test must cover that. If the validation rejects it, stop and report.
5. The runner must not reimplement the engine's selection loop. Preferred: call the production selection
   (`run_assignment_experiments` or the lowest-level engine entry that takes an already-masked model and the
   public hypotheses) so that no target gene identity is needed. If the engine API cannot run on a masked
   model without `mask_gene`, keep the reimplementation but add a parity test on the synthetic fixture that
   proves identical selections and identical cells to the engine for the `fixed`, `random` and `evidence`
   strategies, and record the deviation in the report and in the results README.
6. The runner's own Qwen transport (600 s timeout, temperature 0, seed 0, thinking off) is accepted; the
   alias wrapper (`P<n>` gene aliases, opaque `H<n>`/`C<n>`, no organism or strain) is mandatory and tested.
7. Truth and labels split into separate private calibration and test files — accepted; `evaluate` joins test
   labels only after `calibrate` froze the policy hash.

## 14. Implementation outcome and accepted deviations (2026-10-29, Task N)

Scripts, 100 tests (`test_wp5_*.py`, synthetic fixtures and mocks), `check_wp5_cobra.py` and both job scripts exist
(report `~/tmp/gspa-v3-agy/task-N/report.md`). §13.5 was resolved the preferred way: every arm, including the
exhaustive reference, runs through the production `run_assignment_experiments` loop on the masked public SBML
with a placeholder `mask_gene` that occurs in no GPR; the engine's `public_problem_sha256` must equal the frozen
hash. The supervisor accepts the following deviations, each disclosed in the results:

1. The §7 "unmasked original model" reference row is reported as `not_available` (no public-only stage can
   compute it); it may be added later by a separate private LP stage.
2. The Qwen arm and the transcript audit run on unimatrix only (the endpoint is unreachable from Shaheen);
   Shaheen runs the deterministic arms and their COBRA replay (`--arm-stages deterministic`).
3. Hidden-truth invariance (b) and (c) are realised as scratch-bundle rebuilds whose private files really are
   permuted or flipped; the first-prompt comparison uses an LP-free reconstruction proven equal to the engine's
   request by test.
4. In evidence-guided mode the production builder may add target-blind conditions derived from the hypotheses
   (e.g. `oxygen_removal`); they are shared by every arm, stay unmeasured, and `fixed_panel` uses only the frozen
   menu.
5. The one-experiment / one-call construction budget yields an empty cell list because the production validation
   demands a complete visited × hypotheses rectangle; this was established by code reading and must be confirmed
   by the first cluster `construct` stage before the audit gate is trusted.
6. The random arm reports bootstrap intervals only (fractional per-case means over seeds); the "unknown counts
   against precision" term is structurally zero because a missing truth row fails closed.
7. `abstention` includes `no_consistent_hypothesis`; every class is also counted separately.

Also recorded: the design's Keio input path (`~/tmp/evacuated-from-tmp-2026-10-26/...`) no longer exists on ws;
the file must be fetched from the laptop or unimatrix, or regenerated by `prepare_metabolic_benchmarks.py`.

### §14 addendum (2026-10-29 18:20, supervisor ruling pending Robert's confirmation)
- The `masked_model` alias-scan channel is not among the channels listed in §"Mechanical verification"; the
  implementation added it as a whole-file, case-insensitive token scan. On the real bundle it produced only
  collisions with non-gene content (BioCyc `META:` prefix vs `metA`; reaction-name word "out" vs a dctA synonym;
  metabolite id `mpt`; other genes' refseq synonyms; shared InterPro ids). Ruling: the gating scan of a masked
  model covers the masking boundary of the model row of the table in §4 — every `fbc:geneProduct` element
  (attributes, notes, annotations) and every `fbc:geneProductRef` — with case-sensitive symbols; the whole-file
  scan is kept as an informational, non-gating channel with recorded contexts. A hit inside a geneProduct of
  another gene is still reported (not suppressed).
- Public masked SBML must keep the core namespace as the default namespace (JVM reader compatibility).
- (18:55) Structural scan of the 58 regenerated masked models: 162 target-alias hits in 24 cases, none in the
  target's record or in any geneProductRef; all in other genes' geneProduct annotations: operon/complex names on
  the other subunits (manXYZ ptsX/gptB/mpt, srl/gut sbl/gutA/srlA, malJ, ctr, yidU, arg*, malA, lct) and shared
  InterPro family ids on other B. subtilis proteins. Ruling: an alias that the ORIGINAL unmasked model already
  carries on another gene's geneProduct is a shared identifier; its hits on that other gene are adjudicated
  (channel `masked_model_shared`, non-gating, gene/context kept in private detail.json, counts only in the public
  report). geneProductRef hits, hits in the target record and aliases not shared in the original always gate.
  Result: 0 gating hits, 162 adjudicated. Operon/complex names remaining on co-subunits is intended (the design
  keeps the rest of the network and its annotations; only the target's own record is removed).
