# Frozen plan: reaction context on masked enzyme loci with known truth (acceptance row 5)

Date: 2026-11-02. Status: FROZEN on 2026-11-02 after Robert's approval (section 10); the SHA-256 of this file is
recorded in `benchmark/v3/results/row5-masked-context-v1/PLAN_SHA256` before any arm output is read. Later changes
are appended under "Deviations", dated; the body above that heading is not edited.

Requested by Robert (2026-11-02): "think about Point 5 and how to improve it; 20 may simply be too high a threshold;
ANY improvement would be good, if it is systematic and consistent."

## 1. What has already been seen (stated so that nothing below is post hoc)

- `docs/benchmarks/v3-price-calibration-context.md`, row 5 (462 arms, 22 Price genomes, 6 enzyme predictors, real
  gapsmith gaps): 8 of 148 Price loci are gap-linked (gate 20); the universal, pathway-only and organism-specific modes
  give 0 differing scores because they feed only the learned stage, which abstained on 100 % of candidates (no
  predictor contract on these predictors' claims); the gap prior changes 29 to 75,852 scores per predictor but no
  Price call and no top-1 answer. The mode arms were not run for CLEAN, DeepProZyme and GraphEC (9 to 42 h per cell).
- Earlier: three-mode comparison on the 18-genome cohort (`v3-context-comparison-results.md`, no consistent context
  advantage, 0 gap-linked test positives), observed-gap cohort (`v3-observed-gap-cohort.md`, 521 gap records, 202
  score changes, no labels on them), Stage 2B and follow-ups (`v3-price-stage2b.md`, `v3-stage2b-followups.md`: no
  Price gain; the top-10 candidate limit removes Stage 2B's only harm).
- Those 8 gap-linked Price loci are not re-analysed under a lower gate. They re-enter below only as members of a new,
  larger locus set under a new criterion, and are also reported separately.
- Design-time looks taken for this draft (2026-11-02, `~/tmp/row5-design/`): the experimental labels of 40 genomes
  were joined to their gapsmith pathway-search tables (label-dependent, score-blind), and, to size the room for
  improvement, the already published stage-0 scores of the labelled proteins (rc1 retained files) and the raw DIAMOND
  claim scores were read. No masked arm exists, so no arm outcome was seen.

## 2. Why the previous test could not show an effect, and what a fair test needs

Three separate reasons, each sufficient on its own:

1. Most Price enzymes are not gaps of their own genome's model. Of 152 Price (protein, EC) pairs, 63 are
   "redundant" (another protein of the genome also has a strong gapsmith hit for the reaction), 26 have an EC that
   occurs in no reaction row of the table, 20 sit in pathways gapsmith does not predict, 18 fall under the gap
   builder's exclusions (spontaneous, complex, reference unknown); 25 would become gaps if the protein's own hit were
   removed. Only 8 are gaps as observed.
2. The three context modes act only through the learned stage (`IntegrationWorkflow.runStages`, `consumed_by`), and
   the learned stage scores only claims that carry a predictor contract. The six Price predictors carry none.
3. The gap prior (`GapFillingPrior`) is mode-independent. It rescues one candidate per gap, chosen by score tier, only
   when the candidate's likelihood lies in [0.18, 0.5]. It never looks at where the pathway's other enzymes sit on the
   chromosome.

A fair test therefore needs: (a) many loci whose truth is known and whose reaction is a gap *with pathway context
present*; (b) claims the learned stage accepts, so the modes can act; (c) outcome measures that can register what
context can actually change (call status and which gene gets the gap, not only top-1 of a protein that is already
right); (d) a decision rule that rewards a small but consistent effect and does not need 20 loci to say anything.

## 3. Design: make gaps where the truth is known

### 3.1 Genomes and labels

The 145-genome gapsmith benchmark set (`/scratch/project/k10018/gspa-v3/benchmark-models-v1`, 22 Price + 18 cohort
+ 105 recurrence peers, gapsmith 0.1.0 7b7ff538, the models of the off-growth experiment) is the model source. Only the
40 genomes with experimental labels are used; the 105 peers are recent unstudied isolates without Swiss-Prot
experimental entries.

| Label set | Source | Pairs (protein, complete EC) |
|---|---|---|
| Price-148 | `benchmark/manifests/suite6_price149_manifest.tsv`, 22 genomes | 152 |
| Price Swiss-Prot | `benchmark/v3/offgrowth/data/swissprot_experimental_ecs_price22.tsv` (ECO:0000269 and descendants, exact sequence) | 35 |
| Cohort Swiss-Prot | `~/.cache/gspa-v3-data/training-cohort/labels/labels.tsv`, tier experimental, 18 genomes | 335 |

Unconfirmed is never negative: a call naming another EC on a labelled protein is "unconfirmed"; proteins without
labels are never counted either way.

### 3.2 Which loci can become gaps (design-time count, `~/tmp/row5-design/count_maskable.py`)

For each labelled (genome, protein P, EC E), the (pathway, reaction) groups of the gapsmith reaction table whose EC
field contains E were classified under `PathwayGapBuilder`'s own rules, before and after setting every `good_blast`
row that names P to `no_blast` (the masking operation of 3.3):

| Class | Cohort Swiss-Prot | Price-148 | Price Swiss-Prot | Total |
|---|---:|---:|---:|---:|
| maskable with context (becomes or is a gap; another step of the pathway has a strong hit from another protein) | 87 | 25 | 7 | **119** |
| maskable, no other step present | 15 | 0 | 1 | 16 |
| redundant (another protein also has a strong hit for the reaction) | 65 | 63 | 6 | 134 |
| EC in no reaction row of the genome's table | 122 | 26 | 10 | 158 |
| pathway not predicted by gapsmith | 27 | 20 | 7 | 54 |
| excluded by the gap builder (spontaneous, complex, reference unknown) | 19 | 18 | 4 | 41 |

The 119 maskable loci lie in 28 genomes (12 Price, 16 cohort; the largest: Thermococcus kodakarensis 17, Haemophilus
influenzae 12, Sulfurisphaera tokodaii 11). 107 become gaps only through masking; 12 are gaps already (these include
the previously seen Price loci). Split membership: cohort train 58, cohort calibration 23, cohort test 6, Price
calibration 17, Price test 15. Pathway prediction after masking: for 71 loci the pathway stays above gapsmith's
completeness threshold without the masked reaction and the reaction is not a key reaction ("robust"); 48 are
"marginal" (reported as a stratum; the observed prediction flag is kept, see 3.3). 28 (genome, pathway) pairs hold two
or more maskable loci.

**Held-out loci for the primary test: 61** (cohort calibration + test 29, Price 32), that is, loci on genomes the
learned models were not trained on. The 58 cohort-train loci are a development set (used for 5.3 only, reported
separately, never pooled into the primary test).

### 3.3 Masking operation (label-dependent, score-blind; no engine change)

One cell per masked locus (P, E) in genome G:

- Reaction table: every row whose subject (`stitle` first token) is P and whose status is `good_blast` becomes
  `no_blast`, in every pathway (not only the pathway of E), so that no other reaction of P stays as context for P.
  Rows of other proteins are untouched. The pathway table is unchanged (the prediction flag is kept as observed; the
  "robust / marginal" stratum of 3.2 says how much the flag owed to the masked reaction).
- The observation JSON (`gspa.pathway-search-observations/1`) is built from the masked tables with the row-5 converter
  (`benchmark/v3/price_calibration_context/build_context_inputs.py`, same conflicting-identifier rule); gap derivation
  runs inside the engine as before (`--pathway-observations`).
- The organism reaction set (ORGANISM mode): the draft model's reactions minus the reactions of E that had no other
  gene (so the organism set does not carry the masked enzyme), computed from the draft SBML gene-reaction rules.
- Claims are **not** masked: the question is whether context raises the right EC on the right gene among the claims
  the predictors already make. (The off-growth experiment masked claims as well; that answered a different question.)
- Masking audit per cell: no `good_blast` row names P; the gap record for E's reaction exists in `<out>.gaps.jsonl`
  (`available`); the number of gap records added relative to the unmasked genome is recorded.
- Loci sharing a pathway (28 pathway pairs) are masked one at a time, so each locus keeps the other as context.

### 3.4 Claim sets

| Claim set | Genomes | Learned stage | Source |
|---|---|---|---|
| DIAMOND Swiss-Prot, genus-fenced, contract-bearing | 22 Price | executes | the other session's Step C output (`/mnt/data2/leechuck/gspa-v3-data/price-contract-learned-v1/claims-run/claims/`), read only, used by SHA-256 |
| same | 18 cohort | executes | produced on ws with the same driver (`benchmark/v3/routes_closure/predict_contract_claims.py`, same reference `1d8c4ebb…`, DIAMOND 2.1.11, same parameters) on the cohort proteomes; the existing baseline claims carry no contract |
| six Price predictors (CLEAN, DeepProZyme, EnzBERT, GraphEC, MAPred, TopEC) | 22 Price | abstains (known) | `price-full-matrix-v2`, as in rc1 |

Dense-predictor claim filter (CLEAN, DeepProZyme, GraphEC), score-blind with respect to labels: keep, per protein,
its 10 top-scoring ECs and, per EC, its 64 top-scoring proteins plus ties (the off-growth deviation 3 rule, shown
lossless for top-1 attribution). Measured rates (3.4 s per 1,000 claims in the cohort pass; 50 × stage-0 memory)
make the full matrices infeasible (9 to 42 h per mode arm); the filter brings a CLEAN cell from 45 M to under 0.4 M
claims. Check: the filtered `stage0` reproduces rc1's top-1 and top-3 on every labelled protein of the genome; a
difference is reported, not repaired. The three small predictors run unfiltered.

### 3.5 Arms per cell

| Arm | Command relative to the rc1 production command | Claim sets |
|---|---|---|
| `stage0` | `--enable-priors none` (must equal rc1 / the unmasked stage 0) | all |
| `gaps` | default priors + masked `--pathway-observations` + `--rlc-suggester` + reaction graph and EC aliases ("gap prior only") | all |
| `gaps_learned` | `gaps` + `--learned-model` (contract-3 logistic), no cohort inputs | DIAMOND |
| `universal`, `pathway`, `organism` | `gaps_learned` + cohort inputs (single-genome, clade holdout) + `--cohort-context-mode` | DIAMOND |
| `unmasked_gaps` | `gaps` with the unmasked observation (one per genome, not per locus) | all |
| `ctx_universal`, `ctx_pathway`, `ctx_organism` | `gaps` + the proposed context-weighted gap prior of §5 (only if §5 is approved and built; a new method version) | all |

Engine: the pinned release-candidate JAR `08f8b5a1…` for every arm except the `ctx_*` arms, which need the new build
(its SHA-256 recorded; §5). Cohort inputs for the 18 cohort genomes: the verified context-mode inputs
(`cohort-context-modes-v2-inputs`), as in the three-mode comparison.

### 3.6 Outcomes per masked locus (P, E), per arm, read from the arm's output

Three indicators, each gained / lost / unchanged relative to the comparison arm:

1. **Call**: integrated score of (P, E) ≥ 0.5 (the fixed diagnostic threshold; no predictor obtained a calibrated
   threshold in row 1). This is where the gap prior and the learned stage can act on a true EC that is a candidate
   below threshold (design-time: 26 cohort, 10 Price DIAMOND loci have the true EC in the prior's window).
2. **Protein top-1**: E is P's highest-scoring EC.
3. **Gene for the gap**: P is the highest-scoring protein of the genome for E (reaction-centric; the only indicator
   with room on the dense predictors, where the true EC often sits at the predictor floor among thousands).

Locus outcome: *improved* if at least one indicator is gained and none lost; *harmed* if at least one is lost and
none gained; *mixed* counts as harmed; otherwise *unchanged*. Harm elsewhere, per cell: labelled (P', E') pairs of
the same genome (all three label sets) whose call or top-1 is lost; unconfirmed calls created on labelled proteins are
reported as unconfirmed, and the genome-wide number of score changes and new calls ≥ 0.5 is reported as cost.
Also recorded: the RLC suggestion for E's reaction (does it name P; rank of P), the learned stage's scored /
abstained counts and reasons, and the per-mode local-context signature size for E's reaction.

### 3.7 Decision rule replacing the power gate (pre-registered; no minimum locus count)

Contrasts: C1 `stage0 → gaps` (all claim sets); C2 `gaps → gaps_learned`; C3 each mode versus `gaps_learned`; C4 the
three modes pairwise; C5 `gaps → ctx_<mode>` and the `ctx_*` modes pairwise (only if built). Each contrast is
evaluated on the 61 held-out loci (primary) and reported for all 119 and for the development set separately.

A contrast is called a **systematic, consistent improvement** when all of the following hold on the held-out loci:

1. Direction: one-sided exact sign test on the discordant loci (improved versus harmed), p ≤ 0.05. With 5 discordant
   loci this needs 5 of 5; with 8, 7 of 8; with 12, 10 of 12; with 20, 15 of 20 (exact binomial, α = 0.05).
2. Consistency: improved ≥ harmed within every label set (cohort Swiss-Prot, Price-148, Price Swiss-Prot) and within
   every claim set that has at least 3 discordant loci; no stratum (robust / marginal pathway, already gap / masked
   gap, archaea / bacteria) has harmed > improved + 1.
3. Minimum effect: improved − harmed ≥ 3 loci, and the harm elsewhere (3.6) does not exceed the improved count.

Outcomes: "systematic improvement" (1–3 hold); "no detectable effect" (fewer than 5 discordant loci, or 1 fails
without 1 holding in the harm direction); "inconsistent" (1 holds, 2 or 3 fails); "harm" (1 holds in the harm
direction). Secondary, descriptive: genome-block bootstrap (2,000 resamples, seed 20261102) of net improved loci;
Wilcoxon signed-rank on the rank change of E on P. Nothing is tuned on the held-out loci.

Row-5 reading: the bar ("pathway-only versus universal versus organism-specific ablations on held-out enzyme loci
with actual gaps/context") is met as a *measurement* when C3 and C4 are produced on the 61 held-out loci with the
learned stage executing (scored > 0 on them) and the outcome is stated, whatever it is; a systematic improvement is
an additional finding, not a condition of "measured". If the learned stage still abstains on these claims, the row
stays not met and the reason is named.

### 3.8 Controls

- Unmasked arms (`stage0`, `unmasked_gaps`) per genome: `stage0` must equal rc1's retained scores on the Price
  reference proteins and the unmasked cohort stage 0; `unmasked_gaps` on the 12 already-gap loci must equal the
  earlier row-5 `gaps` arm where the genome and predictor coincide.
- Sham mask (label-free): for each masked locus, a second cell masks a randomly chosen other protein of the same
  genome that has a `good_blast` row in the same pathway for a *different* reaction (seed 20261102, fixed list frozen
  with the plan). Context gains on (P, E) must not appear under the sham mask; the sham cells also give the
  distribution of harm elsewhere under a gap that is not P's.
- Invariance: for every locus and arm, prediction rows equal input claims (no enforcement), and the candidate
  universe is preserved.
- Leakage: the other session's Step L result applies to the Price claims (no near-identical reference entry; 15 Price
  Swiss-Prot labels of Sinorhizobium meliloti are circular with the reference and are reported separately here too).
  The same L1/L2 check is run for the 18 cohort genomes before any score is read (the cohort labels are Swiss-Prot
  entries of these organisms; the genus fence is expected to remove them, and the check proves it).

## 4. Use of the other session's run (plan `tasks/v3-price-contract-learned-2026-11-02.md`, live, not duplicated)

- Reused, read only, by hash: Step C claims for the 22 Price genomes; the Step L leakage result; the row-1 split.
- Not repeated: its Step F (modes with the learned stage on the *unmasked* Price gaps, 8 gap-linked loci). Its F2
  result is the mechanism check this plan depends on: if, with the learned stage executing, the three modes still
  give 0 differing scores on the whole output (possible in single-genome runs where cohort recurrence is unavailable
  and the local features are mode-invariant), then C3/C4 are null by construction for any locus set, and only §5 can
  make the modes matter. The plan then runs C1, C2 and §5 and reports C3/C4 as structurally null with that evidence.
- If its Pass I fails (learned stage does not execute on these claims), this plan waits; it does not build its own
  contract route.
- Its files and results folder (`benchmark/v3/results/price-contract-learned-v1/`) are never edited; this plan writes
  to `benchmark/v3/results/row5-masked-context-v1/` only.

## 5. Proposed GSPA method change (flagged; new method version; off by default)

**Context-weighted gap prior** ("the gap prior reads the reaction context"). Today `GapFillingPrior.beginIteration`
chooses, per gap, the single candidate with the best rescue tier and likelihood among the claims for the gap's EC/GO,
and boosts it by max(1.15 × alphaGap, 0.14 − log-odds). Proposed, behind `--gap-prior-context-mode
{none|universal|pathway|organism}` (default `none`, which is today's behaviour byte for byte):

- For gap reaction R and each candidate protein P' for R's EC within the existing window [0.18, 0.5], compute the
  local context support c(P', R) ∈ [0, 1]: the α^d-weighted fraction of R's context reactions, as `ReactionContext`
  defines them in the chosen mode (universal graph neighbours within k hops; pathway-group co-members; organism-
  restricted neighbours), that have an anchor protein (posterior > τ for the neighbour's function) within W bp of
  P'. This is the density the RLC suggester already computes (`ReactionLocalContextSuggester`, anchor weights and
  `GenomicDensityField`); the change reuses it, it adds no new feature.
- Candidate choice per gap: highest (rescue tier, then c(P', R), then likelihood) instead of (tier, likelihood).
  Boost magnitude: the existing magnitude × (0.5 + c(P', R)), so a candidate with no local context receives half of
  today's boost and one surrounded by the pathway's other enzymes receives 1.5 ×.
- Everything else unchanged: window, tiers, protein-specificity prior, the learned stage, RLC suggestions.
- Parameters (α, k, τ, W) are the RLC suggester's defaults; they are not tuned on any labelled locus. If a parameter
  must be chosen, it is chosen on the 58 cohort-train loci only, before the held-out loci are scored, and recorded.
- Version: GSPA 3.1.0-dev feature, stated as such in every table; the release-candidate JAR runs every non-`ctx_`
  arm. The `ctx_*` arms are compared against `gaps` of the same JAR build as well, so that a build difference is not
  mistaken for the method.

Expected size, said in advance from the design-time room count: on DIAMOND claims the right protein is already the
top protein for E in 46 of 47 cohort loci and 11 of 14 Price loci where it has the claim at all (median one claimer
per EC), so the choice rule can change little there; the gain, if any, is in call status (window loci) and in the
dense predictors' reaction-centric indicator. 40 of 87 cohort loci and 18 of 32 Price loci have no DIAMOND claim for
E on P; no integration-stage context can create a candidate, and these loci can only register through indicator 3 of
the dense predictors or through the RLC suggestion. This is stated so that a small measured effect is read as
expected, not as a disappointment, and so that a large one would be surprising.

Not proposed: widening the window below 0.18 (floor candidates of dense predictors are indistinguishable from
thousands of others; the Stage 2B follow-up showed the harm of touching them), or letting the modes change raw
scores outside a gap.

## 6. Compute and wall-clock (Shaheen k10018, `--exclude=nid00022,nid00497,nid01033`, outputs under
`/scratch/project/k10018/gspa-v3/out/row5-masked-context-v1/`, one package extraction per node)

Measured rates: cohort/feature pass 3.4 s per 1,000 claims; `stage0` and `gaps` 22–115 s on small cells; DIAMOND
claim files 30k–100k claims per genome; filtered dense cells ≤ 0.4 M claims.

| Block | Cells | Arms | Serial estimate | Node-hours (8 cells per 370 GB node, 24 GB heap each) |
|---|---:|---|---:|---:|
| DIAMOND, masked loci | 119 | 6 | ≈ 119 × 18 min ≈ 36 h | ≈ 5 |
| DIAMOND, sham masks | 119 | 2 (`gaps`, best mode) | ≈ 8 h | ≈ 1 |
| Six Price predictors, 32 Price loci | 192 | 2 (`stage0`, `gaps`) | small ≈ 6 h, dense filtered ≈ 11 h | ≈ 3 |
| Unmasked controls and filter check | 40 + 66 | 2 | ≈ 10 h | ≈ 2 |
| `ctx_*` arms (if §5 built), all claim sets | 119 + 192 | 3 | ≈ 24 h + 48 h | ≈ 10 |
| Evaluation | — | — | ws, minutes | — |

Total about 20 node-hours (≈ 4,000 core-hours); wall-clock under half a day once code exists, with the cells
fanned out over 10–20 nodes. Code: masking builder with audit, per-locus input builder, cell runner (extends
`run_context_cell.py`), evaluator and report, tests — one agy-auto task (`~/tmp/gspa-v3-agy/row5-masked/spec.md`),
about one day plus review; §5 is a second task (`row5-ctx-prior`), about one day plus review and Spock tests,
including a byte-identity test for `--gap-prior-context-mode none`.

## 7. Pre-registration steps (in this order)

1. Robert approves or amends this draft; §5 yes/no decided.
2. The masked-locus registry (`loci.tsv`: 119 loci with class, strata, split, sham partner) and the claim-set hashes
   are written; registry SHA-256 recorded.
3. Plan frozen (`PLAN.sha256`), then the other session's Step C/L/F results are read and §4's branch is recorded as a
   dated deviation if it applies.
4. Cells run; the evaluator is run once on the development loci first only if §5 needs a parameter choice (recorded);
   then once on the held-out loci.

## 8. Risks

- The learned stage may still abstain on parts of these claims (input-support manifest covers the contract and the
  context mode; `unseen_reaction_context_mode` is a possible reason for PATHWAY/ORGANISM if the models' manifests list
  UNIVERSAL only). Checked on the first cell; if so, C3/C4 are structurally null and reported as such.
- Single-genome runs have no cohort peers; the modes act only through the local signature. If that signature is
  mode-invariant, C4 is null regardless of loci (the §4 branch).
- Room on DIAMOND claims is small (ranking saturated; 58 of 119 loci have no claim for E on P). An honest "no
  detectable effect" is a likely outcome for C2–C4; C1 and the dense-predictor indicator 3 are where an effect can show.
- Masking optimism: the pathway prediction flag is kept; 48 loci are marginal. Reported as a stratum; criterion 2
  requires consistency across it.
- Circular labels: 15 Price Swiss-Prot labels (Sinorhizobium) and possibly some cohort labels are entries of the
  DIAMOND reference; removed by the fence or reported separately after the L1/L2 check.
- Dense predictors' claim filter changes the candidate universe; the equality check protects the labelled loci, but
  not unlabelled proteins' calls (reported as cost, not as truth).

## 9. Outputs

`docs/benchmarks/v3-row5-masked-context.md`; `benchmark/v3/results/row5-masked-context-v1/` (plan copy, PLAN.sha256,
loci registry and hash, claim-set hashes, per-cell receipts, evaluation tables, SHA256SUMS, LARGE_FILES.tsv); code in
`benchmark/v3/row5_masked_context/` with tests in `benchmark/v3/tests/`; Shaheen wrapper
`benchmark/v3/shaheen/row5-masked-context.sbatch`; if §5 is built, the engine change under `gspa-core` with Spock tests
and a release note naming the new version. `tasks/v3-acceptance.md` is not edited; proposed row text comes with the
report.

## Appendix A. Design-time numbers (2026-11-02, `~/tmp/row5-design/`)

Sign-test power (one-sided exact, α = 0.05, n discordant loci, P(improved) under the alternative):

| n | critical k | p = 0.70 | p = 0.75 | p = 0.80 | p = 0.90 |
|---:|---:|---:|---:|---:|---:|
| 5 | 5 | 0.17 | 0.24 | 0.33 | 0.59 |
| 8 | 7 | 0.26 | 0.37 | 0.50 | 0.81 |
| 12 | 10 | 0.25 | 0.39 | 0.56 | 0.89 |
| 20 | 15 | 0.42 | 0.62 | 0.80 | 0.99 |
| 30 | 20 | 0.73 | 0.89 | 0.97 | 1.00 |
| 40 | 26 | 0.81 | 0.95 | 0.99 | 1.00 |

Reading: a consistent effect (8 or 9 of 10 loci in one direction) is detectable from about 10 discordant loci; a
70/30 effect needs about 30. With 61 held-out loci and, realistically, 10–30 discordant ones per contrast, the test
detects a strong consistent effect and reports "no detectable effect" for a weak one, which is the honest outcome.

Room on the 32 maskable Price loci, six predictors (true EC a candidate but not top-1 at stage 0; in the prior's
window [0.18, 0.5]): MAPred 11 (6 in window), GraphEC 19 (1), DeepProZyme 13 (2), CLEAN 10 (1), EnzBERT 5 (1), TopEC
0. Already top-1: CLEAN 15, EnzBERT 13, MAPred 10, DeepProZyme 9, GraphEC 6. DIAMOND contract claims: 14 of 32 have
the claim, 13 already top-1; cohort DIAMOND homology: 47 of 87, 44 already top-1; family-excluded: 14 of 87.

## 10. Approval and pre-freeze amendments (2026-11-02, before any arm ran)

Robert approved the draft in full on 2026-11-02: (1) the masking design with the consistency rule of 3.7 in place of
the 20-locus gate; (2) the context-weighted gap prior of section 5 as GSPA 3.1.0-dev, off by default; (3) the 58
cohort-train loci may be used for parameter choice and for nothing else. His words: "Lock in the plans, run
implementation with agy, you monitor, control, evaluate."

Amendments made before freezing, from the other session's published result (`docs/benchmarks/v3-price-contract-learned.md`,
step F, read on 2026-11-02; it is their result, not an arm of this plan):

- A1. With gap observations present, the contract-3 learned models decline every candidate ("unseen availability
  profile": a feature-availability profile the model never saw in training), and PATHWAY and ORGANISM are "unseen
  reaction-context mode" for these models even without gaps. On every masked cell the arms `gaps_learned`,
  `universal`, `pathway` and `organism` are therefore identical to `gaps` by construction. The section-4 branch
  applies from the start: these four arms run on the pilot loci only (section 3.5, to document the equality in this
  setting), contrasts C2 to C4 are reported as structurally null with that evidence, and the three-mode ablation is
  carried by the deterministic `ctx_universal`, `ctx_pathway`, `ctx_organism` arms of section 5 (contrast C5), which
  are primary. The criterion of 3.7 is unchanged.
- A2. Pilot, fixed now: 5 masked loci (2 cohort, 3 Price, chosen by the hash order of the locus id), all arms, all
  claim sets, on one node; its receipts decide the heap classes and confirm the masking audit before the fan-out.
- A3. Parameter choice for section 5 (approval 3): the only free parameter is the context window W (bp) around the
  candidate gene; the candidate set is {5,000; 10,000; 20,000}; the rule is the largest net improved count of C5
  (`gaps` to `ctx_universal`) on the 58 development loci, ties to 10,000; anchors tau = 0.3, alpha = 0.5, k = 2 as in
  the RLC suggester. Recorded with the development-loci table before any held-out locus is scored.
- A4. Sham-mask arms: `gaps` and `ctx_universal` (the "best mode" wording of 3.8 is replaced by this fixed choice).
- A5. Results file naming: `PLAN_SHA256` (the coordinator's name) instead of `PLAN.sha256`.

## Deviations (append only, dated)
