# Six follow-up experiments on Stage 2B and off-growth enzyme discovery: pre-registered design (2026-10-31)

Status: design, frozen before any run reads a result. Changes after freezing are recorded in §9 ("Deviations") with a
reason. Requested by Robert (2026-10-31): "Design 1 to 6, all of it; then run it and report the results." The six
experiments were proposed in the Stage 2B investigation (`~/tmp/stage2b-investigation/report.md`, §4) after the
Stage 2B null result (`docs/benchmarks/v3-price-stage2b.md`) and the off-growth pilot
(`docs/benchmarks/v3-offgrowth-discovery.md`, design `tasks/v3-offgrowth-discovery-design-2026-10-30.md`).

No language model is used anywhere. Unconfirmed labels are never negatives.

## 0. Rules that apply to every experiment

**Truth and scoring.** Price-148: the 148 Price et al. (2018) proteins with 151 enzyme numbers in 22 bacteria
(`benchmark/manifests/suite6_price149_manifest.tsv`, one documented exclusion). Swiss-Prot experimental enzymes: the
35 labels on 32 proteins in 6 of the 22 genomes (`benchmark/v3/offgrowth/data/swissprot_experimental_ecs_price22.tsv`).
Price scoring reuses the frozen engine by import, never by copy or edit: `price_reference_metrics`, `ec_identity`,
`v3_metrics`, `evaluate_price_reference.read_scores/references` and the Stage 2B evaluator's `cafa_stats`,
`cafa_fmax`, `project_genome_pairs`, `summarize_corrections_and_harms` (`benchmark/v3/stage2b/evaluate_price_stage2b.py`),
with the same enzyme-identity audit (sha256 `5302494b16ecd3f76fd96d24e48e27c298206b27608be72c5c66f1d197f4f79f`), the
same 148 proteins, genome-block bootstrap (2,000 resamples, seed 13 for Price metrics as in the Stage 2B evaluator;
seed 20261030 for off-growth metrics as in the off-growth design), corrections and harms at the fixed 0.5 pair
threshold, lexical and tie-aware top-1. Off-growth scoring reuses `offgrowth.score` (recovery, reaction re-added,
harm on labelled proteins only, unconfirmed calls, exact one-sided McNemar, Holm).

**Leakage.** Nothing derived from Price-148 or Swiss-Prot labels enters any input, condition list, weight or
threshold. Every threshold and constant is fixed in this document; none is fitted. The Price assay substrate of a
protein (`cognate_exchange_rxn`) is never read. Conditions are per genome, identical for every protein of the genome,
and derived only from label-free, organism-level sources (Biolog plate maps, the genome's own transporter calls).
Labels are read only by the scorers. Each run records the sha256 of this design file in its receipt.

**Baseline and controls.** Every experiment carries the unchanged baseline (today's rule or today's production
output) as an arm, and every experiment that touches labels carries a shuffled/permuted control (fixed seeds 1–5).

**Compute.** Heavy work runs on Shaheen (account k10018, workq, full exclusive nodes, 192 worker processes); outputs
under `/scratch/project/k10018/gspa-v3/out/stage2b-followups-v1/`, per-genome outputs packed in tars. Code lives in
`benchmark/v3/followups/` with plain unittest tests and `benchmark/v3/shaheen/followup-*.sbatch`.

**Reporting.** `docs/benchmarks/v3-stage2b-followups.md` (numbers in tables, plain-word narrative); small result files
in `benchmark/v3/results/stage2b-followups-v1/` with sha256 sums. A failed criterion is reported as "not met" with the
numbers; nothing is re-defined after the fact.

## 1. Competitive gene assignment for off-growth calls

**Question.** When a targeted test re-adds the true reaction, does a competitive assignment rule give it to the hidden
protein more often than "top-scoring candidate" does, without more harm?

**Hypothesis.** The Stage 2B responsibility rule with a synteny covariate, an explicit null gene and paralog
competition raises correct assignments among re-added truth reactions from 23 of 54 to at least 35 of 54.

**Data (stored, no new linear programs).** The off-growth pilot outputs on Shaheen
(`offgrowth-v1/genomes/<genome>.outputs.tar` and `<genome>.swissprot.outputs.tar`: per case, per arm, per test the
added reactions and the test identifier), the union candidate tables (`offgrowth-v1/candidates/`), genome layouts
(`price-full-matrix-v2/layout/<genome>.layout.tsv`: protein, contig, start, end, strand), gapsmith transporter tables
(`benchmark-models-v1/.../<genome>-Transporter.tbl`), KEGG module table (`price-production-v1/reference/kegg_modules.tsv`),
SEED reaction and enzyme-alias tables (`gapsmith-reference-bacteria-v1/reference/`), Biolog panel
(`benchmark/v3/offgrowth/data/biolog_pm1_pm4_panel.tsv`). Cases: the 178 combined cases (148 Price + 32 Swiss-Prot);
the arm of record is A1 ∪ A2 ∪ A3 (union), as in the pilot.

**Arms (assignment rules; the added reactions are identical in all arms).**

| Arm | Rule for an added reaction r with candidate genes H(r) (every candidate with an EC mapped to r) |
|---|---|
| R0 baseline | top-scoring candidate (today's rule, `assign_calls_for_added_reactions`) |
| R1 responsibility | weight w(h) = odds(s_h)^1 · exp(1.5 · C_syn(h, r)); assign to argmax w if w_max > τ_null, else "no gene"; τ_null = 0.05, β1 = 1, β3 = 1.5 (the production constants of `LocusResponsibilityUpdater`) |
| R2 responsibility, stricter null | as R1 with τ_null = 1.0 (the null gene is as plausible as a candidate with score 0.5) |
| R3 responsibility + paralog competition | R1, then: within one arm's call set, a gene may win at most one reaction per enzyme number; reactions are processed in increasing reaction-id order and a gene already holding a reaction of the same enzyme number is skipped in favour of the next candidate by w |
| R4 synteny only | R0 restricted to candidates with C_syn = 1 when any exists, otherwise R0 |

C_syn(h, r) = 1 when, within ±10 loci of h on the same contig (layout order), there is another gene g ≠ h that either
(a) carries a candidate claim (score ≥ 0.25) for a different enzyme number of the same KEGG module as any enzyme number
of r, or (b) is a gapsmith transporter hit for the sole-source substrate of the test that added r (A1 tests only; the
substrate is the test's panel compound). Otherwise 0. The 0.25 floor is the candidate cut already used by the
production `reactions.tsv` export; it is fixed here, not tuned.

**Metrics.** Correct assignments among re-added truth reactions (primary; denominator = cases whose true reaction was
re-added by the union arm = 54), recovery over all 178 cases and the 146 + 32 off-growth cases with genome-block
bootstrap, gene calls, supported / harmful / unconfirmed calls and harm rate on labelled proteins (Price-148 +
Swiss-Prot), "no gene" calls, paired exact one-sided McNemar R1 vs R0 (recovered cases).

**Controls.** Shuffled candidate genes (seeds 1–5) for R1 and R3: recovery must be near zero. Shuffled layout
(gene positions permuted within genome, seeds 1–5) for R1: C_syn must stop helping (recovery falls to R0 ± 2 cases).

**Pre-registered criterion.** Met if R1 or R3 (declared in advance: R1 is primary, R3 secondary) assigns ≥ 35 of the 54
re-added truth reactions to the hidden protein AND its harm rate on labelled proteins is ≤ 0.85 (the union's pilot
value) AND the shuffled-gene control recovers ≤ 2 of 178. Partial: 24–34 correct is reported as an improvement that
misses the bar. Any R arm with harm above 0.85 is reported as harmful.

**Cost.** One Shaheen node, minutes. **Risks.** Paralog clusters in one operon stay ambiguous; expected gains on
near-tie cases only.

## 2. Condition-free targeted tests as a gated evidence source for Stage 2B

**Question.** Can the unblocking (A2) and dead-end (A3) tests, which need no condition list, change integrated enzyme
probabilities for the better when gated by the assignment rule of §1 and a pathway-coherence check?

**Hypothesis.** Reactions that a model must add to let gene-evidenced neighbours carry flux or to drain dead-end
metabolites are usable evidence for the gene that fills them, when gated.

**Data.** The unmasked A2 and A3 runs of the pilot (stored per genome in `genome_results.json`,
`unmasked_arm_results`: every test, its target and its added reactions; label-free: union candidate table and gapsmith
draft of the genome). Per predictor, the stage-0 integrated output of the Stage 2B collection
(`price-stage2b-rc1/collection/attempts/<cell>/.../production/stage0/<genome>.integrated.tsv`, verified by the
collection finalizer). Layouts, KEGG modules, SEED tables as in §1. Gene–reaction rules of the gapsmith draft model.

**Ceiling, stated in advance.** In the unmasked model the truth reaction is present whenever gapsmith assigned it, so
only cases whose true reaction the unmasked union run re-adds can be reached (the pilot's unmasked control re-added
the true reaction in 10 of 178 cases). The +8 top-1 bar below is therefore expected to be unreachable; the
experiment measures what the evidence source does, including harm, not whether it meets a bar it cannot meet.

**Procedure (benchmark-level implementation of the proposed sidecar mode).** For every genome and every added
reaction r of an unmasked A2 or A3 test: candidates H(r) = proteins with a stage-0 claim (score > 1e-5, today's
candidate floor) for an enzyme number mapped to r; responsibility r(h) by the §1 R1 rule using the predictor's own
stage-0 scores; gate G1: r(h) exceeds the null share τ_null / (τ_null + Σw); gate G2 (pathway coherence): r shares a
non-currency metabolite (currency set = the 20 most connected metabolites of the universal SEED table, listed in the
run receipt) with a reaction already assigned by the draft model's gene–reaction rules to a gene within ±10 loci of h.
Update for each gated (h, r): Δlogit = γ_FBA · π_j · r(h) · (1 − p_j) with γ_FBA = 2.2 (production), π_j = 1
(every added reaction of an A2/A3 test carries the demanded flux by construction), p_j = max stage-0 score among H(r)
for the enzyme numbers of r (reaction probability before the evidence). The update is written into a copy of the
stage-0 output (`posterior_prob` only), as `LocusResponsibilityUpdater` would. Arms: stage 0 (baseline), stage 0 +
A2 evidence, stage 0 + A3 evidence, stage 0 + A2 + A3 evidence, and stage 0 + A2 + A3 without gate G2 (to show what
the gate removes). Production wiring (a `--fba-targets` sidecar mode) is deferred until the criterion is met; the
update formula and constants are the production ones.

**Metrics.** Price scoring (§0) per predictor: top-1, top-3, mean reciprocal rank, micro and CAFA F-max with genome
bootstrap, corrections and harms (pair threshold 0.5, lexical top-1, tie-aware top-1), McNemar. Harm on labelled
proteins: harmful calls among gated updates that move a labelled protein's top enzyme to a wrong one (Price-148 +
Swiss-Prot). Held-out Swiss-Prot cases: the 32 proteins scored separately (top-1 before and after).

**Controls.** Shuffled gene identities (seeds 1–5) in the evidence step: corrections must vanish. Stage 0 → stage 0
no-op check (zero corrections and harms).

**Pre-registered criterion.** Met if, for at least four of six predictors, top-1 on the 148 Price proteins rises by
≥ 8 over stage 0 (one-sided exact McNemar p < 0.05), the harm rate among gated updates on labelled proteins is < 0.3,
and CAFA F-max does not fall for any predictor. Failure is reported with the reachable-case ceiling and the harm rate
first.

**Cost.** One node, under an hour (streaming 132 stage-0 outputs). **Risk.** Pilot harm for A2/A3 was 0.77–0.80.

## 3. Condition-matched media without labels: a fixed per-genome condition panel

**Question.** Can the sole-source panel be selected per genome from organism-level, label-free data, keeping the
recovery of the full Biolog panel at a fraction of its tests?

**Hypothesis.** Substrates for which the genome carries a transporter call are the ones whose catabolic enzymes a
sole-source test can expose; a panel restricted to them keeps the full panel's recovery at ≤ 25 % of its tests.

**Data.** Stored A1 per-test results of the pilot (masked runs, 178 cases): each test's panel compound and added
reactions. Genome-specific panels are built before any scoring from (i) the Biolog PM1–PM4 panel (293 usable wells,
`biolog_pm1_pm4_panel.tsv`, sha256 recorded) and (ii) the genome's gapsmith transporter table (substrate compound ids
of transporter hits). BacDive strain growth records were planned as a second organism-level source; the repository's
BacDive manifest covers none of the 22 species and the BacDive API needs credentials that are not available offline,
so BacDive is not used (recorded here, not a deviation). No new linear programs: arm P1 is a subset of the stored A1
tests.

**Arms.** P0 = full Biolog panel (the pilot's A1, 293 wells); P1 = per-genome transporter panel (wells whose SEED
compound has a transporter hit in the genome's gapsmith table; if a genome has fewer than 10 such wells, all its
wells are kept and the genome is flagged); P2 = random panel of the same size as P1 per genome (seeds 1–5, the
size-matched control). Panels are written and hashed (`panels.json`, sha256 in the receipt) before scoring.

**Held-out halves.** No threshold is fitted, so there is no calibration split. As a stability check, the wells are
split by well-column parity (odd / even columns) and P0, P1 are reported on each half as well as on the whole panel.

**Metrics.** Recovery (hidden protein assigned; off-growth stratum primary) and reaction re-added, per arm, for the
178 cases and the 148 Price cases, genome-block bootstrap; tests executed per arm (fraction of P0); harm rate on
labelled proteins; P1 vs P0 paired McNemar (two-sided exact, because P1 ⊂ P0 can only lose).

**Controls.** P2 size-matched random panels (seeds 1–5): P1 must recover more than the mean of P2. Shuffled candidate
genes (seeds 1–5) for P1: near zero.

**Pre-registered criterion.** Met if P1 recovers ≥ P0's recovery (178 cases, union of the P1 tests; equality allowed)
using ≤ 25 % of P0's tests, and P1 recovers ≥ 9 of 148 Price cases (the pilot's A1 value). Reported as not met with
the P1/P0 ratio otherwise.

**Cost.** Minutes on one node (re-scoring). **Risk.** Transporter calls are incomplete (gapsmith's transporter table);
P1 may drop the sugar-phosphate and nucleoside wells that drove the pilot's recoveries.

## 4. Give Stage 2B its covariates

**Question.** The responsibility rule received empty covariate maps in all 132 cells. With a synteny covariate
computed from the genome layout and KEGG modules, does Stage 2B's attribution change, and does it change the Price
metrics?

**Hypothesis.** The empty maps are a configuration defect; with a synteny covariate present the responsibility
concentrates on genes with pathway neighbours, and the Price top-1 changes are reported with tie-aware credit. The
expected Price effect is small (the Price enzymes stay off the growth path); this is a correctness fix.

**Why not GO claims and the GO priors.** The investigation proposed deriving GO claims from the EC claims so that the
genomic-context and coherence priors fire. The Price cells hold up to 45 million EC claims each; expanding them to GO
terms would change the candidate universe the collection finalizer verified and multiply the cell inputs. Instead the
covariate the sidecar expects (`operon_partner_boost_by_gene_ec`) is computed directly for EC claims from the layout
and the KEGG module table, by the sidecar itself behind a new option, so the production path is exercised end to end
except for the JVM launch.

**Data.** For each of the 132 Stage 2B cells: the stored sidecar payload (`fba_input.json`: all gene–EC scores) and
result (`fba_output.json`: repaired reactions, certificates π), the stage-0 output, the genome layout, the KEGG
module table.

**Implementation (production code, behind options whose default keeps today's behaviour).**
`scripts/stoichiometric_locus_refiner.py` gains `--repair-json <file>` (reuse a stored repair result instead of
solving), `--layout-tsv <file>`, `--kegg-modules-tsv <file>`, `--synteny-window 10` (compute C_syn as in §1 (a) for
every (gene, EC) pair with score ≥ 0.25 and write it into `operon_partner_boost_by_gene_ec` when the payload's map is
empty), and the §5 options. The Groovy `StoichiometricLocusRefiner` / `IntegrateCommand` pass new options
`--fba-synteny-from-layout` (uses the run's `--genome-layout`), `--fba-kegg-modules`, `--fba-candidate-min-score`,
`--fba-candidate-top-k` through to the sidecar; defaults off. Full JVM gate and a release-candidate-3 JAR with sha256
and source snapshot staged on Shaheen. The cells are re-attributed by running the sidecar script on the stored
payloads (no new linear programs; the repair result is the stored one), and the updated scores are written into a copy
of the stage-0 output exactly as `LocusResponsibilityUpdater.applyUpdatedScores` does (`posterior_prob` replaced where
it changes by more than 1e-5).

**Arms.** S0 = stage 0; S2B = stored stage 2B output (today); S2B-syn = re-attribution with the synteny covariate;
S2B-syn-k = synteny covariate plus the §5 candidate limit (min score 0.05 and top 10 per enzyme number).

**Metrics.** Covariate counts per cell (non-empty entries), updated loci per cell, Price scoring (§0) for S2B-syn and
S2B-syn-k against S0, and against S2B (paired), per predictor.

**Controls.** Shuffled layout (seed 1): the covariate must change (count differs) and the Price result must not
improve over S2B-syn. Re-attribution with the empty maps must reproduce the stored stage 2B output to within 1e-6
(replay check, all 132 cells).

**Pre-registered criterion.** Met (as a correctness fix) if the covariate map is non-empty in every cell, the replay
check passes in all 132 cells, and no predictor's Price top-1, top-3 or CAFA F-max is worse than S2B's. A Price gain
is not required; any gain is reported with bootstrap intervals.

**Cost.** Payloads are up to 3.4 GB JSON per cell; 132 cells on about 8 nodes, about 2 node-hours.
**Risk.** None beyond compute.

## 5. Stop diluting responsibility for dense predictors

**Question.** For dense predictors (CLEAN, DeepProZyme, GraphEC) every protein with score > 1e-5 is a candidate, so
each repaired reaction's responsibility is spread over thousands of genes and the boosts are tiny tie re-orderings.
Does a candidate limit concentrate the boost on real candidates without harming the Price result?

**Hypothesis.** Limiting candidates per repaired enzyme number to score ≥ 0.05 or the top 10 removes the tie
artefacts and leaves the Price metrics unchanged or better.

**Implementation.** Sidecar options `--candidate-min-score` (default 1e-5 = today) and `--candidate-top-k` (default 0 =
unlimited); candidates below the floor or outside the top k (ties at the k-th score kept) get no responsibility and
no update. Same re-attribution procedure as §4 on the stored payloads.

**Arms.** S0; S2B (stored); K1 = min score 0.05; K2 = top 10 per enzyme number; K3 = both.

**Metrics.** Scores changed per cell; top-call changes at the predictor's floor score (tie artefacts: a protein whose
top enzyme changes while its top score equals its floor score before and after); Price scoring per predictor against
S0 and against S2B.

**Controls.** Replay check (empty options reproduce S2B); the sparse predictors (MAPred, TopEC, EnzBERT) must be
unchanged by K1–K3 apart from candidates genuinely below the floor (reported).

**Pre-registered criterion.** Met if, for every dense predictor, changed scores per cell are < 1,000 in every cell,
tie-artefact top-call changes are 0, and no Price metric (top-1, top-3, CAFA F-max) is worse than S2B's. Reported as
not met otherwise.

**Cost.** Shares §4's node-hours. **Risk.** None.

## 6. RB-TnSeq fitness as a held-out test of condition-specific essentiality predictions

**Question.** Do single-gene knockouts of the genome-scale model on the fixed condition panel predict which genes have
a specific fitness defect on that condition in the public RB-TnSeq data, and does the competitive gene assignment of
§1 improve that prediction over the baseline assignment?

**Hypothesis.** The knockout growth ratio predicts specific phenotypes better than chance (AUC > 0.5 per genome), and
models whose added reactions are assigned by the §1 rule score higher than models assigned by the top-scoring rule.

**Fitness data (test only, never an input).** The public Fitness Browser database (`feba.db`, SQLite, downloaded
from https://fit.genomics.lbl.gov/cgi_data/feba.db on 2026-10-31; sha256 recorded) and the protein sequences
(`aaseqs`, same site). Only data are downloaded; no binaries. Tables used: `Organism`, `Gene`, `Experiment`,
`GeneFitness`, `SpecificPhenotype`. Fitness Browser genes are matched to the 22 proteomes by exact protein-sequence
sha256 against `price-full-matrix-v2/layout/<genome>.sequence_hashes.tsv`; the organism ↔ assembly mapping is the
22-row table written to `benchmark/v3/followups/data/price22_fitness_browser_organisms.tsv` before any scoring, and a
genome is used only if ≥ 80 % of its model genes match a Fitness Browser gene. If the database cannot be obtained or
the match rate fails, the experiment stops and says so.

**Conditions.** The fixed Biolog panel (§3 P0; same file, same hash). A Fitness Browser experiment is matched to a
panel well when its group is a carbon, nitrogen, phosphorus or sulfur source and its condition name maps to the
well's SEED compound by the panel builder's name normalisation (`offgrowth.panel.generate_candidate_name_keys`); the
match table is written and hashed before scoring. Experiments that match no well are listed, not used.

**Predictions.** For each genome and each model variant: for each matched condition, the sole-source medium of the
off-growth design (`build_sole_source_medium`); wild-type growth; if the wild type grows (> 1e-6), single-gene
knockouts of every model gene; prediction score for (gene, condition) = 1 − growth(knockout) / growth(wild type).
Genes without fitness data on that condition are not scored. Model variants: M0 = gapsmith gap-filled model as is
(its own gene–reaction rules); M1 = M0 plus the unmasked union (A1 ∪ A2 ∪ A3) added reactions assigned by the
baseline rule R0; M2 = the same reactions assigned by the §1 R1 rule; M3 (control) = the same reactions assigned to
random model genes (seeds 1–5). Media and solver as in the off-growth design (GLPK through COBRApy; replay not
needed because no gap-filling is solved).

**Truth.** A (gene, condition) pair is positive when the Fitness Browser lists a specific phenotype for the gene in
an experiment matched to that condition (`SpecificPhenotype`, Price et al. 2018 definition) and the fitness is
negative (defect); negative otherwise, among genes with fitness data in that experiment. (Positive fitness in a
condition is not a growth requirement and is treated as negative here; stated as a limit.)

**Metrics.** Per genome and variant: AUC of the prediction score over (gene, condition) pairs restricted to model
genes, and restricted to genes of the added reactions (the only pairs where M1–M3 differ); number of pairs,
positives; pooled AUC across genomes with genome-block bootstrap; precision at the top 50 predicted pairs per genome.

**Controls.** M3 random assignment (seeds 1–5); a permutation of fitness labels within each genome (seeds 1–5) must
give AUC ≈ 0.5.

**Pre-registered criterion.** Met if the pooled AUC of M0 over model genes is > 0.5 with the bootstrap interval
excluding 0.5 in at least 11 of the usable genomes, and the pooled AUC of M2 over added-reaction genes exceeds M1's
with the paired bootstrap interval excluding 0. The first part validates the knockout readout; the second validates
the assignment rule. Each part is reported separately.

**Cost.** About 22 genomes × up to 293 conditions × about 1,500 model genes × 4 variants ≈ 50 million linear programs
upper bound; only conditions where the wild type grows and only matched conditions are scored, so roughly 2–5 million;
one node per genome, well under an hour each. **Risks.** Media and biomass assumptions dominate; condition name
matching may miss experiments; the fitness libraries lack insertions in essential genes, so default-medium
essentials cannot be scored.

## 7. Execution order and receipts

1. Freeze this file; record its sha256 in every run receipt (`design_sha256`).
2. Implement `benchmark/v3/followups/` and the sidecar/Groovy options; run the unit tests; run the full JVM gate;
   build the JAR; stage release-candidate-3 on Shaheen with sha256 and source snapshot.
3. Shaheen jobs: `followup-rescore.sbatch` (§1, §3; one node), `followup-targets.sbatch` (§2; one node),
   `followup-stage2b-rerun.sbatch` (§4, §5; array over cells), `followup-fitness.sbatch` (§6; one node per genome),
   `followup-score.sbatch` (Price scoring of §2, §4, §5 variants).
4. Report.

## 8. Deliverables

Code and tests under `benchmark/v3/followups/`; job scripts `benchmark/v3/shaheen/followup-*.sbatch`; the report
`docs/benchmarks/v3-stage2b-followups.md`; result files in `benchmark/v3/results/stage2b-followups-v1/` with
`SHA256SUMS`; progress entries in `tasks/v3-supervision-plan-2026-10-28.md`.

## 9. Deviations

Recorded after freezing (the file's frozen hash is the one above; this section is the only change). None alters a
criterion, a threshold or an input; they concern reporting and implementation.

1. **Case count (§1, §3).** The stored off-growth outputs hold 180 cases (148 Price + 32 Swiss-Prot); the pilot report
   quoted 178 for its combined report. All 180 stored cases are re-scored; the baseline rule reproduces the pilot's
   23 recoveries and 54 re-added truth reactions exactly.
2. **Experiment 6 readout (§6).** The gapsmith gap-filled model (M0) grows on none of the sole-source wells, so the
   pre-registered readout on M0 cannot be evaluated. The aggregate additionally reports the same per-genome and
   pooled readouts for M1 and M2 (post hoc, labelled as such); the M0 criterion is reported as not testable, not
   re-defined.
3. **Experiment 6 bootstrap (§6).** The pooled genome-block bootstrap and the per-genome pair bootstrap were
   re-implemented with weighted rank sums (numpy) after the first aggregate job proved too slow; the estimates are
   identical to the pure-Python implementation on the two genomes checked.
4. **Experiment 4 (§4).** The sidecar command-line path was exercised on the 22 MAPred cells (small payloads) and the
   remaining 110 cells use the same functions imported in process; results agree exactly where both were run.
5. **Experiment 2 reporting (§2).** The harm rate "among gated updates on labelled proteins" counts updated
   (protein, enzyme number) pairs, as implemented, so one labelled protein with several updated enzyme numbers counts
   several times; the Price-level corrections and harms are reported alongside.
