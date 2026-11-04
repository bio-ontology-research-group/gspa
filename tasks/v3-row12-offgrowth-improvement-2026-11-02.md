# Off-growth enzyme discovery, second round: recover some additional enzymes and pathways at measured precision

Status: **approved and frozen** (2026-11-02). Robert approved the plan in full on 2026-11-02: (1) the Fitness
Browser reannotations as the primary held-out label set (tier A, label-level hold-out on seen genomes, stated
honestly); (2) tier B, including GPU time for predictor inference on the eight cohort proteomes; (3) tier C, the 17
Fitness Browser organisms outside the 145 benchmark genomes. "Lock in the plans, run implementation with agy, you
monitor, control, evaluate." The frozen sha256 of this file is recorded in
`benchmark/v3/results/offgrowth-improvement-v1/PLAN_SHA256` before any tier A/B/C score is read; later changes go
to §13 (dated, append-only). Requested by Robert (acceptance row 12):
"we only need to recover SOME additional EC functions and pathways this way, not all of them."
Exploration that motivates the design: `docs/benchmarks/v3-offgrowth-improvement-exploration.md`.
Previous round: `tasks/v3-offgrowth-discovery-design-2026-10-30.md`, `docs/benchmarks/v3-offgrowth-discovery.md`,
`docs/benchmarks/v3-stage2b-followups.md`. No language model is used anywhere; all rules are deterministic.

## 1. Claim under test

Targeted growth simulations on a gapsmith model, with a specificity gate and a gene-assignment rule that abstains,
recover **some** experimentally supported enzymes and pathway steps that neither sequence-based annotation at its
operating point nor biomass gap-filling delivers, at a labelled precision that is measured on held-out labels and
useful (pre-registered level in §8), and they abstain elsewhere.

## 2. What has been seen (honest statement)

* Seen and used for every design constant below: the pilot outputs on the 22 Price genomes and their outcomes
  against Price-148 and the 35 Swiss-Prot experimental labels (exploration §1–§6).
* Seen but **not scored**: the 615 Fitness Browser EC reannotations on the 22 genomes (exploration §7; file
  written, outcomes never computed), the 79 Swiss-Prot experimental enzyme labels of the eight recurrence-cohort
  bacteria.
* The 22 genomes' models, candidate tables, Biolog panel and the simulations themselves were seen; held-out tier A
  (§5) is therefore held out at the **label** level only. Tier B is held out at the genome level.
* Constants (F, δ, s_min, W, k) are chosen in §4 from the exploration and are fixed here; nothing is fitted on tier
  A or B. The fixed-precision rule of §8 earns its threshold on the seen labels only.

## 3. Data

* Models: gapsmith draft and gap-filled models of the 145 benchmark genomes
  (`/scratch/project/k10018/gspa-v3/benchmark-models-v1/`, October recipe); only bacterial genomes are used.
* Candidate tables, two kinds:
  * **six-predictor** (tier A only): as in the pilot (union of CLEAN, DeepProZyme, EnzBERT, GraphEC, MAPred, TopEC
    stage-1 posteriors, best score per protein and enzyme, dense-table cap unchanged). Five of the six predictors
    exist only as a collaborator's prediction snapshot for the 22 Price proteomes and cannot be run on new genomes.
  * **own-CLEAN** (all tiers): CLEAN (split100 supconH checkpoint, ESM-1b mean embeddings, max-separation rule) run
    by us on orix for the 22 Price, 8 cohort and 17 Fitness Browser proteomes; score = CLEAN's distance-based
    confidence in [0, 1] as produced by the repository's `clean_contract` (`gmm_confidence`), one row per (protein,
    EC) with confidence > 0. Tier A is run with both tables so that the own-CLEAN results of tiers B and C can be read
    against a tier-A reference on the same scale. The §4 constants are applied unchanged to both tables; the
    scale difference is a declared limit, not a tuning opportunity.
* Truth and harm labels (complete four-level ECs only; unconfirmed labels are never negatives):
  * L0 (seen): Price-148, Swiss-Prot experimental labels of the 22 genomes.
  * L1 (tier A, unseen outcomes): Fitness Browser `ReannotationEC` rows joined by exact sequence to the 22
    proteomes, minus L0 proteins: 615 proteins. Label tier: curated reannotation from mutant phenotypes plus
    homology (Price et al. 2018 and later Fitness Browser releases), recorded per row with the `Reannotation`
    comment. Built by a new loader with a receipt (`feba.db` sha256 `e033a1b1…`, query, join counts).
  * L2 (tier B, unseen genomes): Swiss-Prot experimental ECs (ECO:0000269 and descendants, exact sequence) of the
    eight recurrence-cohort bacteria: 79 labels, from `training-cohort/labels/labels.tsv`. These labels trained the
    learned integrator; the candidate tables used here come from the predictors' raw scores with the heuristic
    integration, which that model does not touch; stated as a limit regardless.
  * L3 (tier C, unseen genomes, approved): 293 Fitness Browser EC reannotations (complete ECs) in 17 further
    Fitness Browser organisms (HerbieS, Keio, Miya, Koxy, Putida, Caulo, Cola, Dyella79, Pedo557, Burk376, Btheta,
    DvH, SynE, Bifido, rhodanobacter_10B01, Xantho, Lysobacter_OAE881). Proteomes are the Fitness Browser protein
    sequences (`aaseqs`, protein-coding genes of the `Gene` table), genome layout from the `Gene` table
    (scaffold, begin, end), gapsmith draft and gap-filled models built on Shaheen with the October recipe
    (bacterial reference), own-CLEAN candidate tables (above). Keio is *E. coli* BW25113; its 5 labels are kept.
* Measured growth conditions: the Fitness Browser carbon/nitrogen/phosphorus/sulfur experiments matched to Biolog
  wells (`benchmark/v3/results/stage2b-followups-v1/experiment_panel_matches.tsv`, 17 of 22 genomes); the
  existence of an experiment on compound X is a measured positive growth phenotype. BacDive is not available
  (credentials). Biolog PM1–PM4 mapping as in the pilot (293 wells).
* Reaction–enzyme map, KEGG modules, genome layout (locus order from the proteome FASTA): as in the follow-ups.

## 4. The improved deterministic procedure

Starting point per genome: the masked (or, for the deployment run, unmasked) draft model D and universal set U with
the pilot's weights (1 for a reaction with a candidate gene, 100 otherwise).

**Step 1, baseline minimal fill (new).** Gap-fill D once for biomass on the **minimal medium** M_min: the gapsmith
default medium (`<genome>-medium.csv`) with every organic compound removed except D-glucose (cpd00027) and the two
vitamins thiamin (cpd00305) and riboflavin (cpd00220); water, protons, ammonium, phosphate, sulfate, H2S and the
inorganic ions keep their default uptake bounds. Lactate, the twenty amino acids and chorismate are removed. Sole-source
tests replace, in M_min, glucose (C tests), ammonium (N tests), phosphate (P tests) or sulfate and H2S (S tests) by the
panel compound at the pilot's uptake rate. The added set R_base is recorded, never called, and is shared by every test of
the genome. Model M_base = D + R_base.

**Step 2, tests.**
* T1 differential sole source: for each panel compound, growth on M_base with that compound as the only C (or N, P,
  S) source; gap-fill where growth fails; the fill R_test excludes R_base by construction.
* T2 unblock and T3 dead end: unchanged from the pilot, run on D under the union medium.
* Each fill is solved twice (GLPK and HiGHS through scipy, as in the pilot's replay); a reaction counts as
  **stable** when both solutions contain it.

**Step 3, reaction-level calls (new gate).** A reaction r is called for the genome when some test adds it with a
stable fill of at most **F = 3** reactions (sensitivity: F = 1 and F = 5, reported). Each call records the test,
the fill, and the number of distinct tests that add r. Reactions already carrying gene evidence in D are never
called.

**Step 4, gene assignment with abstention (new rule).** For a called reaction with enzyme numbers E(r), the
candidate protein p* with the highest score over E(r) is named only if all of the following hold:
* score(p*) ≥ **s_min = 0.30**;
* margin over the second-best protein for the same enzyme ≥ **δ = 0.05**, or no second candidate at ≥ s_min;
* p* is not already assigned in D to a reaction whose enzyme it scores higher than E(r) ("free protein"), or E(r)
  is p*'s best-scoring enzyme (mutual best).
Ties within 0.02 of p* that pass the same conditions are co-named (isozymes; at most 2 proteins per reaction).
Otherwise the call stays **reaction-only** (gene abstention) with the ranked candidate list attached.

**Step 5, pathway calls (new, reported separately).** Within one test's fill, reactions that share a non-currency
metabolite form a path; if at least two reactions of a path have named or candidate genes within **W = 10** loci
of each other, the path is a pathway call with its gene cluster. A hidden case is pathway-recovered when its true
reaction is in a pathway call and its protein is in the cluster.

**Step 6, panels.** Two panel arms: P_T, wells whose compound has a gapsmith transporter hit in the genome
(label-free, as in the follow-ups); P_M, wells on which the organism was measured to grow (Fitness Browser
conditions; genomes without data fall back to P_T and are flagged).

Everything above uses only masked, label-free inputs; the masking audit of the pilot is re-run on the new inputs
(baseline fill, panel, stability replay).

## 5. Held-out design

| Tier | Genomes | Hidden cases (truth) | Harm labels | Held out at | Role |
|---|---|---|---|---|---|
| A | the 22 Price genomes | L1: 615 Fitness Browser reannotations (off-growth stratum expected > 95 %) | L0 + L1 | label level | **primary** |
| A0 (paired re-run) | the 22 Price genomes | L0: 180 pilot cases, re-run under the new procedure | L0 + L1 | seen | old-vs-new paired comparison |
| B | 8 recurrence-cohort bacteria | L2: 79 Swiss-Prot experimental ECs | L2 | genome level | confirmatory (own-CLEAN tables) |
| C | 17 Fitness Browser organisms | L3: 293 | L3 | genome level | confirmatory (own-CLEAN tables); pooled with B as the genome-level verdict |

Tier A is run with the six-predictor tables (primary) and the own-CLEAN tables (reference for B and C).

Every constant is fixed before any tier A/B outcome is read. The fixed-precision rule (§8) is applied once, on
L0, before tier A is scored. Tier A genomes are reported both pooled and with the seven Pseudomonas-group genomes
separated from the other fifteen, because about 350 of the 615 labels and the repeated enzyme families lie in that
group.

## 6. Arms (equal simulation budget)

The budget unit is linear programs per genome. The new procedure's budget B_new (baseline fill + P_T tests, two
solves each) defines the common budget; arms that would exceed it are cut by the pilot's fixed rule order.

| Arm | What it does | Budget |
|---|---|---|
| S0 standard gap-filling | biomass fill on the default medium, top-scoring assignment (pilot A0) | own (tiny) |
| S1 old procedure, full | pilot A1 ∪ A2 ∪ A3, top-scoring assignment, no gates | pilot (about 5,700 LPs/case) |
| S1b old procedure at budget | S1 tests cut to B_new by the pilot's rule order | B_new |
| N1 new, transporter panel | §4 with P_T | B_new |
| N2 new, measured-growth panel | §4 with P_M (17 genomes; others as N1, flagged) | ≤ B_new |
| N3 new, random panel | §4 with size-matched random wells, seeds 1–5 | B_new |
| N1−gate, N1−margin, N1−stable | ablations: without the fill gate (F = ∞), without the margin/free rule, without the two-solver stability | B_new |
| Q sequence-only, volume-matched | per genome, call (p, E) when E is p's best enzyme and score ≥ t, with t set so that the number of calls equals N1's gene calls in that genome; also the fixed-precision threshold version | no LPs |
| R shuffled genes | N1 with candidate genes permuted within genome, seeds 1–5 | B_new |

Every arm is run per candidate table (six-predictor on tier A; own-CLEAN on tiers A, B, C). Q is matched to the
arm with the same table.

## 7. Metrics

* **Recovery** (masked, per hidden case, off-growth stratum): gene-level (hidden protein named for a reaction of
  its enzyme), reaction-level (true reaction called), pathway-level (§4 step 5). Genome-block bootstrap intervals;
  paired exact McNemar between arms.
* **Precision and harm** on **distinct labelled decisions** of one unmasked run per genome: a gene call on a protein
  with a label is supported if the called enzyme maps to the same reaction set as a label, harmful otherwise
  (related-substrate harm, first three EC levels equal, is reported as a sub-row but counts as harm). Unlabelled
  calls are unconfirmed. Bootstrap over genomes.
* **Abstention**: share of called reactions left reaction-only; share of candidate reactions not called at all.
* **Volume**: gene calls, reaction-only calls, pathway calls per genome.
* **Cost**: LPs, gap-fills, CPU hours, wall time per genome and per case; fill-size and stability distributions.
* **Additional functions** (the row-12 question): gene calls of N1/N2 on labelled proteins whose called enzyme is
  **not** called by Q at its fixed-precision threshold (or Q abstains), with their precision; reaction-level calls
  whose enzyme has no candidate at ≥ s_min anywhere in the genome ("orphan pathway steps") are counted and listed.

## 8. Pass criteria (pre-registered once frozen)

Fixed-precision level: **0.75** on at least **20** labelled decisions (the pilot's 0.90 was never reached by any
arm with ten decisions; 0.75 is the level at which a call is worth a curator's look). The rule is earned on L0:
an arm whose L0 gene calls reach 0.75 on ≥ 20 decisions may call on tier A; otherwise it abstains on tier A and
the result is negative for that arm. Everything below is evaluated on tier A unless stated.

1. **Measured precision.** N1 (or N2) gene calls on labelled proteins reach precision ≥ 0.75 with ≥ 20 labelled
   decisions in the unmasked run, harm ≤ 0.25, bootstrap lower bound ≥ 0.60.
2. **Some additional recovery.** N1 gene-recovers ≥ 25 of the tier-A hidden cases (about 4 % of ≈ 600 off-growth
   cases), S0 recovers 0, and ≥ 10 of them are enzymes that Q does not call for that protein at its fixed-precision
   threshold.
3. **Better than the old procedure for the same simulations.** On the 180 paired L0 cases, N1 gene recovery is not
   lower than S1b (McNemar, one-sided, α = 0.05) at ≤ 50 % of S1's linear programs, and N1's harm on L0 ∪ L1
   decisions is lower than S1b's by ≥ 0.25 (bootstrap 95 % interval excluding 0).
4. **Pathway level.** Reaction-level recovery ≥ 30 % of tier-A hidden cases with the fill gate on, and ≥ 10
   pathway calls that contain a hidden case's reaction and gene.
5. **Controls.** R recovers ≤ 2 cases per seed; N3 recovers fewer than N1 in ≥ 4 of 5 seeds; the masking audit
   passes for every case.
6. **Genome-level tiers (B + C pooled, own-CLEAN tables, 372 labels).** Secondary verdict: N1 reaches precision
   ≥ 0.75 on ≥ 20 labelled decisions and gene-recovers ≥ 10 hidden cases with S0 at 0; tier A own-CLEAN is reported
   next to it as the same-scale reference. Tier B alone (79 labels) is reported but carries no verdict.

A result that meets 1 and 2 but not 3 is reported as "useful at the measured level but no cheaper than the pilot";
a result meeting none is negative and closes row 12 as "not met, measured". The primary verdict is tier A with the
six-predictor tables; the genome-level verdict is criterion 6.

## 9. Power and expected effect (from the exploration; indicative)

* Tier A hidden cases ≈ 600, of which the hidden protein is rank 1 for its enzyme in about 38 % (≈ 230), the
  reachable set for any top-candidate rule. The pilot's union recovered 21 of its 71 rank-1 cases (30 %); if the
  differential sole-source tests re-add the catabolic reactions they currently bury in 150-reaction fills, 35–50 %
  of rank-1 cases is plausible, i.e. 80–115 gene recoveries before the abstention rule, 40–70 after it. Criterion 2
  (≥ 25) is deliberately modest.
* Precision: the exploration's small-fill calls on labelled proteins had 0.40 (99 / 249 calls, 14 distinct
  decisions) and small-fill mutual-best calls 0.38 per call and 5 / 10 distinct; the margin and free-protein rules
  are expected to remove most of the 30 wrong-gene calls (median margin −0.31 there). With 615 new labels, about
  60–150 labelled decisions per arm are expected in the unmasked runs, giving a precision interval of about ± 0.1.
* Paired comparison (criterion 3): with ≈ 23 S1b recoveries and 40+ N1 recoveries on the 180 L0 cases, the
  one-sided McNemar at α = 0.05 has power > 0.9 for 15 discordant pairs in N1's favour.
* Tier B: 79 labels, about 15–25 expected rank-1 cases; descriptive only.

## 10. Compute (Shaheen, account k10018)

* Masked cases: 615 (tier A) + 180 (A0 re-run) + ≈ 60 usable tier-B cases ≈ 855 cases. Pilot cost was 0.7 CPU h per
  case for the union; the differential T1 fills are smaller but each fill is solved twice: budget 1.0 CPU h per
  case per arm. Arms per case: N1, N2, N3 × 5, three ablations, S1b, R × 5 = 16 solved arms, but N3/R/ablations
  reuse N1's baseline and most LPs; budget 6 CPU h per case → ≈ 5,100 CPU h ≈ 27 exclusive node-hours
  (192 cores). S1 is reused from the pilot for the 180 L0 cases and re-run for tier A at full budget (615 × 0.7 =
  430 CPU h, 2.5 node-hours).
* Unmasked deployment runs: 30 genomes × all arms ≈ 2 node-hours.
* Fan-out: one node per genome per arm group; 30 genomes in parallel; wall-clock per wave 1–4 h; whole run under
  one day. Total ≈ 35 node-hours ≈ 6,700 physical core-hours (13,400 as Slurm counts threads), well inside the
  10 M core-hour allocation. Exclude nid00022, nid00497, nid01033; extract packages once per node; check `kuq`
  (inode quota) before the fan-out.
* Own-CLEAN tables: ESM-1b embeddings and CLEAN inference for 47 proteomes (22 + 8 + 17, about 200,000 proteins)
  on orix (`pi-hohndor` or `freecycle`), well under 5 GPU-hours, using the existing CLEAN checkpoint and the
  embedding workers of the Empty Quarter GPU set. The Price raw matrices of the other five predictors are a
  collaborator's snapshot and cannot be regenerated, which is why tiers B and C use own-CLEAN only.
* Tier C: 17 gapsmith models (about 80 s each, one node) from Fitness Browser proteomes, plus 293 masked cases
  (≈ 300 CPU h, 2 node-hours).

## 11. Implementation and provenance

Code under `benchmark/v3/offgrowth/` (extend, do not fork): `arms.py` (baseline minimal fill, differential T1,
two-solver stability, fill-size recording), `score.py` (reaction gate, assignment rule with abstention, pathway
calls, distinct-decision scoring, volume-matched sequence-only comparator), a new `fitness_browser_labels.py`
(L1 loader with receipt), `panel.py` (measured-growth panel from the follow-up tables), `report.py` (tables of
§7), unit tests for every rule on synthetic models. Shaheen scripts under `benchmark/v3/shaheen/offgrowth2-*.sbatch`.
Written by agy-auto (Gemini) from specs in `~/tmp/gspa-v3-agy/offgrowth2/`, reviewed before submission. Outputs to
`/scratch/project/k10018/gspa-v3/out/offgrowth-v2/`; results and `SHA256SUMS` to
`benchmark/v3/results/offgrowth-v2/`; report `docs/benchmarks/v3-offgrowth-improvement.md`.

Freeze procedure: (1) finalise §4 constants and §8; (2) `sha256sum` of this file into the run receipt; (3) run the
L0 fixed-precision check; (4) only then score tier A/B. Any later edit is a dated entry in §13.

## 12. Risks

* **The differential fill may not expose the catabolic route.** If M_base already contains alternative routes from
  the universal set, the compound-specific fill may be empty (growth without any new reaction) or still large. The
  fill-size distribution after step 1 is the first thing to inspect; a pilot on two genomes (one *Pseudomonas*,
  *S. meliloti*) runs first and is reported in §13 before the fan-out.
* **Label quality of L1.** Fitness Browser reannotations are curated from phenotypes plus homology, some with
  hedged comments ("possible", "probably"); rows with such comments are kept but reported as a separate stratum.
  Many L1 enzymes may already be in Swiss-Prot and in the predictors' training data, which inflates both the arm
  and the sequence-only comparator equally; the "additional functions" metric is defined relative to Q for this
  reason.
* **Repeated families.** Seven closely related *Pseudomonas* genomes hold more than half of L1; pooled counts
  overweight them, hence the split reporting in §5.
* **Gene-assignment ceiling.** Recovery cannot exceed the rank-1 share (≈ 38 %); isozyme pairs (1.3.8.1) stay out
  of reach unless both paralogs pass the co-naming rule.
* **Reaction-level precision is not labelled.** No negative reaction labels exist; reaction-level calls are reported
  with recovery and volume only, and the Biolog negatives that would allow a test are not available.
* **Tier B depends on external predictor runs** and on labels that trained another component; it is confirmatory
  only.
* **Budget matching.** Equal LP counts do not equalise information (one differential test is cheaper and more
  specific than one pilot test); cost is reported alongside and the S1 full arm is kept for the unmatched view.

## 13. Deviations

Dated, append-only. None of the entries below was written after reading a tier A, B or C score.

1. **2026-11-02, vitamins in sole-source media (implementation).** The pilot's element-replacement function strips
   every compound containing the tested element, which would remove thiamin and riboflavin from C, N and S test
   media. The implementation restores the two vitamins at their minimal-medium bounds after replacement, as §4
   step 1 prescribes; the unit test that had encoded the stripping was corrected. No run had been scored.
2. **2026-11-02, two-genome pilot of the differential fill (plan §12, risk 1).** On *S. meliloti* and
   *Pseudomonas* N2C3 (3 masked cases each plus the unmasked run; six-predictor tables) the baseline minimal fill
   adds 141 and 73 reactions, T2/T3 fills have a median of 1–2 reactions, N/P/S differential sole-source fills a
   median of 2–3, but C-source differential fills are bimodal: 13 of 109 and 3 of 91 C tests fill ≤ 3 reactions,
   while 74 and 67 fill ≥ 40. The large C fills share a block of 29 reactions present in ≥ 90 % of them
   (carbamate kinase, homoisocitrate, methylmalonyl-CoA, corrinoid and CO-dehydrogenase steps): a candidate-
   supported detour that the 1 : 100 weighting prefers over a single unsupported reaction. Targeted re-runs of
   first-round recoveries show the intended effect where the catabolic step has a candidate: mannose isomerase
   returns in fills of 1 (*Phaeobacter*) and 3 (*Acidovorax*) instead of 169 and 127; the *Azospira*
   methylmalonyl-CoA epimerase recovery of the first round sits inside the detour block (fills 24–49) and the
   *Paraburkholderia* mannose-isomerase fills of size 1 are unstable (the two solvers pick different duplicate
   reactions). Decision: the frozen procedure (gate on the raw stable fill, F = 3) is kept as primary. Cost per
   masked case is about 70–85 s wall on one node (about 7,000 CPU s, 8,000–9,000 LPs), within budget.
3. **2026-11-02, one additional reported scoring variant (not primary).** `N1_blocksub`: identical to N1 except that
   fill size is measured after removing the run's common block, defined label-free per run as the reactions
   present in at least half of that run's T1 fills with ≥ 20 reactions. It is reported next to the ablations and
   carries no pass criterion; it exists to show whether the detour block hides specific calls.
4. **2026-11-02, sequence-only comparator on seen labels (pre-registered step, §8).** On L0 (180 labelled proteins
   with a candidate) the six-predictor best-EC rule has exact top-1 precision 0.233 and reaches 0.75 at no
   threshold with ≥ 20 decisions, so `Q_fixed` abstains on every tier; the volume-matched Q is unaffected. Recorded in
   `benchmark/v3/results/offgrowth-improvement-v1/q_fixed_threshold_L0_six.txt` before any held-out score.
5. **2026-11-02, label files.** Tier A label rows: 643 (615 proteins) after removing Price-148 and Swiss-Prot
   proteins and partial ECs; tier C: 318 rows (293 proteins); tier B: 79 rows in 7 genomes (*M. magneticum* has no
   Swiss-Prot experimental EC). Hedged rows (66 tier A, 24 tier C) are kept and reported as a stratum.
6. **2026-11-03, implementation fixes found while running (no rule, constant or metric changed).** (a) The
   evaluation built every genome's scoring context in the parent process and then forked workers, which exhausted
   memory on 22 genomes; it now scores one genome per worker and builds each genome's permuted candidate maps once
   instead of once per run. A test checks that the two paths give byte-identical tables. (b) The measured-growth
   panel arm silently ran as a copy of the transporter-panel arm in the first complete scoring pass, because the
   organism map was read under the wrong column name; after the fix every tier A run uses a real measured panel and
   the whole tier was scored again from scratch. No held-out number had been read when (a) was fixed; (b) was found
   by reading the fallback flag in the cost table, and the corrected numbers are the ones reported.
7. **2026-11-03, tier A reported before tiers B and C.** Tier A needs only the existing six-predictor candidate
   tables; tiers B and C wait on CLEAN inference for 47 proteomes. The tier A verdict is reported on its own and
   criterion 6 is left open.
8. **2026-11-03, tiers A-clean, B and C as run.** CLEAN was run by us on 47 proteomes on orix through the
   repository's sidecar with GMM confidence and no score floor; candidate tables keep its calls with confidence
   above zero and a complete four-level enzyme number, which is the plan's §3 definition and drops about a quarter
   of the raw calls. Two proteins were scored from sanitised sequences (one ambiguous residue, thirteen internal
   stop codons replaced by X) because the language model cannot read those characters; both are recorded in the
   inference receipt and their rows are kept unchanged. Tiers B and C have no separate seen label set, so their
   threshold and test labels are the same file: the arms fit nothing, but the sequence-only comparator's fixed
   threshold is fitted on the labels it is scored on and is therefore reported as descriptive only. Tier B's eight
   bacteria are not Fitness Browser organisms, so its measured-panel arm falls back to the transporter panel in
   every run. One tier B genome carries no label row and contributes only its unmasked run.

