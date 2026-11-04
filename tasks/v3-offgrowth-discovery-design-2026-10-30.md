# Finding enzymes that growth on one medium cannot reveal: pre-registered design (2026-10-30)

Status: design, frozen before any run. Changes after the first result are recorded in §11 with a reason.
Requested by Robert (2026-10-30): "Do it. I suggest to not use Qwen yet, but do something deterministic, rule-based instead."

## 1. Question

Gap-filling a draft metabolic model towards biomass on one medium adds only the reactions that growth on that medium
needs. Most enzymes (alternative substrates, salvage, detoxification, secondary metabolism) are invisible to it.
We ask:

1. Do targeted gap-filling tests — sole-source growth, unblocking of evidenced reactions, and production of dead-end
   metabolites — recover experimentally assigned enzymes that biomass gap-filling misses?
2. How many simulations does each test need per genome?
3. Does a fixed, rule-based selection of tests keep most of the recovery at a fraction of the simulations?

No language model is used anywhere in this experiment.

## 2. Truth set

**Price-148**: the 148 reannotated enzymes in 22 bacteria of Price et al. 2018 (*Nature* 557:503, mutant phenotypes from
RB-TnSeq), as prepared in `benchmark/v3/prepare_price.py` (manifest `benchmark/manifests/suite6_price149_manifest.tsv`;
148 of 149 rows after the documented exclusion; complete four-level ECs only). Proteomes: the staged
`price-raw-sources/genomes/<assembly>/protein.faa` (hashes in `price-full-matrix-v2/layout/report.json`).

**Strata** (used only for reporting, never as an input to any arm): a truth reaction is *on the growth path* if deleting it
from the unmasked, gap-filled model reduces biomass on the default gapsmith medium below 1 % of the wild type; otherwise
*off the growth path*. The primary analysis is the off-growth stratum.

Unconfirmed labels are never negatives (§6).

## 3. Masking

For each truth protein P with EC E (one masked case per (P, E)):

- remove every GSPA claim of P for E and for ECs that map to the same reactions;
- in the gapsmith draft model, remove P from the gene–reaction rules of every reaction mapped from E; a reaction left with
  no gene loses its genomic evidence (it stays in the universal database, so gap-filling can add it back);
- the masking audit checks that no arm input (model, claims, evidence tables, rule features) still links P to E.

Masking reuses the existing masked-assignment builder and its audit where possible (`MaskedAssignmentBuilder`,
`audit_wp5_masking.py`), extended to the new inputs.

## 4. Models and engine

Draft and gap-filled models from the October Empty Quarter recipe (gapsmith 0.1.0 binary 7b7ff538…, DIAMOND 2.1.24,
reference `data_merged` 25bad8d2…), built for the 22 genomes by the Shaheen rebuild in benchmark mode
(`benchmark/emptyquarter_mags/shaheen/`, `--mode benchmark`).

**Candidate genes** ("GSPA candidates" below): for each genome, the union of the six predictors' release-candidate
integration outputs (collection `price-production-rc1`, arm `stage1` = production integration with default priors, all
retained candidates), keyed by (protein, EC) with score = the maximum over the six predictors. Retention floor: the
engine's own retention (every retained candidate counts). Masking (§3) removes the truth (P, E) pairs from this table.

Targeted gap-filling is a weighted linear program on the gapsmith draft plus the gapsmith universal SEED reaction set:
minimise the weighted sum of added reaction fluxes subject to a target flux ≥ ε (ε = 0.01 of the target's maximum in the
universal-completed model, fixed here). Reaction weights: 1 for a reaction with a GSPA candidate gene above the retention
floor, 100 otherwise (gapsmith's default weighting scheme in spirit; the values are fixed here, not tuned). Solver: HiGHS
through COBRApy (the Python 3.10 test environment), with GLPK replay on a fixed 5 % sample of tests as an independent check
(agreement criterion: same added-reaction set, objective within 1e-6 relative).

## 5. Arms (all deterministic)

| Arm | Tests | Purpose |
|---|---|---|
| A0 biomass | biomass gap-fill on the default gapsmith medium (the production model) | baseline: what we do now |
| A1 sole source | for every compound of the frozen source panel (§5.1): growth with it as the only C (or N, P, S) source; gap-fill where growth fails | substrates |
| A2 unblock | every reaction with gene evidence in the masked draft that cannot carry flux under the union medium (flux variability) is a target; gap-fill until it can | disconnected pathways |
| A3 dead ends | every metabolite produced but not consumed (or consumed but not produced) by gene-evidenced reactions gets a sink (or source) target; gap-fill | end products, salvage |
| A4 rules@B | the A1–A3 tests ordered by the fixed rules of §5.2, cut at budget B ∈ {10 %, 25 %} of that genome's A1–A3 test count | fewer simulations |

### 5.1 Source panel

Compounds of the Biolog PM1–PM4 plates (carbon, nitrogen, phosphorus, sulfur sources) that map to a ModelSEED compound with
an exchange and transport route in the gapsmith universal database. The mapping table is built and frozen before any run
(file and sha256 recorded in the run receipt). Unmapped compounds are listed, not dropped silently.

### 5.2 Rules for A4 (fixed order; ties by test identifier)

1. Tests of a source compound that a transporter annotation of the genome (GSPA claim or gapsmith transporter table) names.
2. Tests whose target lies in a KEGG module that is ≥ 50 % evidenced in the masked genome (the module's missing steps are
   what the test could fill).
3. Unblock and dead-end targets with at least one GSPA candidate gene (above the retention floor) for a reaction on the
   cheapest fill path in the universal-completed model (a structural lookup, no LP).
4. All remaining A1–A3 tests.

Within a rule, tests are ranked by the number of distinct candidate genes supporting them (more first). The rules use
only label-free, masked inputs.

## 6. Calls and scoring

- **Call.** A reaction added by any test of an arm is assigned to the GSPA candidate gene with the highest score for an EC
  mapped to that reaction (retained weak candidates included). No candidate gene: the reaction is called with no gene
  (counted as "reaction only").
- **Recovered.** Masked case (P, E) is recovered by an arm if the arm calls P for a reaction mapped from E.
- **Harm.** A call of a reaction to a protein whose experimentally supported function (Price-148 truth, or Swiss-Prot
  experimental EC for that protein) maps to a different reaction. Calls on proteins without experimental support are
  neither correct nor harmful; they are reported as unconfirmed.
- **Cost.** Number of linear programs and gap-fills, CPU seconds and peak memory per genome and arm.

Metrics per arm: recovery (primary: off-growth stratum), harm rate among calls on experimentally supported proteins,
unconfirmed call count, cost. Genome-block bootstrap (2,000 resamples, seed 20261030) for intervals. No threshold is fitted,
so there is no calibration split; all 22 genomes are evaluated once.

## 7. Pre-registered comparisons

- H1: A1 ∪ A2 ∪ A3 recovers more off-growth cases than A0 (paired by case, one-sided exact McNemar).
- H2: each of A1, A2, A3 separately versus A0 (same test), Holm-corrected across the three.
- H3: A4@25 % keeps ≥ 80 % of the A1 ∪ A2 ∪ A3 off-growth recoveries (descriptive with bootstrap interval).
- Harm: any arm whose harm rate exceeds A0's by more than 10 percentage points is reported as harmful regardless of recovery.

A negative result is reported as such. The sample (148 cases, fewer off-growth) may be too small for H2; that is stated,
not hidden.

## 8. Controls

- Shuffled candidate genes: rerun A1–A3 call assignment with GSPA candidate genes permuted within each genome (fixed seeds 1–5);
  recovery must drop to near zero, else the call rule leaks.
- Unmasked run: every arm on the unmasked model, to separate "the test can reach this reaction" from "masking removed it".
- Solver replay agreement (§4).

## 9. Execution

Pilot on the 22 genomes on Shaheen (one node per genome), recording cost per test type. The pilot results answer question 2
and size any later run on the Empty Quarter (17,055) and AllTheBacteria genomes (species representatives first).

## 10. Deliverables

Scripts under `benchmark/v3/offgrowth/` (panel builder, masking, arm runner, scorer, report) with plain unittest tests;
Shaheen job scripts; a report `docs/benchmarks/v3-offgrowth-discovery.md` (numbers in tables, plain-word narrative, negative
results stated).

## 11. Deviations

All recorded before the pilot results were read, except D6, which concerns reporting only.

1. **Call attribution (§3, §6).** Read literally, §3 removes the hidden protein's claims from the candidate table and §6
   assigns added reactions with that same table, so no arm could ever recover a case. Arms are solved on the fully
   masked inputs (model, candidate weights, transporter table, rule features; masking audit passes), and added
   reactions are then assigned with the unmasked candidate table. The literal version is also reported
   ("masked-candidate attribution", zero by construction).
2. **Solver (§4).** COBRApy on Shaheen has no HiGHS interface (it needs osqp, not installed). Primary solver GLPK; the
   replay solves the same linear programs, extracted as sparse matrices from the configured GLPK model, with
   scipy.optimize.linprog (HiGHS), on the fixed 5 % sample. The replay agreement criterion is unchanged.
3. **Candidate table for dense predictors (§4).** For genomes whose union exceeds 250,000 (protein, enzyme) pairs the
   table keeps, per enzyme, the 64 top-scoring proteins and all proteins tied at the top score, plus every protein's
   best enzyme. Attribution uses only the top candidate and reaction weights only need one candidate, so neither is
   changed. An earlier version also kept every claim of the truth proteins; that label-dependent input was removed before
   any run.
4. **Harm labels (§6).** Swiss-Prot experimental ECs for the 22 genomes did not exist and were built
   (`benchmark/v3/offgrowth/swissprot_labels.py`: reviewed entries, ECO:0000269 and descendants, exact-sequence match;
   35 labels on 32 proteins). Genome runs scored harm with Price-148; the primary report re-derives every call's class
   with Price-148 plus these labels; both reports are kept.
5. **Source panel (§5.1).** Built from Biolog's plate maps (Rev E, sha256 9cef31cd…); 35 manual name mappings with a
   reason each, one corrected by the supervisor (β-glycerol phosphate = glycerol 2-phosphate, cpd01908). 293 of 379
   substrate wells usable; the rest listed.
6. **Reporting.** Recovery counts for "reaction re-added" are reported next to protein recovery, as a diagnostic.
7. **(2026-11-02) Full benchmark panel, candidate source.** For the full-panel run (addendum
   `tasks/v3-offgrowth-full-panel-2026-11-02.md`, frozen sha256 a7be7111… before any scored check) the six-predictor
   union of §4 exists only for the 22 Price genomes and this design names no fallback; all 145 genomes therefore use
   the genus-fenced Swiss-Prot DIAMOND baseline adapter (`homology_baseline.py`, training parameters) with a repaired
   fence (taxid plus Swiss-Prot lineage instead of the first word of the organism name) and a leakage check. Hidden
   cases add Swiss-Prot experimental ECs of all 145 genomes (D4 rule). Everything else in §3–§8 and D1–D6 is unchanged.
