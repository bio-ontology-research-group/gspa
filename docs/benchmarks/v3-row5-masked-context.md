# Reaction context on masked enzyme loci with known truth (acceptance row 5)

Status: complete (2026-11-03). Plan frozen before any arm ran: `tasks/v3-row5-context-improvement-2026-11-02.md`
(SHA-256 `ed8a65ea…` in `benchmark/v3/results/row5-masked-context-v1/PLAN_SHA256`); every later change is a dated
entry under its "Deviations". Results folder: `benchmark/v3/results/row5-masked-context-v1/` (`SHA256SUMS`,
`LARGE_FILES.tsv`; evaluation tables under `evaluation-{development,heldout,all}/`, the pilot's learned-stage
arms under `evaluation-pilot-learned/`). Raw cell outputs stay on Shaheen under
`/scratch/project/k10018/gspa-v3/out/row5-masked-context-v1/` (`cells-main/`, `cells/` = pilot) and the retained
scores on ws under `/mnt/data2/leechuck/gspa-v3-data/row5-masked-context-v1/`.

## Question

Can reaction context, the gap-filling prior and the universal, pathway-only and organism-specific context modes,
raise the right enzyme number on the right gene when that enzyme's reaction is a gap of the genome's own metabolic
model? The earlier row-5 test on the 22 Price genomes could not answer this: only 8 of 148 Price enzymes were gaps of
their model, the modes act only through the learned stage, and the learned stage declines every candidate once gap
observations are present (`v3-price-calibration-context.md`, `v3-price-contract-learned.md`).

## Design in brief

* Truth: experimentally supported (protein, EC) pairs in the 40 benchmark genomes that have gapsmith models
  (22 Price genomes: Price-148 and Swiss-Prot experimental ECs; 18 cohort genomes: Swiss-Prot experimental ECs).
* Masking: for one labelled enzyme at a time, every `good_blast` row of the gapsmith pathway-search table that names
  that protein is set to `no_blast`, so its reaction becomes a gap of a predicted pathway whose other steps are still
  present. Claims are not masked. One cell per masked locus and claim set; a sham cell masks an unlabelled protein of
  another reaction of the same pathway; an unmasked cell per genome and claim set.
* Claim sets: contract-bearing DIAMOND Swiss-Prot claims (genus-fenced; the Price claims of the other session,
  cohort claims from the same driver) and the six Price enzyme predictors (the dense ones filtered, label-free, to
  each protein's top-10 ECs and each EC's top-64 proteins; the filtered stage 0 reproduces the release stage 0 on the
  labelled loci, see the audit table).
* Arms: stage 0 (no priors); default priors with the masked gaps on the release engine (`gaps`); the same on the
  3.1.0-dev engine with the EC-namespace repair of the gap prior (`gapsfix`); the repair plus the context-weighted gap
  prior in universal, pathway-only and organism-specific mode (`ctx_*`, window 5, 10 and 20 kb); the learned-stage
  mode arms on the pilot loci only (plan A1: with gap observations present the contract-3 models decline every
  candidate, so these arms are identical to `gaps` by construction).
* Outcome per locus and claim set: call at 0.5, protein top-1, and "gene for the gap" (the protein is the genome's
  top-scoring protein for that EC); a locus is improved when an indicator is gained and none lost, harmed when one is
  lost (mixed counts as harmed). Harm elsewhere counts other labelled pairs of the genome that lose a call or top-1.
* Decision rule, frozen: one-sided exact sign test on discordant loci (p ≤ 0.05), improved ≥ harmed in every label
  set, claim set with ≥ 3 discordant loci and stratum, net gain ≥ 3 loci and harm elsewhere not above the gain;
  outcomes "systematic improvement", "no detectable effect", "inconsistent", "harm". The 58 development loci
  (cohort train genomes) chose the window only (all three windows tied; 10 kb by the tie rule); the 61 held-out loci
  (cohort calibration and test, all Price) decide.

## What the pilot found before any contrast was evaluated

The first context arm of the pilot reported that the context support was never consulted. The cause is a latent
defect of the release engine: integrated EC candidates are keyed `EC:x.x.x.x` while gap records carry bare
`x.x.x.x`, so the gap-filling prior's targeted rescue never matches an EC candidate and only GO terms mapped from
gap ECs are reached. This explains the earlier row-5 result (the gap prior changed GO rows but no Price EC call). A
repair flag was added to the 3.1.0-dev engine (`--gap-prior-ec-namespace-fix`, off by default; the default stays
byte-identical to 3.0.1), two arms were added (`gapsfix`, `sham_gapsfix`), every context arm runs with the repair on,
and the contrasts were amended before any outcome was read (plan deviation of 2026-11-02). The repair is a method
change on the same footing as the context weighting and is reported as such.

## Loci, cells and audits

| Quantity | Value |
|---|---|
| Labelled (protein, EC) pairs on the 40 genomes | 521 |
| Maskable loci with pathway context | 119 (cohort Swiss-Prot 87, Price-148 25, Price Swiss-Prot 7; 28 genomes) |
| of which gaps only through masking / gaps already | 107 / 12 |
| Development / held-out loci | 58 / 61 |
| Cells (locus × claim set; cohort loci: DIAMOND; Price loci: DIAMOND and six predictors) | 311 (held-out 253) |
| Arms per cell (main run) | 19 |
| Cells complete, candidate universe preserved, masked gap derived, `gaps_ctxjar` byte-identical to `gaps` | 311 / 311 / 311 / 311 |
| Already-gap loci whose unmasked and masked gap records agree | 58 / 58 (held-out cells) |
| Cohort leakage (291 labelled proteins vs the DIAMOND reference) | 0 exact, 0 near-identical, 0 same-organism entries |
| Compute | staging 12 min; pilot 23 cells, one node; main run 12 nodes × 1 h 10 min; evaluation 3 jobs |

## Results on the held-out loci (61 loci, 253 cells)

Improved / harmed count loci × claim set; "harm elsewhere" counts other labelled pairs of the same genome; the
interval is the genome-block bootstrap of the net improved count (2,000 resamples).

| Contrast | Improved | Harmed | Discordant | p (one-sided sign test) | Harm elsewhere | Unconfirmed calls created | Net, 95 % interval | Verdict |
|---|---:|---:|---:|---:|---:|---:|---|---|
| C1 release gap prior: `stage0 → gaps` | 0 | 0 | 0 | 1 | 0 | 0 | [0, 0] | no detectable effect |
| C1b repair: `gaps → gapsfix` | 47 | 2 | 49 | 2.2e-12 | 0 | 20 | [30, 61] | **systematic improvement** |
| C5 universal: `gapsfix → ctx_universal` (10 kb) | 1 | 6 | 7 | 0.99 (harm 0.06) | 14 | 0 | [−10, 0] | no detectable effect |
| C5 pathway-only: `gapsfix → ctx_pathway` | 1 | 2 | 3 | 0.88 (harm 0.5) | 2 | 0 | [−4, 2] | no detectable effect |
| C5 organism-specific: `gapsfix → ctx_organism` | 1 | 5 | 6 | 0.98 (harm 0.11) | 13 | 0 | [−8, 0] | no detectable effect |
| C5 pathway-only vs universal | 4 | 0 | 4 | 0.0625 | 0 | 0 | [1, 8] | no detectable effect |
| C5 organism-specific vs universal | 1 | 0 | 1 | 0.5 | 0 | 0 | [0, 3] | no detectable effect |
| C5 organism-specific vs pathway-only | 0 | 3 | 3 | 1 (harm 0.125) | 11 | 0 | [−6, 0] | no detectable effect |
| S1 sham mask: `sham_gaps → sham_gapsfix` | 11 | 0 | 11 | 4.9e-4 | 0 | 20 | [3, 21] | see below |
| S1 sham mask: `sham_gapsfix → sham_ctx_universal` | 0 | 2 | 2 | 1 | 13 | 0 | [−5, 0] | no detectable effect |
| U1 unmasked versus masked `gaps` | 0 | 0 | 0 | 1 | 0 | 0 | [0, 0] | no detectable effect |

The repair (C1b) by claim set and label set, held-out loci:

| Slice | Cells | Improved | Harmed | p | Verdict |
|---|---:|---:|---:|---:|---|
| DIAMOND Swiss-Prot claims | 61 | 22 | 0 | 2.4e-7 | systematic improvement |
| CLEAN (filtered) | 32 | 6 | 0 | 0.016 | systematic improvement |
| GraphEC (filtered) | 32 | 8 | 0 | 0.0039 | systematic improvement |
| DeepProZyme (filtered) | 32 | 4 | 0 | 0.0625 | no detectable effect |
| MAPred | 32 | 5 | 2 | 0.23 | no detectable effect |
| EnzBERT | 32 | 2 | 0 | 0.25 | no detectable effect |
| TopEC | 32 | 0 | 0 | 1 | no detectable effect |
| cohort Swiss-Prot / Price-148 / Price Swiss-Prot | 29 cells / 175 / 49 | 16 / 25 / 6 | 0 / 2 / 0 | | consistent |
| archaea / bacteria | 6 / 247 | 3 / 44 | 0 / 2 | | consistent |
| gap already / gap through masking | 58 / 195 | 11 / 36 | 0 / 2 | | consistent |
| pathway prediction marginal / robust | 122 / 131 | 26 / 21 | 1 / 1 | | consistent |

Genome-wide cost of the repair on the held-out cells: 2,889 more calls at or above 0.5 and 10,792 changed scores
over 253 cells (about 11 new calls per cell), against 47 corrected loci and 20 unconfirmed calls on labelled
proteins (`evaluation-heldout/tables.md`, section 2).

The context modes on the development loci (58, DIAMOND claims; the window choice set) and on all 119 loci
(sensitivity, not primary):

| Contrast | Development: improved / harmed, p (harm) | All loci: improved / harmed, p (harm) |
|---|---|---|
| C1b repair | 22 / 0, improvement p 2.4e-7 | 69 / 2, improvement p 1.1e-18 |
| C5 universal | 0 / 7, 0.0078 (harm) | 1 / 13, 9.2e-4 (harm) |
| C5 pathway-only | 0 / 6, 0.016 (harm) | 1 / 8, 0.020 (harm) |
| C5 organism-specific | 0 / 6, 0.016 (harm) | 1 / 11, 0.0032 (harm) |
| C5 pathway-only vs universal | 1 / 0 | 5 / 0, improvement p 0.031 |
| C5 organism-specific vs universal | 1 / 0 | 2 / 0 |
| C5 organism-specific vs pathway-only | 0 / 0 | 0 / 3 |

Learned-stage arms (pilot loci only, 5 loci × claim sets): `evaluation-pilot-learned/tables.md`. With gap
observations present the learned stage scored no candidate in any mode, so `gaps_learned`, `universal`, `pathway`
and `organism` equal `gaps` row for row (0 differing scores in every pilot cell), as the other session's step F
had shown on the unmasked Price gaps.

## What this shows

The release engine's gap prior does nothing for enzyme numbers: on 253 held-out cells with a real, derived gap at the
masked locus, `stage0` and `gaps` give the same call, the same top-1 and the same top gene everywhere. This is the
EC-namespace defect, not a lack of signal. Once gap records and EC candidates are compared in one namespace, the same
prior corrects 47 of the 253 held-out cells and harms 2, with no harm to other labelled pairs of the same genome, and
the gain is consistent across label sets, claim sets, domains, gap kinds and pathway strata. By the frozen rule this
is a systematic, consistent improvement, on held-out loci, from reaction context in the plain sense: the gap of a
predicted pathway lifts the right enzyme number on the right gene above the call threshold. The sham control
separates the mechanism from the masking: on loci that become gaps only through masking, masking an unrelated
protein of the same pathway corrects nothing (0 of 195 cells), while the true mask corrects 36; on loci that are gaps
without masking, the repaired prior acts under either mask (11 of 58), as it should.

The context weighting of that prior, which is where the three modes act, does not help. Where it changes an outcome
it mostly removes a call: by construction a candidate with no anchored pathway neighbour within the window receives
half of the rescue boost and falls back below 0.5. On the held-out loci this is "no detectable effect" for all three
modes and for every mode-versus-mode contrast; on the development loci and on all 119 loci the direction is harm at
the rule's significance level. Pathway-only context is the least harmful mode and differs from universal context
only in its favour (4 to 0 held-out, 5 to 0 over all loci), but not enough on the primary loci to count as an
improvement. The modes differ among themselves by one to a few loci because they differ in how often a gap has any
anchored context at all; none of them finds candidates that the repaired prior misses.

## Verdict against the pre-registered criteria

| Contrast (held-out loci) | Verdict |
|---|---|
| C1 release gap prior | no detectable effect (structurally null for EC candidates) |
| C1b EC-namespace repair of the gap prior | systematic improvement |
| C5 universal, pathway-only, organism-specific context on top of the repair | no detectable effect, direction harm (harm on the development and all-loci scopes) |
| C5 mode versus mode | no detectable effect; pathway-only ≥ universal in every scope |
| Learned-stage mode arms (pilot) | structurally null (0 scored candidates with gap observations) |

Row 5 reading (plan 3.7): the pathway-only versus universal versus organism-specific ablation is now measured on
61 held-out enzyme loci (253 cells) that carry actual, derived gaps with pathway context, in the deterministic
gap prior, with unmasked and sham controls; the measured answer is that none of the three modes improves the
assignment and that the mode choice changes a handful of loci. The one systematic, consistent improvement found is
the repair of the gap prior itself, a 3.1.0-dev method change that is off by default and recommended for the next
release after a Price and cohort regression run of the full pipeline with the flag on.

## Limits

* The learned-stage modes could not be tested with gap observations: the contract-3 models decline every candidate
  whose availability profile includes gap features. A model trained with gap features is needed for that arm.
* Masking keeps gapsmith's pathway prediction flag; 46 of the 119 loci are "marginal" (the pathway's completeness
  would drop below the threshold without the masked reaction). The repair's gain is present in both strata.
* Dense predictors were filtered to top-10 ECs per protein and top-64 proteins per EC; outcomes on their unlabelled
  proteins are not evaluated. TopEC emits one EC per protein and cannot gain.
* The sham control is uninformative for the 12 already-gap loci (the true gap persists under any mask); the control
  rests on the 107 masked-gap loci.
* Cost is reported, not judged: the repair adds about 11 calls at or above 0.5 per genome cell, none of them on a
  labelled protein with a different label, and 20 unconfirmed calls on labelled proteins over 253 cells.
* One gspa-core test (`ExtractorHardeningSpec`, a 20-second wall-clock bound unrelated to this change) failed during
  the rebuild under a load average above 80 on ws; its rerun is recorded in the plan's deviations.
