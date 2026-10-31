# Finding enzymes that growth on one medium cannot reveal: pilot on the 22 Price bacteria

Status: pilot complete (2026-10-30). Pre-registered design: `tasks/v3-offgrowth-discovery-design-2026-10-30.md`
(deviations in its §11). Code: `benchmark/v3/offgrowth/` (20 tests), Shaheen jobs `benchmark/v3/shaheen/offgrowth-*.sbatch`.
Results: `benchmark/v3/results/offgrowth-pilot-v1/` (`report-price148-swissprot/` is the primary report; `report-price148/`
counts harm against the Price-148 labels only). No language model is used anywhere.

## Question

Gap-filling a metabolic model to grow on one medium only adds what that growth needs. We asked whether targeted
gap-filling tests recover experimentally assigned enzymes that biomass gap-filling misses, how many simulations that
takes, and whether fixed rules can choose a small subset of the tests.

## Set-up in brief

* Truth: the 148 Price et al. (2018) reannotated enzymes in 22 bacteria. For each case the protein's evidence for its
  enzyme is hidden from the model, the candidate table and every rule input (masking audit passes for all cases).
* Models: gapsmith draft and gap-filled models of the 22 genomes (October Empty Quarter recipe, rebuilt on Shaheen).
* Candidate genes: union of the six enzyme predictors' outputs from the release-candidate Price collection, best score
  per protein and enzyme.
* Tests (arms): biomass gap-filling on the default medium (A0, current practice); growth on each mapped Biolog PM1–PM4
  compound as sole C, N, P or S source (A1; 293 of 379 substrate wells usable, the rest listed with the reason in the
  report); making blocked gene-evidenced reactions carry flux (A2); giving dead-end metabolites an outlet (A3); their
  union; and the fixed rule order cut at 10 % and 25 % of each genome's tests (A4).
* A reaction added by a test is assigned to the top-scoring candidate gene for an enzyme mapped to it. A case is
  recovered when the hidden protein is assigned its true reaction. Harm counts only calls on proteins with experimental
  labels (Price-148 or Swiss-Prot experimental EC, 35 labels on 32 proteins); all other calls are unconfirmed.

## Results

### Recovery (148 cases; 146 off the default growth path)

| Test | Reaction re-added | Hidden protein recovered | Recovery, 95 % interval (genome bootstrap) |
|---|---:|---:|---|
| A0 biomass, default medium | 0 | 0 | 0.000 |
| A1 sole sources | 28 | 9 | 0.061 [0.024, 0.104] |
| A2 unblock reactions | 27 | 14 | 0.095 [0.049, 0.143] |
| A3 dead-end metabolites | 23 | 12 | 0.081 [0.030, 0.134] |
| A1 + A2 + A3 | 51 | 21 | 0.142 [0.086, 0.197] |
| A4 rules, 10 % of tests | 20 | 8 | 0.054 [0.023, 0.091] |
| A4 rules, 25 % of tests | 25 | 11 | 0.074 [0.044, 0.107] |

### Pre-registered comparisons (off-growth cases)

| Comparison | Test only | A0 only | p (exact McNemar, one-sided) | Holm-adjusted | Outcome |
|---|---:|---:|---:|---:|---|
| H1: A1 + A2 + A3 vs A0 | 21 | 0 | 4.8e-07 | — | supported |
| H2: A1 vs A0 | 9 | 0 | 0.0020 | 0.0020 | supported |
| H2: A2 vs A0 | 14 | 0 | 6.1e-05 | 0.00018 | supported |
| H2: A3 vs A0 | 12 | 0 | 0.00024 | 0.00049 | supported |
| H3: A4 at 25 % keeps ≥ 80 % of the union's recovery | 11 of 21 = 0.52 [0.40, 0.67] | | | | not met |

### Calls and harm (Price-148 + Swiss-Prot labels)

| Test | Gene calls | Supported | Harmful | Unconfirmed | Harm rate, 95 % interval |
|---|---:|---:|---:|---:|---|
| A0 | 9,543 | 0 | 22 | 9,521 | 1.00 |
| A1 | 52,693 | 55 | 474 | 52,164 | 0.90 [0.69, 0.99] |
| A2 | 67,994 | 135 | 550 | 67,309 | 0.80 [0.72, 0.96] |
| A3 | 43,945 | 102 | 336 | 43,507 | 0.77 [0.69, 0.94] |
| A1 + A2 + A3 | 105,539 | 160 | 892 | 104,487 | 0.85 [0.78, 0.93] |
| A4 at 10 % | 40,063 | 55 | 387 | 39,621 | 0.88 [0.66, 0.99] |
| A4 at 25 % | 52,889 | 72 | 471 | 52,346 | 0.87 [0.70, 0.97] |

Counts add up over the 148 cases, so the same call on a nearly identical model is counted once per case.

### Controls

| Control | Result |
|---|---|
| Shuffled candidate genes (seeds 1–5), union | 0 of 148 recovered for every seed |
| Literal masked-candidate attribution | 0 recovered in every arm (by construction) |
| Unmasked model and candidates | union recovers 9 cases whose true reaction is missing from the gapsmith model even without masking |
| Solver replay (GLPK vs HiGHS, 8,912 gap-fills) | objective agrees in every replay (largest relative difference 3.7e-10); the same set of added reactions in 4,227 (47 %) |

### Cost (all 148 cases)

| Test | Tests run | Linear programs | CPU hours |
|---|---:|---:|---:|
| A0 | 148 | 740 | 0.3 |
| A1 | 43,364 | 128,712 | 29.5 |
| A2 | 153,181 | 621,769 | 59.3 |
| A3 | 63,093 | 97,958 | 18.4 |
| A1 + A2 + A3 | 259,638 | 848,439 | 107.1 |
| A4 at 10 % | 26,031 | 69,344 | 17.8 |
| A4 at 25 % | 64,945 | 144,290 | 31.8 |

About 1,750 tests (0.7 CPU hours) per case and 4.9 CPU hours per genome for the full union; one 192-core node finished
each genome in 4 to 40 minutes.

### Extension: Swiss-Prot experimental enzymes as additional hidden cases

Approved by Robert on 2026-10-30, after the pilot: the 32 proteins with Swiss-Prot experimental enzyme numbers in 6 of
the 22 genomes (no overlap with Price-148) were hidden and tested in the same way (reports `report-swissprot-cases/`,
`report-combined/`). All 32 lie off the default growth path.

| Test | Swiss-Prot cases: reaction re-added | Swiss-Prot cases: protein recovered | Both sets (178 cases): protein recovered |
|---|---:|---:|---:|
| A0 biomass | 0 | 0 | 0 |
| A1 sole sources | 2 | 1 | 10 |
| A2 unblock | 2 | 2 | 16 |
| A3 dead ends | 1 | 1 | 13 |
| A1 + A2 + A3 | 3 | 2 | 23 |
| A4 at 25 % | 1 | 1 | 12 |

The two recovered Swiss-Prot enzymes are both in *Sinorhizobium meliloti* (UDP-glucose 6-dehydrogenase 1.1.1.22,
glutamate-1-semialdehyde 2,1-aminomutase 5.4.99.2). On the 178 combined cases H1 and all three H2 comparisons hold
(23 vs 0, one-sided exact McNemar p = 1.2e-07; Holm-adjusted H2 p ≤ 0.001) and H3 is again not met (12 of 23 = 0.52).
Classic, well-studied enzymes are recovered at a lower rate than the Price reannotations (2 of 32 vs 21 of 148).

## What this shows

Biomass gap-filling on one medium recovered none of the hidden enzymes, which is expected: removing any one of them
leaves the gapsmith model growing (only 2 of 148 cases lie on the default growth path). Targeted tests recovered 21, in
11 of the 22 genomes, and every pre-registered comparison against biomass gap-filling holds. The recovered enzymes are
mainly in sugar and nucleoside catabolism (rhamnose isomerase and rhamnulose-phosphate aldolase, mannose isomerase,
uridine nucleosidase, uronate dehydrogenases), the kind of function a sole-source or unblocking test can expose. The
shuffled-gene control gives zero, so the gene assignment does not leak the answer.

The effect is small and the calls are mostly not usable as they stand: the right reaction comes back in 51 cases but
the hidden protein wins the gene assignment in only 21, and most calls on proteins with experimental labels are wrong.
The assignment rule (top-scoring candidate for any enzyme mapped to the reaction) is too crude for promiscuous
predictor scores. The fixed rules do not keep enough of the recovery at a quarter of the tests (H3 not met). The linear
programs have alternative optima: the two solvers agree on every objective but choose a different reaction set in about
half of the replays, so which gene is called can depend on the solver.

## Limits

* Pilot on 22 genomes and 148 cases; recovery intervals are wide.
* Harm uses only experimental labels (Price-148 plus 35 Swiss-Prot labels); 99 % of calls are unconfirmed and are not
  counted either way.
* Results depend on the gapsmith universal database, the Biolog-to-ModelSEED mapping (77 % of substrate wells usable)
  and the default medium.
