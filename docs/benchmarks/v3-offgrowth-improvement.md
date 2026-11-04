# Off-growth enzyme discovery, second round: measured on held-out curated labels

Status: tier A complete (2026-11-03). Pre-registered plan: `tasks/v3-row12-offgrowth-improvement-2026-11-02.md`,
frozen at sha256 `36af054b8e3a010f4c86dbe58660098158ccf86449b9553e85de78e2ef404f83`
(`benchmark/v3/results/offgrowth-improvement-v1/PLAN_SHA256`) before any held-out score was read; deviations in
its §13. Exploration that motivated the design: `docs/benchmarks/v3-offgrowth-improvement-exploration.md`.
First round: `docs/benchmarks/v3-offgrowth-discovery.md`. Code: `benchmark/v3/offgrowth2/` (24 unit tests),
Shaheen job `benchmark/v3/shaheen/offgrowth2-genome.sbatch`. Results:
`benchmark/v3/results/offgrowth-improvement-v1/`. No language model is used anywhere; unconfirmed labels are never
negatives.

## What was asked and what was done

The first round showed that targeted growth simulations recover enzymes that biomass gap-filling misses, but its
calls were mostly wrong and its label set was small. This round changes three things at once: sole-source tests are
run **differentially** from a baseline minimal-medium fill, so a compound's fill is its own catabolic steps rather
than a rebuild of amino-acid biosynthesis; a reaction is called only when some test adds it inside a **small,
solver-stable fill**; and a gene is named only when the sequence ranking is **unambiguous**, otherwise the call
stays at reaction level. The arms are scored against 643 curated enzyme reannotations of the Price laboratory
(Fitness Browser) that were never used in the first round, in the same 22 bacteria.

## Set-up

| Item | Value |
|---|---|
| Genomes | the 22 Price bacteria, gapsmith draft and gap-filled models (`benchmark-models-v1`, October recipe) |
| Masked cases | 823 (148 Price-148, 32 Swiss-Prot experimental, 643 Fitness Browser reannotations); 800 off the default growth path, 23 on it |
| Masking audit | passes for every case; 0 cases skipped |
| Candidate genes | union of the six enzyme predictors, best score per protein and enzyme (as in the first round) |
| Labels for the pre-registered threshold (seen) | Price-148 + Swiss-Prot experimental (L0) |
| Labels for the test (held out, outcomes never read before the freeze) | Fitness Browser EC reannotations, 643 rows on 615 proteins (L1) |
| Simulation budget | one exclusive node per genome; 22.4 node-hours for the wave, about 8,400 linear programs and 7,000 CPU seconds per run |

## Results

### Recovery (masked cases, off the growth path)

| Arm | Gene recovery, 95 % interval | Reaction recovery | Pathway recovery | L0 cases (178) | L1 cases (622) |
|---|---|---|---:|---|---|
| S0 biomass gap-filling | 7/800 = 0.009 [0.00, 0.02] | 11/800 = 0.014 | 0 | 0 | 7 |
| S1 first-round procedure | 241/800 = 0.301 [0.28, 0.33] | 342/800 = 0.427 | 0 | 23 | 218 |
| S1b same at the new budget | 231/800 = 0.289 [0.26, 0.32] | 321/800 = 0.401 | 0 | 23 | 208 |
| **N1 new procedure, transporter panel** | **91/800 = 0.114 [0.09, 0.14]** | 204/800 = 0.255 | 20/800 | 5 | 86 |
| N2 new procedure, measured-growth panel | 90/800 = 0.113 [0.09, 0.14] | 202/800 = 0.253 | 19/800 | 5 | 85 |
| N1 without the fill gate | 124/800 = 0.155 | 286/800 = 0.357 | 0 | 12 | 112 |
| N1 without the margin and free-protein rule | 161/800 = 0.201 | 204/800 = 0.255 | 20/800 | 11 | 150 |
| N1 without the stability requirement | 102/800 = 0.128 | 237/800 = 0.296 | 32/800 | 8 | 94 |
| N1 at fill ≤ 1 | 49/800 = 0.061 | 138/800 = 0.172 | 0 | 3 | 46 |
| N1 at fill ≤ 5 | 101/800 = 0.126 | 224/800 = 0.280 | 26/800 | 7 | 94 |
| N1 with the common block removed | 94/800 = 0.117 | 210/800 = 0.263 | 22/800 | 5 | 89 |

The first round's 23 recoveries on its own 178 cases are reproduced exactly by S1.

### Precision and harm on distinct labelled decisions (unmasked run, one per genome)

| Arm | Held-out labels: supported | harmful | precision, 95 % interval | harm | Seen labels: precision (decisions) | Unconfirmed calls |
|---|---:|---:|---|---:|---|---:|
| S0 | 4 | 5 | 0.444 [0.00, 0.80] | 0.56 | 0.000 (3) | small |
| S1 | 84 | 181 | 0.317 [0.27, 0.36] | 0.68 | 0.179 (67) | large |
| S1b | 77 | 165 | 0.318 [0.27, 0.36] | 0.68 | 0.188 (64) | large |
| **N1** | **25** | **21** | **0.543 [0.41, 0.68]** | **0.46** | 0.800 (5) | 1,112 |
| N2 | 25 | 21 | 0.543 [0.41, 0.68] | 0.46 | 0.800 (5) | 1,112 |
| N1 without the fill gate | 35 | 35 | 0.500 [0.42, 0.59] | 0.50 | 0.375 (16) | larger |
| N1 without the margin and free-protein rule | 48 | 79 | 0.378 [0.30, 0.45] | 0.62 | 0.233 (30) | larger |
| N1 without the stability requirement | 30 | 24 | 0.556 [0.44, 0.67] | 0.44 | 0.556 (9) | 1,180 |
| N1 at fill ≤ 1 | 7 | 11 | 0.389 [0.19, 0.61] | 0.61 | 0.750 (4) | smaller |
| N1 at fill ≤ 5 | 28 | 24 | 0.538 [0.44, 0.67] | 0.46 | 0.667 (6) | 1,160 |
| Sequence evidence alone, volume-matched | 11 | 5 | 0.688 [0.53, 0.88] | 0.31 | 0.000 (1) | — |

Of N1's 21 harmful calls, 13 name an enzyme in the same sub-subclass as the protein's labelled enzyme.
Abstention: N1 leaves 85 % of its called reactions without a gene and does not call 55 % of the reactions its
admitted tests add.

### Controls

| Control | Result |
|---|---|
| Shuffled candidate genes (seeds 1–5), N1 | 0, 1, 0, 0, 0 gene recoveries |
| Size-matched random panels (seeds 1–5), N1 | 87, 89, 88, 88, 87 against N1's 91 |
| Masking audit | 823 of 823 cases pass; 0 skipped |
| Solver stability | 88 % of a fill's reactions survive the second solver; a run's differential sole-source tests have a median fill of 11 reactions and a median stable fill of 9 (medians over the 845 runs) |

### The functions the simulations add

45 of N1's gene calls land on a labelled protein that the volume-matched sequence-only comparator does not call
for that enzyme; 24 of them agree with the curated label.

| Enzyme | Genomes | Pathway it belongs to |
|---|---:|---|
| 2.7.1.92 5-dehydro-2-deoxygluconokinase | 5 | myo-inositol catabolism |
| 5.3.1.7 mannose isomerase | 3 | mannose catabolism |
| 1.13.11.5 homogentisate 1,2-dioxygenase | 2 | tyrosine catabolism |
| 1.13.11.27 4-hydroxyphenylpyruvate dioxygenase | 2 | tyrosine catabolism |
| 2.7.1.58 2-dehydro-3-deoxygluconokinase | 2 | Entner–Doudoroff from gluconate |
| 1.2.1.38 N-acetyl-glutamyl-phosphate reductase | 2 | arginine biosynthesis |
| 1.14.16.1 phenylalanine 4-monooxygenase, 1.4.1.9 leucine dehydrogenase, 3.5.1.9 arylformamidase, 3.5.1.25 N-acetylglucosamine-6-phosphate deacetylase, 3.5.4.5 cytidine deaminase, 1.2.1.88 proline dehydrogenase, 1.1.1.3 homoserine dehydrogenase, 5.3.1.14 rhamnose isomerase | 1 each | amino-acid, amino-sugar and nucleoside catabolism |

The wrong additional calls repeat across genomes too: dihydrouracil dehydrogenase (7), N-carbamoylputrescine
amidase (5) and L-rhamnose mutarotase (4) account for three quarters of them.

The 20 pathway-level recoveries are the clearest cases: in *Shewanella oneidensis* MR-1 and *Stutzerimonas
stutzeri* RCH2 the tyrosine route is recovered as a unit (4-hydroxyphenylpyruvate dioxygenase and homogentisate
1,2-dioxygenase on adjacent genes, with phenylalanine 4-monooxygenase in the first), and in *Sinorhizobium
meliloti* the inositol route (myo-inosose-2 dehydratase with the 5-dehydro-2-deoxygluconokinase above).

## The same genomes with our own candidate tables, and the two genome-level tiers

Everything above uses the six enzyme predictors' union as the candidate table, a snapshot that exists only for the
22 Price proteomes. To test genomes those predictors never saw, we ran CLEAN ourselves on 47 proteomes (1.8 GPU
hours on one H200 node; 212,722 proteins; GMM confidence, no score floor) and built candidate tables from its
output. The same 22 genomes were then re-run with these tables, together with the eight recurrence-cohort bacteria
(tier B, Swiss-Prot experimental labels) and 17 further Fitness Browser organisms outside the 145 benchmark
genomes (tier C, curated reannotations). Both new tiers are held out at the **genome** level.

| Run | Genomes | Masked cases (off-growth) | Candidate rows per genome | Shaheen |
|---|---:|---|---:|---|
| Tier A, six predictors | 22 | 823 (800) | 347,283 | 22.4 node-hours |
| Tier A, own CLEAN | 22 | 823 (800) | 3,579 | 26.4 node-hours |
| Tier B, own CLEAN | 8 | 75 (60) | 3,579 | 3.7 node-hours |
| Tier C, own CLEAN | 17 | 318 (280) | 3,579 | 11.9 node-hours |

### Recovery and precision with the own-CLEAN tables

| Run | Arm | Gene recovery | Reaction recovery | Precision on the held-out labels, 95 % interval | Harm | Gene abstention |
|---|---|---|---|---|---:|---:|
| Tier A | S0 biomass | 1/800 | 1/800 | 6/6 = 1.000 | 0.00 | 0.37 |
| Tier A | S1 first round | 179/800 = 0.224 | 308/800 = 0.385 | 71/106 = 0.670 [0.59, 0.75] | 0.33 | 0.73 |
| Tier A | **N1 new** | 97/800 = 0.121 | 195/800 = 0.244 | **24/29 = 0.828 [0.70, 0.96]** | 0.17 | 0.87 |
| Tier A | N2 measured panel | 92/800 = 0.115 | 187/800 = 0.234 | 24/28 = 0.857 [0.72, 1.00] | 0.14 | 0.87 |
| Tier A | Sequence alone | — | — | 26/29 = 0.897 [0.76, 1.00] | 0.10 | — |
| Tier B | S1 first round | 8/60 = 0.133 | 9/60 = 0.150 | 3/3 = 1.000 | 0.00 | 0.77 |
| Tier B | **N1 new** | 3/60 = 0.050 | 4/60 = 0.067 | 1/1 = 1.000 | 0.00 | 0.89 |
| Tier C | S1 first round | 107/280 = 0.382 | 165/280 = 0.589 | 58/72 = 0.806 [0.72, 0.87] | 0.19 | 0.73 |
| Tier C | **N1 new** | 42/280 = 0.150 | 85/280 = 0.304 | **16/21 = 0.762 [0.62, 0.94]** | 0.24 | 0.88 |
| Tier C | N2 measured panel | 39/280 = 0.139 | 81/280 = 0.289 | 14/18 = 0.778 [0.64, 0.92] | 0.22 | 0.88 |
| Tiers B and C pooled | **N1 new** | 45/340 = 0.132 | 89/340 = 0.262 | **17/22 = 0.773** | 0.23 | 0.88 |

Controls hold in every run: shuffled candidate genes recover at most one case in any seed, size-matched random
panels stay below the transporter panel, and the masking audit passes for all 1,216 cases with none skipped.

### The candidate table, not the simulation rule, sets precision

The two candidate tables differ by a factor of about a hundred in density, and that single difference moves
precision more than any arm or ablation in this benchmark.

| Candidate table | Rows per genome | N1 precision | First-round rule precision | Biomass gap-filling precision |
|---|---:|---|---|---|
| Six predictors' union | 347,283 | 25/46 = 0.543 | 84/265 = 0.317 | 4/9 = 0.444 |
| Our own CLEAN | 3,579 | 24/29 = 0.828 | 71/106 = 0.670 | 6/6 = 1.000 |

With the sparse table every arm, including the first round's crude rule and plain biomass gap-filling, calls
precisely; with the dense table none of them does. The gate and the abstention rule still help within each table
(0.828 against 0.670 on tier A, 0.762 against 0.806 on tier C), but they are second-order next to the choice of
what a gene may be called for at all. The dense table's promiscuity is the same effect the first round reported as
"the assignment rule is too crude for promiscuous predictor scores"; it is better described as the table being too
promiscuous for any assignment rule.

### Functions the simulations add on the sparse table

With the own-CLEAN tables, 17 of N1's tier A calls on labelled proteins name an enzyme the volume-matched
sequence-only comparator does not call, and 13 of them agree with the curated label: pyruvate dehydrogenase E1
component (3 genomes), N-acetyl-glutamyl-phosphate reductase and 2-dehydro-3-deoxygluconokinase (2 each),
homogentisate 1,2-dioxygenase, isocitrate lyase, gluconokinase, nitrite reductase, butyryl-CoA dehydrogenase and
acetylesterase. The butyryl-CoA dehydrogenase is notable: that family defeated every rule in the first round
(eight re-added reactions, no correct gene), and the sparse table resolves one of them. On tiers B and C no call
beats the sequence-only comparator, which on those tiers is fitted and scored on the same labels (below).

## Verdict against the pre-registered criteria (§8 of the plan)

Criteria 1 to 5 are judged on tier A with the six-predictor tables, which the plan names as primary; criterion 6 on
tiers B and C pooled.

| Criterion | Required | Measured (primary) | Outcome |
|---|---|---|---|
| 1. Measured precision | ≥ 0.75 on ≥ 20 decisions, harm ≤ 0.25, lower bound ≥ 0.60 | 0.543 on 46 decisions, harm 0.46, lower bound 0.41 | **not met** (met on the own-CLEAN table: 0.828 on 29, harm 0.17, lower bound 0.70) |
| 2. Some additional recovery | ≥ 25 gene recoveries, biomass gap-filling at 0, ≥ 10 beyond sequence-only | 91 recoveries, biomass gap-filling at 7, 24 beyond sequence-only | partly met |
| 3. Better than the first round at the same budget | recovery not lower than S1b, harm lower by ≥ 0.25 | 91 against 231; harm 0.46 against 0.68 | **not met** |
| 4. Pathway level | reaction recovery ≥ 30 %, ≥ 10 pathway calls containing a hidden case | 25.5 %, 20 pathway calls | partly met |
| 5. Controls | shuffled ≤ 2 per seed, random panels below N1 in ≥ 4 of 5 seeds, masking audit passes | 0–1 per seed, below in 5 of 5, audit passes on all 1,216 cases | met |
| 6. Genome-level tiers | ≥ 0.75 on ≥ 20 decisions and ≥ 10 gene recoveries, biomass gap-filling at 0 | 0.773 on 22 decisions, 45 gene recoveries, biomass gap-filling at 1 | met except the last part |

**Final verdict: the improved procedure is not accepted at its pre-registered bar.** On the primary table it
recovers some additional functions and abstains as designed, but its precision stays at about half, and its recall
at about a third, of what the criteria asked. The genome-level tiers, run on genomes no part of this work had
seen, do reach the precision level and do recover functions, so the design transfers; they reach it with a
candidate table on which the first round's cruder rule reaches it too. The honest reading is that this round
measured the right thing and found that the lever is elsewhere: for calls that are worth a curator's time, the
candidate table must be sparse, and given a sparse table the simulations add a modest, measurable set of functions
that sequence evidence alone does not name.

## What the numbers say about the mechanism

The differential design does what it was meant to do on single reactions: the sole-source fills shrink from
150–220 reactions to a median of 11 per run, and where the catabolic step has a candidate gene the enzyme comes back
inside a fill of one to three reactions. The gate built on that is what buys the precision. It is also what costs
the recall, and the ablations locate the cost precisely: removing the gate adds a third more recoveries and loses
precision, while removing the margin and free-protein rule adds three quarters more recoveries and loses more
precision still. The first round's rule names a gene for every added reaction and so wins on count; two thirds of
those names are wrong where a label exists.

The held-out labels also change the picture of the first round. The curated reannotations are enzymes the Price
laboratory assigned **because** a mutant had a specific growth phenotype, so they are exactly the functions a
growth simulation can reach: the first round's procedure recovers 35 % of them against 13 % of its own label set.
Biomass gap-filling, which recovered nothing in the first round, recovers seven of them.

Measured growth conditions do not help. With the measured-growth panel built from 926 matched Fitness Browser
carbon, nitrogen, phosphorus and sulfur experiments, the panel is a third of the size of the transporter panel and
the result is the same to within one case, which says the informative wells are already selected by the
transporter annotation.

## Limits

* Tier A is held out at the **label** level only: the 22 genomes, their models, their candidate tables and the
  simulations themselves were seen in the first round. Tiers B and C, which are held out at the genome level, are
  not yet run.
* The curated reannotations carry phenotype and homology evidence together; 66 of the 643 rows have a hedged
  curator comment and are kept, not dropped.
* Seven closely related *Pseudomonas* genomes hold about half of the held-out labels; the split tables in
  `arms.tsv` report them separately, and N1's recovery is lower there than in the other fifteen genomes.
* Reaction-level and pathway-level calls have no negative labels, so they are reported with recovery and volume
  only.
* Precision rests on 46 labelled decisions on the primary table and 29 on the own-CLEAN table; the intervals are
  wide, and the comparison with sequence evidence alone rests on 16 and 29 decisions.
* On tiers B and C the threshold and test label sets are the same file, because those tiers have no separate seen
  label set. The sequence-only comparator's fixed threshold is therefore fitted on the labels it is scored on and
  is reported as descriptive only; the arms themselves fit nothing.
* Tier B has no measured growth conditions (its eight bacteria are not Fitness Browser organisms), so its
  measured-panel arm falls back to the transporter panel in every run and is a duplicate there. In tier C only one
  organism lacks measured conditions.
* Two proteins were scored by CLEAN from sanitised sequences, with one unreadable character each replaced by X
  (an ambiguous residue in one *Haemophilus influenzae* protein, thirteen internal stop codons in one
  *Xanthomonas campestris* protein). Their candidate rows are kept unchanged.
* The own-CLEAN candidate tables keep only calls with a confidence above zero and a complete four-level enzyme
  number, as the plan defines them; this drops about a quarter of CLEAN's raw calls, which are the zero-confidence
  and preliminary-number rows.
* Tier C's organisms were chosen because they carry curated reannotations, so their labels are again enzymes
  assigned from growth phenotypes; tier B's are Swiss-Prot experimental enzymes and are a different, harder kind of
  target, which its small counts cannot separate.

## Cost and outputs

| Item | Value |
|---|---|
| Shaheen jobs | 28 (22 genomes plus six pilot and tier C model builds), 22.4 node-hours on exclusive 192-core nodes |
| Linear programs | 7,093,806 across the 845 runs; 1,642 CPU hours |
| Scoring | 7 minutes on the workstation, 22 arms per run, 2,000-resample genome bootstrap |
| Outputs | `/scratch/project/k10018/gspa-v3/out/offgrowth2-v1/` (per-genome tars); tables in `benchmark/v3/results/offgrowth-improvement-v1/` |
| All four runs together | 94 Shaheen jobs, 64.4 node-hours (about 12,400 physical core-hours), 1,216 masked cases; 1.8 GPU hours on orix for the CLEAN inference |
