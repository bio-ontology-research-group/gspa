# Six follow-up experiments on Stage 2B and off-growth enzyme discovery

Status: complete (2026-10-31). Pre-registered design: `tasks/v3-stage2b-followups-design-2026-10-31.md` (sha256
`5a468ece…`, frozen before any run; deviations in its §9). Code: `benchmark/v3/followups/` (16 unit tests), Shaheen
jobs `benchmark/v3/shaheen/followup-*.sbatch`, outputs `/scratch/project/k10018/gspa-v3/out/stage2b-followups-v1/`.
Results: `benchmark/v3/results/stage2b-followups-v1/` (`SHA256SUMS`). Production options added for experiments 4 and 5
are off by default and are in release candidate 3 (`/project/k10018/gspa-v3/release-candidate-3/`, JAR sha256
`8391d73f…`, full JVM gate 1,583 tests, 0 failures, 6 skipped). No language model is used anywhere; unconfirmed labels
are never negatives.

## Summary

| Experiment | Pre-registered criterion | Outcome |
|---|---|---|
| 1. Competitive gene assignment | ≥ 35 of 54 re-added truth reactions assigned to the hidden protein, harm ≤ 0.85, shuffled-gene control ≤ 2 | Not met: 25 of 54 (baseline 23), harm 0.851 |
| 2. Unblocking / dead-end tests as gated evidence | top-1 + 8 on Price-148 for ≥ 4 of 6 predictors, harm among gated updates < 0.3, no CAFA F-max loss | Not met: no predictor gains top-1 (CLEAN −2), harm 0.97–0.999 |
| 3. Per-genome transporter panel | recovery ≥ full panel at ≤ 25 % of its tests, ≥ 9 of 148 Price cases | Not met: 9 of 10 at 40 % of the tests, 8 of 148 |
| 4. Synteny covariate for Stage 2B | covariate non-empty in every cell, replay reproduces stage 2B, no Price metric below stage 2B | Met as a correctness fix; no Price change in any predictor |
| 5. Candidate limit for dense predictors | < 1,000 changed scores per cell, 0 tie artefacts, no Price metric below stage 2B | Met for the top-10 limit (the 0.05 floor is inert): CLEAN's Stage 2B harm is undone, no gain |
| 6. RB-TnSeq fitness as a held-out test | knockout readout above chance in ≥ 11 genomes; competitive assignment better than baseline | Not met: the gap-filled model grows on none of the sole-source wells; readouts at chance |

## 1. Competitive gene assignment for off-growth calls

The stored off-growth outputs (180 cases: 148 Price and 32 Swiss-Prot, union of the sole-source, unblocking and
dead-end tests) were re-assigned under five rules without new linear programs. R0 is today's rule (top-scoring
candidate). R1 is the production responsibility rule (sequence odds times exp(1.5 × synteny flag), null gene at
0.05); the synteny flag is 1 when a gene within ten loci carries a candidate for another enzyme of the same KEGG module
or is a transporter of the test's sole source. R2 is R1 with the null gene at 1.0, R3 adds paralog competition (one
reaction per enzyme number per gene), R4 restricts R0 to synteny-flagged candidates.

| Rule | Correct among 54 re-added truth reactions | Recovered, 180 cases (95 % interval) | Gene calls | No-gene calls | Supported | Harmful | Harm rate (95 % interval) | Gained / lost vs R0 (one-sided p) |
|---|---:|---|---:|---:|---:|---:|---|---|
| R0 baseline | 23 | 23 (0.074–0.176) | 130,366 | 0 | 190 | 1,013 | 0.842 (0.776–0.930) | — |
| R1 responsibility | 25 | 25 (0.085–0.187) | 130,366 | 0 | 200 | 1,145 | 0.851 (0.780–0.936) | 2 / 0 (0.25) |
| R2 strict null | 22 | 22 (0.074–0.169) | 64,251 | 66,115 | 169 | 606 | 0.782 (0.712–0.894) | 2 / 3 (0.81) |
| R3 paralog competition | 25 | 25 (0.085–0.187) | 130,366 | 0 | 155 | 1,147 | 0.881 (0.816–0.954) | 2 / 0 (0.25) |
| R4 synteny only | 22 | 22 (0.071–0.169) | 130,366 | 0 | 152 | 1,146 | 0.883 (0.821–0.946) | 2 / 3 (0.81) |

| Control | Recovered |
|---|---|
| R1, shuffled candidate genes (seeds 1–5) | 0, 0, 0, 0, 0 |
| R3, shuffled candidate genes (seeds 1–5) | 0, 0, 0, 0, 0 |
| R1, shuffled genome layout (seeds 1–5) | 25, 21, 21, 24, 22 |

The responsibility rule gains two cases (a Phaeobacter sugar-phosphate enzyme and a Shewanella enzyme, both near
ties where the hidden protein scored 0.52–0.53 against a paralog at 0.57–0.59) and loses none, far from the 35
required; the gain is not distinguishable from chance and the shuffled-layout control recovers as many cases, so the
synteny flag itself contributes nothing measurable. Harm stays at the pilot's level. The strict null gene removes half
of the calls and lowers harm to 0.78 but also loses three recoveries. The 0.26-versus-0.57 cases named in the
investigation (acyl-CoA dehydrogenase families) stay with the paralog under every rule. Cost: one node, six minutes.

## 2. Condition-free targeted tests as a gated evidence source

The stored unmasked unblocking (A2) and dead-end (A3) runs add 757–1,097 reactions per genome. For each predictor
cell, every candidate of an added reaction (stage-0 score above the sidecar's 1e-5 floor) receives the production
update (γ = 2.2, π = 1, responsibility by rule R1, reaction probability = best candidate score) when it passes the null
gate and the pathway-coherence gate (the reaction shares a non-currency metabolite with a reaction the draft model
assigns to a gene within ten loci). The sparse updates were scored on the 148 Price proteins exactly as in experiments
4 and 5.

| Predictor | Gated updates, A2 ∪ A3 (per 22 genomes) | Updates on labelled proteins: supported / harmful | Top-1 stage 0 → A2 → A3 → A2 ∪ A3 → without gate | Top-3 stage 0 → A2 ∪ A3 → without gate | Top-5 stage 0 → A2 ∪ A3 → without gate | Gained / lost top-1 (A2 ∪ A3) | CAFA F-max stage 0 → A2 ∪ A3 |
|---|---:|---|---|---|---|---|---|
| CLEAN | 3,027,330 | 4 / 5,327 | 60 → 58 → 59 → 58 → 58 | 87 → 81 → 78 | 95 → 84 → 79 | 0 / 2 | 0.409 → 0.409 |
| DeepProZyme | 2,732,479 | 4 / 4,632 | 46 → 46 → 46 → 46 → 46 | 90 → 90 → 90 | 101 → 101 → 101 | 0 / 0 | 0.380 → 0.380 |
| EnzBERT | 112,096 | 3 / 153 | 53 → 53 → 53 → 53 → 52 | 76 → 77 → 66 | 84 → 83 → 73 | 0 / 0 | 0.447 → 0.445 |
| GraphEC | 2,982,714 | 3 / 5,325 | 18 → 18 → 18 → 18 → 18 | 37 → 37 → 37 | 52 → 52 → 52 | 0 / 0 | 0.216 → 0.216 |
| MAPred | 28,514 | 3 / 80 | 50 → 50 → 50 → 50 → 46 | 95 → 94 → 94 | 99 → 99 → 99 | 1 / 1 | 0.388 → 0.388 |
| TopEC | 1,848 | 0 / 3 | 5 → 5 → 5 → 5 → 5 | 5 → 5 → 5 | 5 → 5 → 5 | 0 / 0 | 0.035 → 0.035 |

| Control | Result |
|---|---|
| Shuffled candidate genes, seed 1, A2 ∪ A3 (top-1 per predictor) | CLEAN 57, DeepProZyme 44, EnzBERT 53, GraphEC 18, MAPred 48, TopEC 5 (losses of 2–3 in three predictors) |
| Stage 0 → stage 0 no-op | 0 corrections, 0 harms for every predictor |
| Reachable ceiling | the unmasked union re-added the true reaction of 10 of 180 cases in the pilot |

No predictor gains top-1 (MAPred trades one gain for one loss; CLEAN loses two), so the +8 bar is missed by every
predictor, as the ceiling made likely. The harm among gated updates on labelled proteins is 0.97–0.999: a condition-free
test that makes a blocked reaction carry flux or drains a dead end demands a reaction, and the update then goes to
whichever genes carry any score for its enzyme numbers, which for the dense predictors are thousands of genes per
reaction. The coherence gate removes 95 % of the updates and all of the larger losses seen without it (MAPred 50 → 46,
EnzBERT top-3 76 → 66), but what remains still does not help. The shuffled-gene control loses as much as the real
assignment, which confirms that the surviving updates carry no gene-specific signal on these proteins. Cost: one node,
46 minutes for the evidence, two node-hours for the scoring.

## 3. Condition-matched panels without labels

For each genome the panel P1 keeps the Biolog wells whose compound has a gapsmith transporter hit in that genome
(105–124 of 293 wells per genome; `panels.json` sha256 `0b920012…`, written before scoring). P2 are size-matched random
panels. Everything is re-scored from the stored sole-source test results with today's assignment rule; no new linear
programs. BacDive could not be used: the repository's BacDive manifest covers none of the 22 species and the BacDive
API needs credentials.

| Panel | Tests (fraction of P0) | Recovered, 180 cases (95 % interval) | Price-148 recovered | Reaction re-added | Harm rate | Lost vs P0 (two-sided p) |
|---|---:|---|---:|---:|---:|---|
| P0 full Biolog panel | 52,740 (1.00) | 10 (0.021–0.090) | 9 | 30 | 0.859 | — |
| P1 transporter panel | 21,172 (0.40) | 9 (0.020–0.079) | 8 | 21 | 0.830 | 1 (1.0) |
| P2 random, same size (seeds 1–5) | 21,172 (0.40) | 6, 5, 7, 4, 7 (mean 5.8) | 5, 4, 6, 3, 6 | 13–19 | 0.83–0.92 | 3–6 |
| P0 odd / even columns | 0.49 / 0.51 | 6 / 8 | 5 / 7 | 16 / 21 | 0.85 / 0.87 | 4 / 2 |
| P1 odd / even columns | 0.20 / 0.20 | 6 / 7 | 5 / 6 | 9 / 17 | 0.75 / 0.84 | 4 / 3 |
| P1, shuffled candidate genes (seeds 1–5) | — | 0 in every seed | 0 | 21 | 1.0 | 10 |

The transporter panel keeps nine of the ten recoveries (losing one Sinorhizobium case) and beats every size-matched
random panel, so transporter calls do select the informative wells. It needs 40 % of the tests, not the 25 % required,
and recovers 8 Price cases, one short of the pilot's 9; the criterion is not met. Either panel half alone loses
recoveries, which shows how thin the signal is: ten recoveries spread over thirty re-added reactions. Cost: one node,
six minutes.

## 4. Synteny covariate for Stage 2B and 5. Candidate limit for dense predictors

The 132 stored Stage 2B cells were re-attributed with the production sidecar code (imported; the stored repair result
is reused, no linear program is solved) under seven arms, and the sparse score updates were scored against stage 0 and
the stored stage 2B with the frozen Price metrics. The replay arm (empty covariates, original candidate floor)
reproduces the stored stage 2B scores to within 1e-6 in all 132 cells, and the command-line path of the sidecar
(`--repair-json`, `--layout-tsv`, `--kegg-modules-tsv`, `--candidate-min-score 0.05`, `--candidate-top-k 10`) gives the
same updates as the function-level run in the 22 MAPred cells where it was exercised.

| Predictor | Score floor | Repaired reactions per cell (median) | Changed scores per cell, stage 2B as stored (min / median / max) | Cells ≥ 1,000 | Top-call changes (tie artefacts) | Synteny covariates per cell (min / median / max) | Changed scores with top-10 limit (min / median / max) | Top-call changes with top-10 limit |
|---|---:|---:|---|---:|---:|---|---|---:|
| CLEAN | 0.0775 | 22 | 11,682 / 44,419 / 86,386 | 22 | 37,544 (124) | 135 / 239 / 330 | 50 / 150 / 300 | 1,797 |
| DeepProZyme | 0.0775 | 8 | 11,687 / 33,678 / 56,098 | 22 | 0 (0) | 219 / 300 / 583 | 40 / 80 / 220 | 24 |
| GraphEC | 0.0894 | 114 | 37,898 / 112,535 / 212,076 | 22 | 1 (0) | 2 / 197 / 561 | 240 / 920 / 970 | 48 |
| EnzBERT | 0.0775 | 26 | 899 / 1,766 / 4,729 | 20 | 653 (0) | 131 / 206 / 331 | 126 / 170 / 340 | 1,406 |
| MAPred | 0.0803 | 5 | 1 / 54 / 141 | 0 | 6 (0) | 360 / 551 / 1,227 | 1 / 42 / 108 | 5 |
| TopEC | 0.106–0.124 | 113 | 3 / 6 / 12 | 0 | 0 (0) | 103 / 167 / 359 | 3 / 6 / 12 | 0 |

Two facts follow before any Price scoring. The 0.05 candidate floor changes nothing for any predictor, because every
predictor's own lowest exported score (0.0775–0.124) already lies above it; the dilution comes from thousands of
candidates at or just above the predictor floor, which only the top-10 limit removes (changed scores fall from tens
of thousands to 40–970 per cell, under the 1,000 bar in every cell, with zero tie artefacts). The synteny covariate is
non-empty in every cell (2–1,227 gene–enzyme pairs), but with unlimited candidates it changes the weights without
changing which scores move, so its changed-score counts equal the stored stage 2B's.

Price scoring (exact enzyme numbers, 148 proteins; arms with identical rows are merged):

| Predictor | Arm | Top-1 | Top-3 | Top-5 | Mean reciprocal rank | Top call changed vs stage 0 | Gained / lost vs stage 0 | CAFA F-max (95 % interval) |
|---|---|---:|---:|---:|---:|---:|---|---|
| CLEAN | stage 0 | 60 | 87 | 95 | 0.516 | — | — | 0.409 (0.330–0.506) |
| CLEAN | stage 2B as stored = replay = synteny = synteny with shuffled layout = floor 0.05 | 59 | 80 | 84 | 0.489 | 10 | 0 / 1 | 0.409 (0.330–0.506) |
| CLEAN | top-10 limit (K2) = floor + top-10 (K3) = synteny + top-10 | 60 | 87 | 95 | 0.516 | 0 | 0 / 0 | 0.409 (0.330–0.506) |
| DeepProZyme | all arms | 46 | 90 | 101 | 0.479 | 0 | 0 / 0 | 0.380 (0.335–0.458) |
| EnzBERT | stage 0 | 53 | 76 | 84 | 0.443 | 0 | — | 0.447 (0.357–0.540) |
| EnzBERT | stage 2B as stored = replay = synteny = shuffled layout = floor 0.05 | 53 | 75 | 82 | 0.441 | 0 | 0 / 0 | 0.447 (0.356–0.540) |
| EnzBERT | top-10 limit (K2, K3, synteny + top-10) | 53 | 75 | 84 | 0.442 | 0 | 0 / 0 | 0.447 (0.357–0.540) |
| GraphEC | all arms | 18 | 37 | 52 | 0.217 | 0 | 0 / 0 | 0.216 (0.104–0.319) |
| MAPred | all arms | 50 | 95 | 99 | 0.494 | 0 | 0 / 0 | 0.388 (0.321–0.472) |
| TopEC | all arms | 5 | 5 | 5 | 0.034 | 0 | 0 / 0 | 0.035 (0.007–0.075) |

Micro F-max, pair-level corrections and harms at the 0.5 threshold are identical across arms for every predictor
(`summary_exact.tsv`). Criterion 4 is met as a correctness fix: the covariate is present in every cell, the replay
reproduces the stored stage 2B, and no Price metric is below stage 2B's; there is no Price gain, because the synteny
weights only redistribute responsibility among candidates that already received it, and with unlimited candidates the
Price enzymes still see nothing (the shuffled-layout control is identical to the real layout). Criterion 5 is met for
the top-10 limit: changed scores fall under 1,000 in every cell, tie artefacts disappear, and the only Stage 2B
harm on Price-148, CLEAN's loss of one top-1 and seven top-3 enzymes, is undone (CLEAN returns to stage 0 exactly;
EnzBERT recovers its top-5). The 0.05 floor is inert for the reason above. The limit removes harm; it does not create
a gain, since the repaired reactions still do not touch the Price enzymes.

## 6. RB-TnSeq fitness as a held-out test of condition-specific essentiality

Fitness data: the public Fitness Browser database and protein sequences (`feba.db`, 10.4 GB, sha256 `e033a1b1…`;
`aaseqs`; downloaded 2026-10-31, data only). All 22 assemblies were matched to a Fitness Browser organism by exact
protein-sequence hash (table `benchmark/v3/followups/data/price22_fitness_browser_organisms.tsv`, written before
scoring); 1,400 carbon/nitrogen/phosphorus/sulfur-source experiments were matched to Biolog wells by compound name
(1,495 experiment–well pairs) and 337 were not (listed in `experiments_unmatched.tsv`). Five genomes fell below the 80 % model-gene match gate and are
excluded (one is a different strain from its Fitness Browser namesake: 5 % match; four are at 75–79 %).

| Model variant | Pairs (gene × experiment), model genes | Positives | Pooled AUC (95 % interval) | Pairs, added-reaction genes | Positives | Pooled AUC (95 % interval) | Permuted-label AUC (mean of 5 seeds) |
|---|---:|---:|---|---:|---:|---|---:|
| M0 gapsmith gap-filled model | 0 | 0 | — (wild type grows on no sole-source well) | 0 | 0 | — | — |
| M1 + union added reactions, baseline assignment | 811,283 | 2,022 | 0.511 (0.498–0.527) | 328,318 | 759 | 0.533 (0.507–0.562) | 0.506 / 0.500 |
| M2 + same reactions, responsibility assignment | 818,805 | 2,041 | 0.512 (0.500–0.528) | 335,840 | 778 | 0.535 (0.507–0.566) | 0.500 / 0.508 |
| M3 random assignment (seeds 1–5) | 575,592 | 1,637 | 0.500–0.519 | 92,627 | 374 | 0.545–0.557 | 0.499–0.520 |

| Criterion part | Result |
|---|---|
| Knockout readout above chance in ≥ 11 of 17 usable genomes (M0) | not testable: M0 wild type grows on 0 of 4–64 matched wells in every genome; M1 above chance in 1 genome, M2 in 3 |
| Responsibility assignment better than baseline on added-reaction genes (M2 − M1) | +0.002 (−0.007 to +0.013): not met |
| Cost | 4.69 million linear programs, 16 CPU hours, about one minute per genome on one node |

The gapsmith gap-filled model is filled for its default medium only and cannot take up most Biolog compounds, so the
pre-registered readout on it cannot be evaluated (reported as a limit, not re-defined). The variants that carry the
union of the targeted tests' added reactions do grow, but their knockout predictions agree with the specific
phenotypes only at chance level: pooled AUC 0.51 with the interval touching 0.5, and random gene assignment (M3)
scores as well as either rule on the added-reaction genes. The fitness data therefore give no support to the gene
assignments, and no support for the responsibility rule over the baseline.

## What this shows

None of the six follow-ups turns the Stage 2B null result or the off-growth pilot's weak gene assignment into a gain
on the Price enzymes. The two mechanism fixes (synteny covariate present, candidate limit for dense predictors) are
correct and are kept, off by default, in release candidate 3: the top-10 limit removes the only Stage 2B harm seen on
Price-148 without adding anything. The competitive assignment rule, the transporter-selected panel and the
condition-free evidence source all do what they were designed to do in a small way (two more cases, nine of ten
recoveries at 40 % of the tests, 95 % fewer updates) and none reaches its pre-registered bar; the harm rates on labelled
proteins stay at 0.78–0.999. The RB-TnSeq fitness data, used only as a test, give no support to the gene assignments
produced by any rule, and the gapsmith gap-filled models cannot be tested on sole-source wells at all. The limits on
what growth-based tests can say about these enzymes are structural, not a matter of weights.

## Cost and outputs

| Item | Value |
|---|---|
| Shaheen jobs | 37 (`gspa-followup-*`, 15730704–15730883; two cancelled by me and resubmitted, one failed on a missing module and resubmitted) |
| Node-hours (exclusive 192-core nodes) | 9.7 (about 1,850 physical-core-hours; 3,700 as Slurm counts hardware threads) |
| Linear programs solved | 4.69 million (experiment 6 only; experiments 1–5 solve none) |
| Outputs on Shaheen | `/scratch/project/k10018/gspa-v3/out/stage2b-followups-v1/` 7.4 GB (sparse update records, diagnostics, scoring reports); `/project/k10018/gspa-v3/fitness-browser-v1/` 399 MB (prepared fitness tables); `/project/k10018/gspa-v3/release-candidate-3/` 46 MB |
| Results in the repository | `benchmark/v3/results/stage2b-followups-v1/` (73 MB, `SHA256SUMS`) |
| Fitness Browser download | `feba.db` 10.4 GB and `aaseqs` 118 MB in `~/tmp/fitbrowser/` on the workstation (not copied into the repository) |

## Limits

* Experiments 1–3 re-score stored outputs; the added reactions themselves are the pilot's and depend on its solver,
  medium and universal database.
* Experiment 2 is a benchmark-level implementation of the proposed sidecar mode (same update formula and constants);
  production wiring is deferred because the criterion is not met.
* Experiments 4 and 5 recompute only the attribution; the repair is the stored one.
* Experiment 6 scores 17 of 22 genomes; specific phenotypes are rare (0.25 % of pairs), condition-name matching misses
  23 % of experiments, and the fitness libraries cannot test genes that are essential in rich medium.
