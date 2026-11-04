# Finding enzymes that growth on one medium cannot reveal: full benchmark panel

Status: complete (2026-11-02). Pre-registered design `tasks/v3-offgrowth-discovery-design-2026-10-30.md` (deviations
D1–D7); addendum for this run `tasks/v3-offgrowth-full-panel-2026-11-02.md`, frozen before any scored check (sha256
`a7be7111…`, `benchmark/v3/results/offgrowth-full-v1/plan.sha256`). Pilot: `docs/benchmarks/v3-offgrowth-discovery.md`.
Results: `benchmark/v3/results/offgrowth-full-v1/` (`report/report.md` is the runner's full report, `report/tables/`
the secondary tables, `SHA256SUMS`, large archives listed in `LARGE_FILES.tsv` and kept on Shaheen). No language model
is used anywhere.

## What was run

* Models: the gapsmith draft and gap-filled models of all 145 benchmark genomes (22 Price, 18 recurrence cohort,
  105 recurrence peers; 79 bacteria, 66 archaea), built on 2026-10-30 with the pilot's recipe; 145 of 145 present and
  matching their proteome hashes, so none had to be built.
* Hidden cases: Price-148 (148 cases) plus every protein of the 145 proteomes with a Swiss-Prot experimental complete
  EC, found with the pilot's unchanged rule (reviewed entries of the organism, experimental evidence codes, exact
  full-sequence match): 425 labels on 372 proteins in 33 genomes. On the 22 Price genomes the labels are identical to
  the pilot's 35. In all, 520 cases in 49 genomes; the other 96 genomes have no experimentally labelled protein and
  could not be scored, so no tests were run on them.
* Candidate genes: the six-predictor union exists only for the 22 Price genomes, so every genome used the Swiss-Prot
  DIAMOND homology baseline (the training adapter and parameters; score = identity × coverage), EC claims only.
  Its reference was rebuilt with a repaired genus fence (below). The 22 Price genomes were re-run with this source;
  the pilot gives the comparison with the six-predictor union on them.
* Tests, masking, call rule, statistics and controls exactly as pre-registered (A0 biomass; A1 sole sources; A2
  unblocking; A3 dead ends; union; A4 rules at 10 % and 25 % of tests; GLPK with HiGHS replay on 5 % of fills;
  genome-block bootstrap; shuffled genes, unmasked run).
* Shaheen: one candidate job, 49 genome jobs on 49 nodes in parallel (2.6 to 38.8 minutes each, 9.5 node-hours),
  one aggregation job.

## Leakage and the repaired fence

The old fence compared only the first word of the Swiss-Prot organism name. The new one also drops an entry when its
taxid is a panel organism or species taxid, or when any name in its Swiss-Prot lineage is a fenced genus.

| Reference | Entries | Dropped by taxid | Dropped by genus | of which only by lineage | Hidden proteins with an identical entry | with a ≥ 95 % identical entry | Highest identity |
|---|---:|---:|---:|---:|---:|---:|---:|
| old training reference | 82,364 | — | — | — | 13 | 13 | 100 % |
| repaired fence (used) | 82,233 | 13,539 | 44,165 | 1,048 | 0 | 0 | 89.1 % |

The 13 leaks in the old reference were all *Sinorhizobium meliloti* proteins whose Swiss-Prot entries are named
"Rhizobium meliloti". An independent pass over the flat file found no remaining entry with a panel taxid or fenced
genus (`inputs/fence_check.json`).

## Results

### Recovery of hidden enzymes (520 cases; 466 off the default growth path)

| Test | Reaction re-added (off-growth) | Protein recovered (off-growth) | Recovery, 95 % interval | Protein recovered (54 on-growth cases) |
|---|---:|---:|---|---:|
| A0 biomass, default medium | 2 | 2 | 0.004 [0.000, 0.011] | 22 |
| A1 sole sources | 52 | 24 | 0.052 [0.032, 0.071] | 32 |
| A2 unblock | 44 | 16 | 0.034 [0.020, 0.048] | 30 |
| A3 dead ends | 38 | 15 | 0.032 [0.018, 0.046] | 21 |
| A1 + A2 + A3 | 83 | 34 | 0.073 [0.049, 0.097] | 32 |
| A4 rules, 10 % | 39 | 17 | 0.036 [0.019, 0.057] | 32 |
| A4 rules, 25 % | 48 | 20 | 0.043 [0.024, 0.063] | 32 |

### By subset (union A1 + A2 + A3, off-growth cases)

| Subset | Genomes | Cases | Recovered | Recovery, 95 % interval |
|---|---:|---:|---:|---|
| Price genomes | 22 | 178 | 7 | 0.039 [0.011, 0.065] |
| Recurrence cohort | 17 | 251 | 24 | 0.096 [0.062, 0.128] |
| Recurrence peers | 10 | 37 | 3 | 0.081 [0.000, 0.231] |
| Bacteria | 30 | 241 | 14 | 0.058 [0.030, 0.086] |
| Archaea | 19 | 225 | 20 | 0.089 [0.052, 0.123] |
| Price-148 cases | 22 | 146 | 4 | 0.027 [0.007, 0.053] |
| Swiss-Prot cases | 33 | 320 | 30 | 0.094 [0.062, 0.126] |

All arms: `report/tables/subsets.tsv`; per genome: `report/tables/per_genome.tsv`.

### Pre-registered comparisons (off-growth cases)

| Comparison | Test only | A0 only | p (exact McNemar, one-sided) | Holm-adjusted | Outcome |
|---|---:|---:|---:|---:|---|
| H1: A1 + A2 + A3 vs A0 | 32 | 0 | 2.3e-10 | — | supported |
| H2: A1 vs A0 | 22 | 0 | 2.4e-07 | 7.2e-07 | supported |
| H2: A2 vs A0 | 14 | 0 | 6.1e-05 | 1.2e-04 | supported |
| H2: A3 vs A0 | 13 | 0 | 1.2e-04 | 1.2e-04 | supported |
| H3: A4 at 25 % keeps ≥ 80 % of the union | 20 of 34 = 0.59 [0.42, 0.75] | | | | not met |

### Calls and harm (Price-148 and all Swiss-Prot experimental labels)

| Test | Gene calls | Supported | Harmful | Unconfirmed | Harm rate, 95 % interval | Harm gate (> A0 + 0.10) | Distinct labelled decisions | of which harmful |
|---|---:|---:|---:|---:|---|---|---:|---:|
| A0 | 25,961 | 116 | 273 | 25,572 | 0.70 [0.54, 0.88] | — | 44 | 17 (0.39) |
| A1 | 93,172 | 383 | 1,565 | 91,224 | 0.80 [0.74, 0.85] | exceeded | 154 | 76 (0.49) |
| A2 | 95,093 | 267 | 1,287 | 93,539 | 0.83 [0.69, 0.90] | exceeded | 140 | 74 (0.53) |
| A3 | 64,511 | 253 | 941 | 63,317 | 0.79 [0.57, 0.90] | no | 107 | 51 (0.48) |
| A1 + A2 + A3 | 143,357 | 509 | 2,301 | 140,547 | 0.82 [0.72, 0.86] | exceeded | 216 | 117 (0.54) |
| A4 at 10 % | 72,373 | 287 | 1,149 | 70,937 | 0.80 [0.75, 0.87] | no | 124 | 55 (0.44) |
| A4 at 25 % | 87,620 | 310 | 1,313 | 85,997 | 0.81 [0.75, 0.87] | exceeded | 140 | 63 (0.45) |

The first columns add calls over the 520 cases, as the design specifies, so one call on near-identical models counts
once per case; the last two count each (genome, protein, reaction) once.

### Cost (520 cases)

| Test | Tests run | Linear programs | Gap-fills | CPU hours |
|---|---:|---:|---:|---:|
| A0 | 520 | 2,600 | 520 | 0.9 |
| A1 | 152,360 | 451,298 | 126,614 | 103.9 |
| A2 | 395,166 | 1,564,808 | 260,349 | 176.5 |
| A3 | 164,546 | 267,952 | 93,644 | 54.7 |
| A1 + A2 + A3 | 712,072 | 2,284,058 | 480,607 | 335.1 |
| A4 at 10 % | 71,447 | 202,909 | 63,349 | 53.3 |
| A4 at 25 % | 178,183 | 417,364 | 131,690 | 100.2 |

### Controls

| Control | Result |
|---|---|
| Shuffled candidate genes (seeds 1–5), union | 0 of 520 recovered for every seed |
| Literal masked-candidate attribution | 0 in every arm (by construction) |
| Unmasked model and candidates, union | 18 off-growth cases recovered (true reaction absent from the gapsmith model even without masking) |
| Solver replay, 27,027 fills | objective within 1e-6 in 26,962 (99.8 %); 65 fills in 6 genomes differ by up to 8.9 % (mostly one nitrogen-source test), GLPK higher; same added-reaction set in 10,999 (40.7 %) |

### Comparison with the pilot on the 22 Price genomes (same 180 cases)

| Test | Pilot (six-predictor union) | This run (DIAMOND baseline) | Both | Pilot only | This run only | p (exact McNemar, two-sided) |
|---|---:|---:|---:|---:|---:|---:|
| A1 | 10 | 4 | 2 | 8 | 2 | 0.11 |
| A2 | 16 | 5 | 4 | 12 | 1 | 0.0034 |
| A3 | 13 | 2 | 2 | 11 | 0 | 0.00098 |
| A1 + A2 + A3 | 23 | 7 | 5 | 18 | 2 | 0.00040 |
| A4 at 25 % | 12 | 3 | 2 | 10 | 1 | 0.012 |

The re-added reactions barely change (union 54 in the pilot, 43 here). What drops is gene attribution. The
DIAMOND baseline gives a candidate for the true enzyme to only 55 of the 146 off-growth Price-148 proteins, and
makes it the top candidate for 25 (`report/tables/candidate_rank.tsv`). For the Swiss-Prot cases the figures are
152 and 120 of 320.

## Checks

| Check | Result |
|---|---|
| C1 models | pass: 145/145, hashes match |
| C2 labels | pass: receipt with every query; Price-genome labels identical to the pilot's |
| C3 fence | pass: 0 entries with a panel taxid, 0 with a fenced genus |
| C4 leakage | pass: 0 hidden proteins with an identical or ≥ 95 % identical entry |
| C5 completion | pass: 49/49 genome jobs exit 0; all 520 cases in all seven arms, the unmasked run and the masked-candidate run |
| C6 masking audit | pass for every case (the runner stops on any failure) |
| C7 shuffled control | pass: 0 recovered for every seed |
| C8 solver replay | objective criterion pass (99.8 % ≥ 99 %); the design's same-reaction-set criterion fails (40.7 %), as in the pilot |

## What this shows

The pilot's main finding holds on the larger panel. Biomass gap-filling recovers almost none of the off-growth
enzymes (2 of 466). The targeted tests recover 34, in 20 genomes. Every pre-registered comparison against biomass
gap-filling is supported, and the shuffled-gene control gives zero. The rule-based quarter of the tests keeps 59 % of
the recoveries, short of the pre-registered 80 % (H3 not met, as in the pilot).

The effect stays small, and the calls are not usable as they stand. Recovery is 7 %, and most gene calls on labelled
proteins are wrong (54 % of distinct labelled decisions in the union). Four arms exceed the pre-registered harm gate,
which is measured against a biomass baseline whose own harm rate is lower than in the pilot. With the homology
candidates, recovery on the Price reannotations falls from 23 to 7 of 180 cases. The masked reactions are found
about as often as before, but the homology baseline rarely ranks the right gene first for these reannotated
functions. The gene assignment, not the simulation, is the limit.

The added cases also change what "off-growth" contains. Many recovered Swiss-Prot enzymes are biosynthetic rather
than catabolic, for example archaeal lipid synthesis (2.5.1.41, five genomes), homoserine O-acetyltransferase
(2.3.1.31, four genomes), enolase and RuBisCO. They are off the growth path only because the gapsmith biomass and
model do not need them, especially in archaea. Catabolic recoveries remain: mannose isomerase (5.3.1.7) in three
Price genomes, *N*-acetylneuraminate lyase (4.1.3.3), 2-dehydro-3-deoxygluconokinase (2.7.1.45) and alcohol
dehydrogenase (1.1.1.1).

## Limits

* Hidden cases are experimental labels only. 96 of the 145 genomes have none, and Swiss-Prot experimental labels
  concentrate on 17 model organisms (303 of 372 proteins).
* One candidate source (homology transfer from a genus-fenced Swiss-Prot) is used throughout. It is much sparser than
  the six-predictor union, so the Price-genome numbers here are not comparable to the pilot's except in the paired
  table above.
* The Biolog panel and the default medium are bacterial; for archaea the sole-source tests are a poor model of real
  growth, and 225 of the 466 off-growth cases are archaeal.
* Harm counts only labelled proteins; 98 % of calls are unconfirmed and are not counted either way.
* The linear programs have alternative optima: reaction sets agree in 41 % of replays, and 0.2 % of GLPK fills were
  not optimal (objective up to 8.9 % above HiGHS).

## Reproduce

Inputs on ws (`benchmark/v3/offgrowth/full_panel.py genomes | labels | reference | hidden-faa`), then on Shaheen
`benchmark/v3/shaheen/offgrowth-full-candidates.sbatch` (leakage check, final reference, candidate tables),
`offgrowth-genome.sbatch` (unchanged, `GENOMES=<49 ids> EXTRA_CASES_TSV=<labels> CANDIDATES_DIR=… OUT_ROOT=…
REPO_ROOT=… WALLTIME=12:00:00 … --submit`), `offgrowth-full-aggregate.sbatch` (aggregate plus
`full_panel_tables.py`) and `full_panel_replay.py` (C8). Jobs: candidates 15756966, genomes 15756986–15757034,
aggregate 15759626, replay check 15760157 (`slurm/sacct.txt`).
