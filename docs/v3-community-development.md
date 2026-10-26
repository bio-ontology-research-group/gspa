# v3 community benchmark (item 9): measured Keio auxotroph pairs

`benchmark/v3/run_community_benchmark.py` runs the 1,035 measured Wintermute & Silver 2010 pairs
(`wintermute2010_cases.jsonl`, see `docs/benchmarks/v3-metabolic-data.md`) through the production
community engine (`scripts/community_metabolism.py`, which builds on `metabolic_experiments.py`).

**Status: implemented and mechanistically tested; no real run has been made.** The full run is the
root's Slurm job. Nothing below is a result.

## What is tested, and what is not

Tested: whether a *static balanced-growth* model of two gene-deletion strains sharing a medium pool
can grow together where neither grows alone, compared with the measured co-culture label. Not tested:
4-day fold change, lag, density, competitive fitness, exchanged metabolites, transporters, or gene /
function recovery. A predicted rescue assigns no enzyme and says nothing about whether the
transporter exists. Model: iML1515 as a **development transportability baseline** (all 46 strains
map by exact locus; the authors used iAF1260).

## Fixed protocol (hashed, set before any outcome)

The protocol lives in `PROTOCOL` in the runner; its SHA-256 is in every report and pinned by a test, so a
change is a reviewed edit, not drift.

| Item | Value |
|---|---|
| Medium | measured M9 components from the case inputs; glucose uptake **10 mmol/gDW/h** (a normalisation assumption, not inferred from 20 g/l); ammonium, phosphate, other printed inorganic ions, H2O, H+ at ceiling 1000 |
| Trace metals | the list in `inputs.medium.assumed_not_in_paper` at ceiling 1000 (assumed required, not measured). **Sensitivity arm: trace absent** |
| Oxygen | 0, 2 and 20 mmol/gDW/h, three equal arms; none selected by outcome |
| Organics | none beyond glucose |
| Maintenance, biomass | model-native `R_ATPM` lower bound (6.86, read from the model); the model's single active biomass reaction |
| Community | balanced growth, common rate mu, fractions 0.5/0.5, objective max mu |
| Pool | every native extracellular species (unique SBML species id) is a pool metabolite; every native external exchange is a port, support `predicted` (model-derived, not experiment-curated) |
| Removal | every pool metabolite may leave up to 1000 (waste accumulates in the medium); removal supplies nothing |
| Internal bounds | kept; non-port single-metabolite reactions get no uptake |
| Deletions | locus tag -> gene product -> `Model.blocked_by_genes` (GPR-aware: complexes block, isozymes do not); a locus not in the model is an error, never guessed |
| Growth call | mu >= 1e-6 |

Six arms: oxygen {0, 2, 20} x trace {present (primary), absent (sensitivity)}.

## Per pair, four matched LPs: two monocultures, community, control (same pool, supply, removal, deletions)

1. each deletion strain alone (fraction 1),
2. the community,
3. **control**: the same pair with uptake from the pool forbidden for every metabolite that is not a
   medium nutrient or H2O, H+, CO2, O2.

The control is recorded explicitly per pair (`control_comparison`: category, community and control mu,
number of closed uptake ports) and summarised per arm (`control_comparison.categories`, also among
registered positives). It only tests whether growth *needs* the transfer of the forbidden metabolites.
H2O, H+, CO2, O2 and the medium nutrients remain exchangeable, so **a rescue whose control still grows
(`growth_persists_in_control`) cannot be ascribed to the blocked transfer**, and a rescue that is lost is
conditional on the model, not evidence of an actual exchange. The registered prediction is unchanged; a
secondary `prediction_control_gated` demotes such positives to negative (unknown if the control did not
run) and is evaluated separately (`evaluation_control_gated`).

**Prediction.** `positive` iff the community grows and both monocultures do not. One or two viable
monocultures are flagged separately (`one_viable_monoculture`, `growth_without_partner`) and are not
positive. `no_community_growth` is negative. Infeasibility under the fixed maintenance counts as
no growth in the primary analysis (the model cannot satisfy maintenance, so growth is not predicted),
but a **strict** variant treats it as `unknown`; solver errors and an exhausted budget are always
`unknown`. All 1,035 pairs stay in every denominator, including failures and unmapped strains.

**Labels** (paper thresholds on mean day-4 fold growth): synergy > 50, no_synergy < 8, intermediate kept
separate and never forced to binary. They are joined after each arm's simulations finish;
`simulate_pair` receives `inputs` only (a test flips every label and checks that inputs, digests and all
simulation fields are unchanged).

## Cost control

Monocultures and communities are cached on a hash of the exact condition (model, deletions as blocked
reactions, supply, removal, fractions, threshold, maintenance); a hit returns what the LP would give.
`build_lp` reuses COO triplets cached as an attribute of each `Model` (no module-level reference, so the
cache is freed with the model; keyed on the identity of `Model.S`, which is never modified; tests compare
the cached LP with a cold build and with the legacy entry-by-entry construction, and check the model is
garbage-collected) and `build_organism` lets 46 deletion variants share one
parsed model; LP semantics are unchanged (the existing 26 community tests pass). Each arm has a
solver-call cap (default: planned need), and the report records calls used, cache hits and misses,
infeasible and error counts, build and solve seconds, peak RSS, CPU count, library versions and Slurm
variables. Cases that fall past a cap become `unknown`, not no growth, and the run is marked `partial_run`.

Planned LPs per arm: 46 monocultures + 2 x 1,035 = about 2,120; six arms about 12,700. One LP build takes
about 0.15 s on the real model (4,087 x 6,087 for a pair); solve time is not measured here.

## Running

```bash
PY=/tmp/gspa-v3-testenv/bin/python
D=/tmp/gspa-v3-metabolic-data
$PY benchmark/v3/run_community_benchmark.py --cases $D/derived/wintermute2010_cases.jsonl \
    --model $D/raw/iML1515.xml --out <new dir>                      # all six arms
$PY benchmark/v3/run_community_benchmark.py ... --arms o2_20_trace_present   # one arm per Slurm array task
$PY -m unittest benchmark/v3/tests/test_run_community_benchmark.py          # fixture LPs only, no network
```

`--arms` (at least one known, distinct arm), `--max-cases` and `--max-solver-calls` (both at least 1) mark
the run `partial_run`; invalid selections are rejected before anything is written.

**Execution versus coverage.** `status` is `complete` only for the default invocation when every case of
every arm ran all four LPs to a solver answer and no primary prediction is unknown. Otherwise it is
`incomplete_execution` (failed cases, LPs skipped by the budget or ending in a solver problem, unknown
predictions) or `partial_run` (a deliberate subset). `report.execution` and each arm give `n_executed`,
`n_failed` (exceptions, unmapped strains, digest or model-hash mismatch), `n_incomplete_simulation`
(an LP was skipped or inconclusive) and `execution_complete`; `coverage` gives how many cases carry a
prediction. Exit code: 0 complete; 2 clean but partial; 1 otherwise.

Outputs in `<out>`: `report.json` (protocol, input / model / engine / community-module / runner / protocol
hashes, git commit, resources, per-arm summaries, limitations), `protocol.json`, `model_audit.json`,
`summary_<arm>.json`, and `cases_<arm>.jsonl` (every pair: strains, loci, the four simulation outcomes with
condition hashes and cache flags, both predictions, label, error text for failures).

Per arm the summary reports: coverage (executed, failed, unmapped, with a prediction), growing
monocultures (strains and pair counts), community growth, infeasible and solver-problem counts, control
growth, and the confusion of the three labels against positive / negative / unknown, in strata
(all, both monocultures not growing, one viable, both viable, unknown state), plus precision, recall,
specificity, balanced accuracy and MCC for synergy vs no_synergy (point counts and proportions, recall and
specificity also counting unknown as a miss), and the observed synergy rate among predicted positives and
negatives. All six arms are reported side by side; none is a "best" arm.

**Uncertainty.** The 1,035 pairs are dyads of 46 shared strains, so pair-level binomial (Wilson) intervals
would treat dependent observations as independent and are not reported. `strain_resampling` (primary,
strict and control-gated predictions) gives, for each metric, a leave-one-strain-out range with the
vertex-jackknife SE and the three most influential strains, and a vertex-bootstrap percentile interval
(1,000 replicates, fixed seed, strains resampled, pair weight m_a x m_b, self-copy pairs dropped; settings
hashed as `analysis_sha256`). These are descriptive sensitivity summaries, **not independent confidence
intervals**; with 46 strains they are rough. They do not cover model-choice or assumption uncertainty.

## Review of the label extraction (done on the acquired file, no model involved)

1,035 unique pairs of 46 strains, each with 4 day-4 readings; every `inputs_sha256` re-verifies; labels
178 / 430 / 427 (synergy / no_synergy / intermediate) agree with the stated thresholds in all 1,035
rows (0 conflicts); the 17.2% synergy fraction matches the paper's 17%; fold growth runs 2.9 to 104
(median 9.2); no self pairs among the 1,035. Known data caveat carried from the preparation step: three
self-pair controls exceed 8-fold on the mean, which the paper does not report.

## Findings from checking the model (no LP solved)

* All 46 loci map to iML1515 gene products; the model hash equals the case hash; 331 external exchanges
  map one-to-one to extracellular species.
* **leuB blocks no reaction** in iML1515 under GPR evaluation (a redundant isozyme in the model). Pairs
  containing leuB are simulated as that strain being unaffected by the deletion; the measured strain is
  an auxotroph, so these pairs are expected model misses. They are in `n_strains_blocking_no_reaction`,
  not removed.
* Native `R_ATPM` is 6.86, not the 8.39 of the auxotroph manuscript protocol; the fixed rule here is the
  model-native value.

## Gaps and honest limitations

* **No real run.** Expected cost is small (about 12.7k LPs) but unmeasured.
* **iAF1260 not run.** The authors' model was not acquired (the task allowed no network). A side-by-side
  with an iAF1260 file needs its own cases (locus mapping and hash); the runner takes any model whose
  hash equals the case hash, so that is a new cases file, not a code change.
* Trace metals, oxygen, glucose uptake and free waste removal are assumptions. The three oxygen arms and
  the trace-absent arm quantify sensitivity; they do not validate the assumptions. Trace absent is likely
  infeasible for every pair in this model; if so, those arms report all-negative with infeasibility counts.
* Static balanced growth with equal fractions cannot represent dynamics, ratios or competition.
* Port support is `predicted`; no transporter is experimentally curated, and the unmapped / `closed`
  ports are the model's own.
* Intermediate pairs (41% of the data) are measured but not binary; they appear only in the confusion
  and in `intermediate_predicted_positive`.
* The control cannot forbid transfer of H2O, H+, CO2 or O2; the nutrient list is in each arm summary.
* Development set only: no held-out split, and the protocol was fixed without reading results, but the
  pair set was chosen by the experimenters' auxotroph design, so rates are not population estimates.
* Monk 2017, Harcombe 2014 data and measured exchanged metabolites or transporter knockouts remain
  unavailable (see `docs/benchmarks/v3-metabolic-data.md`).
