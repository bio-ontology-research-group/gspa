# Paired coculture comparison and dependence sensitivity

The analysis implementation is tested and its specification is frozen for all
1,320 cocultures and 7,920 case/arm plans. Full simulation, independent numerical verification and the frozen analysis
have now completed. See the [complete results](v3-nutrient-community-results.md);
no whole audit requirement is accepted.

## What is compared

`benchmark/v3/analyze_nutrient_community.py` wraps the existing receipt-verified
endpoint evaluator. Each species, medium/oxygen arm, and primary/no-nitrogen
scope is reported separately. Bulk-growth comparisons include all three exchange
policies and constant-positive/constant-negative baselines. Paired rescue uses
the permissive community and its matched permissive monocultures; the other
community policies do not have matched monoculture executions and therefore do
not supply rescue comparisons. The rescue endpoint remains the bulk-OD proxy
defined in the [execution protocol](v3-nutrient-community.md).

Every comparison uses the same complete cases for its two methods. Reports
include unknown observations and predictions, scorable coverage, confusion
counts, corrections, newly introduced errors, and paired accuracy differences.
Missing measurements or failed solves cannot become negative outcomes. Duplicate
observation identities, duplicate matched monocultures, and coculture identities
that disagree with the frozen inputs fail validation.

## Dependence and interpretation

The 15 pairs within each environment reuse six selected strains. Neither wells
nor pairs are independent draws of organisms. The analysis uses 2,000 seeded
resamples and reports two **descriptive percentile ranges**:

1. Resample complete environments with replacement, keeping all pairs in each
   environment together. This description is conditional on the fixed strain
   panel.
2. Independently resample environments and strain identities. An observed pair
   receives the environment multiplicity multiplied by both strain
   multiplicities. The same strain weight applies in either pair position;
   same-strain pairs are never invented.

The same weights apply to both methods in each comparison. A resample with no
complete cases has no accuracy value; it is omitted from percentile calculation
and counted explicitly. Fewer than two environments produce no range. All six
leave-one-strain-out comparisons are also reported by dropping every pair that
contains the excluded strain, without refitting or changing thresholds.

Shared members induce dependence between different pairs, beyond repeated
observations of the same pair; see Aronow, Samii and Assenova (2015),
[Cluster-Robust Variance Estimation for Dyadic Data](https://doi.org/10.1093/pan/mpv018).
Multiplicative factor weights are motivated by Owen and Eckles (2012),
[Bootstrapping data arrays of arbitrary order](https://doi.org/10.1214/12-AOAS547).
Neither reference establishes nominal coverage for this particular shared-strain
scheme with six selected strains. These are not calibrated confidence intervals,
significance tests, or evidence of accuracy on unseen genotypes. The four
no-nitrogen environments are especially sparse. Related nutrient combinations
may also share effects beyond the environment groups. Technical-replicate
measurement error and uncertainty in the OD threshold are not propagated here.

## Freeze and verification

The analysis freeze was written at `2026-10-27T06:31:45.198283+00:00`. It binds
the analysis source, execution source, execution protocol, and outcome-free case
catalogue. It reads neither observations nor predictions. This is an analysis
checkpoint while execution is running, not a claim that the entire study was
prospectively preregistered: earlier single-organism outcomes and preliminary
community numerical checks were already examined. No coculture accuracy result
was used to choose these comparisons, and every registered arm remains reported.

The final execution plan is the unchanged `nutrient-community-v1/plan-v3` cache
directory. The analysis freeze is `nutrient-community-v1/analysis-plan-v1.json`.
Use `freeze --plan PATH --out FILE` to create a fresh specification;
`analyze --plan PATH --run PATH --observations FILE --analysis FILE --out PATH`
verifies that specification, validates observations, and requires the complete
merged execution receipt. Existing freeze and output paths are never overwritten.
The report binds the underlying endpoint evaluation by hash; promotion is `none`.

Seven focused tests cover shared environment weights, shared strain weights in
both pair positions, undefined resamples, paired denominators, correction/error
direction, unknowns, sparse designs, order invariance, identity rejection, and a
synthetic freeze/run/analyze path. One initially failing test exposed mutable
aliasing between the returned freeze and the in-memory specification. The freeze
now makes a deep copy; tampering cannot silently alter the validation target.
The initial failure log is retained alongside the successful rerun.
The full benchmark suite passes 652 tests with 13 optional skips. These checks
verify analysis mechanics and integrity, not empirical predictive validity.

The analysis source, tests, frozen specification, source-reference receipts and
test logs are archived in
`benchmark/v3/results/nutrient-community-analysis-v1/`. The earlier execution
checkpoint remains unchanged. Full numerical and empirical results belong in a
separate completed-run archive when available.

## Completed execution

Job 10167 waited for all original array 10150 tasks to terminate, then used the
completed shard-1 replay and all five other shards. It first performs the complete
independent collection audit. Endpoint analysis is permitted only if all 39,600
comparisons are present and the sole disagreement is exactly the previously
recorded solver unknown, with no endpoint, replay or identity discrepancies.
Any additional discrepancy stops execution before observation scoring. The
original unknown remains unchanged in the merged predictions.

The job uses the existing frozen analysis specification and NumPy 2.2.6/Python
3.10.19 on a Slurm compute node. All arms remain reported; no outcome-based
selection or model promotion was introduced. The job exited zero; the
[complete results](v3-nutrient-community-results.md) include independent endpoint
recomputation and retain the original numerical unknown.
Execution package: `/data/hohndor/gspa-v3/community-completion-v1/`.
