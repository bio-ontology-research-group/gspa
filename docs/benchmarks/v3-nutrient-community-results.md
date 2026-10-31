# Complete nutrient coculture results

All 1,320 cocultures and six frozen arms have completed: 7,920 case/arm plans,
39,600 independent LP comparisons and 15,840 GPR comparisons. The complete local
and remote collection audits agree: **39,599 numerical agreements and the one
original production solver unknown**, with no missing shards, additional
numerical discrepancies, comparable-endpoint disagreements or replay differences.
Every GPR comparison agrees. Numerical completion does not establish biological
cross-feeding, transport or enzyme attribution.

Shard 4 completed in 2:10:15; shard 5 also completed with exit zero. The original
shard-1 interrupted prefix and its independent replay remain preserved. The
original unknown remains `inconclusive` with null growth; the independent optimal
zero-growth solution does not replace it. See the [numerical audit](v3-nutrient-community-audit.md).

Job 10167 applied the [frozen dependence-aware analysis](v3-nutrient-community-analysis.md)
only after verifying the complete collection. It exited zero and retained every
arm and unknown. A separate standard-library checker reconstructs observations
from the frozen observation table and plans, then independently verifies all
7,920 endpoint records, 192 method summaries and 240 paired summaries across
48 groups. All eight completion output hashes match. This second check does not
repeat the resampling or re-estimate raw assay means.

## Primary assays

Each species contributes 600 primary cocultures per arm, across 40 environments
and six reused strains. The bulk endpoint is mean OD600 > 0.02; it does not show
that both organisms grew. The rescue proxy additionally requires both matched
monocultures to be below that threshold. Accuracy values below are percentages.
`printed` uses the recorded medium; `additional_inorganic` adds unprinted nutrients
as a sensitivity analysis. Oxygen values are model uptake-rate assumptions, not
measured rates. Neither arm nor oxygen is chosen after inspecting outcomes.

### Bulk-growth proxy

| Species | Medium | Oxygen | Permissive | Native unsupplied | No unsupplied uptake | Always growth | Never growth |
|---|---|---:|---:|---:|---:|---:|---:|
| Bs | additional_inorganic | 0 | 30.67 | 30.67 | 30.67 | 69.33 | 30.67 |
| Bs | additional_inorganic | 2 | 67.67 | 30.67 | 30.67 | 69.33 | 30.67 |
| Bs | additional_inorganic | 20 | 74.33 | 30.67 | 30.67 | 69.33 | 30.67 |
| Bs | printed | 0 | 30.67 | 30.67 | 30.72* | 69.33 | 30.67 |
| Bs | printed | 2 | 32.50 | 30.67 | 30.67 | 69.33 | 30.67 |
| Bs | printed | 20 | 32.50 | 30.67 | 30.67 | 69.33 | 30.67 |
| Ec | additional_inorganic | 0 | 52.83 | 47.83 | 47.83 | 46.17 | 53.83 |
| Ec | additional_inorganic | 2 | 51.17 | 47.50 | 47.50 | 46.17 | 53.83 |
| Ec | additional_inorganic | 20 | 51.17 | 47.50 | 47.50 | 46.17 | 53.83 |
| Ec | printed | 0 | 53.83 | 53.83 | 53.83 | 46.17 | 53.83 |
| Ec | printed | 2 | 53.83 | 53.83 | 53.83 | 46.17 | 53.83 |
| Ec | printed | 20 | 53.83 | 53.83 | 53.83 | 46.17 | 53.83 |

Bs = *B. subtilis*; Ec = *E. coli*. The starred value uses 599 scorable cases,
with the original unknown still counted among 600 total cases. Every other
primary method summary has 600 scorable cases. Paired comparisons always use
the same complete cases for both methods.

### Paired-rescue proxy

| Species | Medium | Oxygen | Permissive | Always rescue | Never rescue |
|---|---|---:|---:|---:|---:|
| Bs | additional_inorganic | 0 | 43.33 | 56.67 | 43.33 |
| Bs | additional_inorganic | 2 | 68.67 | 56.67 | 43.33 |
| Bs | additional_inorganic | 20 | 68.33 | 56.67 | 43.33 |
| Bs | printed | 0 | 43.33 | 56.67 | 43.33 |
| Bs | printed | 2 | 45.67 | 56.67 | 43.33 |
| Bs | printed | 20 | 45.67 | 56.67 | 43.33 |
| Ec | additional_inorganic | 0 | 62.50 | 21.83 | 78.17 |
| Ec | additional_inorganic | 2 | 61.50 | 21.83 | 78.17 |
| Ec | additional_inorganic | 20 | 61.50 | 21.83 | 78.17 |
| Ec | printed | 0 | 78.17 | 21.83 | 78.17 |
| Ec | printed | 2 | 78.17 | 21.83 | 78.17 |
| Ec | printed | 20 | 78.17 | 21.83 | 78.17 |

Only the permissive community has the matched monoculture executions needed for
the rescue comparison. The two restrictive policies are not assigned invented
rescue scores.

## Interpretation and sensitivity

Under the printed medium, Bs bulk accuracy is 30.67–32.50%, below always-growth
at 69.33%; its rescue accuracy is 43.33–45.67%, below always-rescue at 56.67%.
Ec community predictions are no growth throughout the printed-medium arms and exactly matches
the never-growth/never-rescue baselines. Thus the printed-medium results do not
establish a useful discrimination of these endpoints.

Adding unprinted inorganic nutrients changes the result substantially. For Bs
at oxygen 20, permissive bulk accuracy is 74.33%, five percentage points above
always-growth: 30 corrections and no newly introduced errors relative to that
constant baseline. Both descriptive resampling ranges for the difference are
[0, 12.5] percentage points, and all leave-one-strain-out differences are +5 points.
This is an outcome under a sensitivity medium, not a validated recipe correction
or an independently selected best model.

For Bs rescue at added-inorganic oxygen 2 and 20, differences from always-rescue
are +12.00 and +11.67 points. These include 166 corrections/94 new errors and
145 corrections/75 new errors, respectively. Environment-only ranges exclude
zero, but the shared-environment-and-strain ranges are [-11.68, 45.77] and
[-10.77, 45.58] points. The stronger apparent result therefore depends on the
assumed dependence structure; neither range is a calibrated confidence interval.
No significance, unseen-strain performance or model promotion is claimed.

Ec added-inorganic rescue accuracy falls to 61.50–62.50%, below never-rescue at
78.17%. Differences are -15.67 to -16.67 points, with 61–64 corrections but
155–164 new errors. All six leave-one-strain-out differences are negative;
shared-strain resampling ranges still include zero. All negative results and
all 240 registered paired contrasts remain in the report.

## No-added-nitrogen controls

Each species contributes 60 separate controls per arm (four environments).
Observed positive bulk endpoints occur in 13 Bs and 22 Ec controls; the rescue
proxy is positive in 13 Bs and 17 Ec controls. These measured bulk signals must
not be silently converted to no-growth labels.

Bs predicts no growth/rescue under every arm and policy: 47/60 correct for each
endpoint. Ec printed-medium arms also predict none: 38/60 correct for bulk and
43/60 for rescue. Under all three Ec added-inorganic oxygen arms, permissive
bulk predictions give 22 true positives and 23 false positives (37/60 correct);
permissive rescue gives nine true positives and nine false positives (43/60
correct). Each restrictive policy gives three false-positive bulk calls and
no true positives (35/60 correct). These four-environment controls support only
very limited descriptive sensitivity, not population inference.

## Remaining acceptance work and evidence

The growth model still uses predicted exchange ports, a fixed 1:1 biomass ratio,
and conditional balanced growth. The assay measures bulk OD. No measured exchanged
metabolite, transporter knockout, member-specific growth, recovered enzyme or
false enzyme-assignment effect is verified here. Those original requirements
remain open; no whole audit item is accepted.

The complete evidence archive is
`benchmark/v3/results/nutrient-community-complete-v1/`. It retains all six completed
shards (using the independent shard-1 replay), the interrupted history, frozen
inputs, complete numerical audit, merged predictions, all 48 analysis groups,
source and execution receipts, and the independent endpoint check. Earlier
incomplete checkpoints remain unchanged.
