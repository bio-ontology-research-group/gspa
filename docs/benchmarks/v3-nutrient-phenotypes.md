# Matched nutrient and supplementation phenotypes

This benchmark tests the production metabolic experiment engine against the
known-genotype, measured endpoint phenotypes audited in
[the nutrient dataset](v3-nutrient-environments.md). It separates successful
execution and numerical agreement from biological agreement. The deletion
identities are supplied experimental inputs; the model's GPR assignments are
not independent gene-function truth.

## Frozen design

The catalogue contains all 1,144 single-organism assays: 1,040 primary assays
and 104 no-added-nitrogen controls. Each runs under six predeclared arms, giving
6,864 conditions. The primary cases include 480 unsupplemented auxotroph assays,
480 matched supplemented assays and 80 wild-type assays. Cocultures remain in
the source catalogue for a separate community evaluation.

`run_nutrient_phenotypes.py prepare` reads the outcome-free inputs and the
model preflight. It binds the input file, model files, implementation and every
execution plan by SHA-256. `run` cannot read an observation file. Only the
separate `evaluate` command joins measured OD600 to completed predictions.
No case is selected by measured growth, replicate availability or agreement
with a model.

The two models are iML1515 for the EcNR1/MG1655-derived experiment and iYO844 for
the *B. subtilis* 168 experiment. Their exact model hashes and the strain
difference are recorded in the dataset audit. Canonical gene identities are
used to evaluate complete GPR expressions. In particular, the E. coli leuB
deletion retains the model's alternative route; serA removes all three
GPR-dependent reactions. Neither case is forced into a single-reaction loss.

## Medium and assumptions

Every arm uses the printed MSgg ingredients and the case's named carbon,
nitrogen and amino-acid supplement. “No added carbon” and “no added nitrogen”
retain basal tryptophan and thiamine where the model represents them. Organic
nitrogen sources can also provide carbon. Ammonium sulfate supplies both its
nitrogen and sulfur components.

The six arms combine oxygen uptake ceilings of 0, 2 and 20 mmol/gDW/h with:

- `printed`: only the printed ingredients, plus water/proton balance and the
  stated oxygen assumption.
- `additional_inorganic`: the same medium plus sulfate, Fe(II), Cu(II), cobalt,
  molybdate and Ni(II), wherever the model represents their exact exchange.
  This is an explicitly unprinted supply sensitivity, including sulfur; it
  must not be described as the experimentally measured recipe.

All arms are reported. None is chosen by its observed agreement. Fixed uptake
ceilings are 10 for each named carbon/nitrogen source, 1 for an amino-acid
supplement (5 for serine), 0.1 for basal tryptophan, 0.01 for thiamine, 10 for
MOPS, and 1,000 for basal inorganic components and water/protons. These are
modeling assumptions in mmol/gDW/h, not rates measured in the paper or inferred
from the reported concentrations. If several ingredients supply one exchange,
the largest applicable ceiling is used; ingredient rates are not added.

Every missing exact exchange is recorded as an unrepresented ingredient.
Thiamine is not replaced with thiamine monophosphate. MOPS is supplied where
its exact exchange exists and remains explicitly unrepresented otherwise.
The frozen binding audit records MOPS as unrepresented in iML1515, and chloride
and thiamine as unrepresented in iYO844. These omissions are retained in every
applicable plan, including chloride supplied by the named glucosamine salt.
The engine closes uptake through every other single-metabolite reaction,
including undeclared sinks and demands. Exchange direction comes from
stoichiometry. Internal transport bounds and model-native ATP maintenance
(6.86 for iML1515; 9 for iYO844) remain intact.

## Production and independent execution

Auxotroph cases call `metabolic_experiments.run_experiments` with one explicit
genotype hypothesis, fixed selection, two LP calls, and no optional FVA. This
runs a matched model wild-type control and the genotype deletion in the same
medium. Those simulated controls are retained but are not counted as additional
independent measured wild-type assays. Actual wild-type assay cases call
`Runner.run(baseline_item())`, because the experiment selector intentionally
does no work for an empty candidate menu. No dummy candidate is invented.

The primary prediction is positive when the applicable optimal biomass flux
is at least 1e-6. An infeasible maintenance-constrained model is recorded as
such and predicts no growth; solver errors and unrun cases remain unknown.
The complete production records and solver-call ledger are retained for each
case. No annotation, confidence or learned feature is promoted by this driver.

`check_nutrient_cobra.py` independently parses the original SBML and GPRs with
COBRApy, reconstructs the bounds from the frozen manifest, deletes the canonical
genes, and optimizes with GLPK. It compares every WT/deletion flux and every
GPR-disabled reaction set to production SciPy/HiGHS results. It does not use
the production stoichiometric matrix, GPR evaluator, knockout list or solution.
Numerical agreement is an implementation check, not biological validation.
The complete planned comparison contains 13,200 LPs: one per WT assay and two
per deletion assay, across all six arms.

One registered prediction is unknown: `NUTRIENT2026:Bs:Glcn + Gly:R+R` under
`printed:oxygen=0`. Its argH knockout returned SciPy status 4 / HiGHS model
status 15 (unknown), whereas its matched wild type had a feasible zero-growth
solution. A separate numerical diagnosis reproduced the error with the default
and dual-simplex solvers; disabling presolve or using `highs-ipm` produced an
optimal zero-growth solution with mass-balance residual below 3.2e-10 and bound
violation below 5.7e-14. These diagnostic calls read no phenotype observations
and do not replace the frozen prediction. The registered two-call budget
contains no retry allowance. Unknown remains unknown in all phenotype scores.

The original independent replay completed all 13,200 comparisons but failed
its strict gate: besides the unknown production solve, two finite GLPK biomass
values differed from HiGHS by 3.76e-6 and 1.21e-6. Fresh GLPK models reproduced
the HiGHS values, identifying dependence on the reused solver basis. The revised
oracle resets GLPK to a standard basis before every solve, uses a primal
feasibility tolerance of 1e-9, and checks finite primal values, mass balance and
reaction bounds independently. Solutions with residuals or bound violations
above 1e-6 are unverified. The objective comparison tolerance remains 1e-6
absolute and relative. Production predictions and the frozen protocol remain
unchanged; both original failures and numerical diagnoses are preserved.

The revised full replay verifies all 13,200 disabled-reaction sets and 13,198
flux comparisons. Its two unverified comparisons are the retained production
solver-unknown case and a finite WT biomass difference of 1.31e-6 for
`NUTRIENT2026:Bs:GNAc + NoC:M+Met` under
`additional_inorganic:oxygen=20`. Both finite values predict growth. The strict oracle therefore exits 1, not 0; neither
the original run nor the final replay is reported as an all-comparisons pass.
Final numerical residuals and the complete discrepancy record are preserved in
`numerical-verification.json`.

The remaining finite comparison is tolerance-sensitive in GLPK: fresh solves
give 0.71949508714 at default tolerances, 0.71949639601 at primal/dual tolerances
of 1e-9, and 0.71949508713 at 1e-10. Four separate HiGHS configurations give
approximately 0.71949508713, matching production. This diagnostic explains why
small residuals alone cannot establish objective agreement. It does not replace
the full replay or retrospectively turn its failed comparison into a pass.
The full revised replay's maximum mass-balance residual is 6.13e-10 and maximum
bound violation is 5.88e-11. Both quantitative discrepancies and unknown outcomes
remain visible, even when the binary growth classification is unaffected.

## Scoring and limitations

Observed growth uses the raw-reconstructed mean OD600 > 0.02, retaining the
paper's monoculture threshold. This binary endpoint is compared with feasible
growth, not with predicted growth rate. The 96-hour assay, yield, lag, temperature
and regulation are not simulated by steady-state FBA.

Reports separate primary assays from no-nitrogen controls, both species, assay
kind, strain and environment. Paired supplementation rescue requires growth
with supplementation and no growth without it, in the same species, strain,
environment and arm. An unavailable member makes the paired result unknown.
The raw-derived experiment-level SEM is retained in the scored records.
Shared strains and environments mean the counts are not independent biological
replicates; no independent-sample confidence interval is claimed.

This benchmark does not discriminate missing annotation from missing sequence,
establish transporter function, validate a newly proposed enzyme assignment,
or evaluate coculture exchange. Those remain separate acceptance requirements.

## Measured results

The registered production run completed all 6,864 conditions in 2,215.8 seconds:
1,889 growth predictions, 4,974 no-growth predictions and one unknown. All
conditions, including all no-nitrogen controls, remain in the scored output.
The following table shows primary assays only. Each species has 520 individual
assays and 240 matched supplementation pairs per arm.

| Medium arm | Oxygen ceiling | B. subtilis assay agreement | E. coli assay agreement | B. subtilis rescue agreement | E. coli rescue agreement |
|---|---:|---:|---:|---:|---:|
| Printed | 0 | 49.7%* | 44.0% | 20.1%* | 20.0% |
| Printed | 2 | 63.5% | 44.0% | 45.8% | 20.0% |
| Printed | 20 | 64.0% | 44.0% | 46.3% | 20.0% |
| Additional inorganic | 0 | 49.8% | 76.3% | 20.0% | 52.5% |
| Additional inorganic | 2 | 81.9% | 77.5% | 67.1% | 53.3% |
| Additional inorganic | 20 | 84.8% | 77.1% | 69.2% | 52.5% |

*One unknown prediction gives evaluable denominators of 519 assays and 239
pairs. Its measured observation remains available. The unknown prediction is
neither counted as agreement nor as a negative prediction.

Under the printed medium, E. coli has no positive growth prediction in any
primary assay, at any oxygen ceiling. This includes 37 measured growing WT
conditions out of 40. B. subtilis with oxygen 2 or 20 predicts growth for only
3 of the 37 measured growing WT conditions under the printed medium. These
WT failures make a causal deletion interpretation unreliable in those arms.

The additional-inorganic arms improve individual-assay agreement, but the
paired rescue result remains weak. For each species, 192 of 240 measured pairs
show supplementation rescue. A constant positive prediction therefore achieves
80% accuracy on the complete observed pair catalogue. No model arm exceeds
that baseline. The archived `constant-baselines.json` also reports both constant
classes for individual assays, controls and every arm, without fitting a model
or changing any prediction. Full TP/FP/FN/TN counts and unknowns are retained in
`evaluation.json` by species, strain, environment, assay kind and control status.

These results do not establish accurate biological rescue inference. They expose
dependence on unmeasured supply/oxygen assumptions and disagreement with the
measured phenotype, despite extensive implementation tests. No sensitivity arm,
annotation, recovery method or learned model is promoted on these results.

## Reproduction

```bash
python benchmark/v3/run_nutrient_phenotypes.py prepare \
  --inputs "$DATA/derived-v3/inputs.jsonl" \
  --registration "$DATA/models/registered-preflight.json" --out "$PLAN"
python benchmark/v3/run_nutrient_phenotypes.py run \
  --plan "$PLAN" --out "$RUN" --ec-model "$EC_MODEL" --bs-model "$BS_MODEL"
python benchmark/v3/check_nutrient_cobra.py \
  --plan "$PLAN" --run "$RUN" --ec-model "$EC_MODEL" --bs-model "$BS_MODEL" \
  --out "$ORACLE"
python benchmark/v3/run_nutrient_phenotypes.py evaluate \
  --plan "$PLAN" --run "$RUN" --observations "$DATA/derived-v3/observations.jsonl" \
  --out "$EVALUATION"
```

Output directories/files must be new. Execution rejects changed plans, source
or models. Scoring checks complete prediction coverage and input/output hashes.
The Slurm wrapper stages immutable inputs onto job-local scratch, runs without
the observations, and copies outputs and exit status back before scratch is
removed. The current wrapper reserves two hours for production plus the revised
oracle and excludes node007 after its observed launch failure. The archived
original job retains its original wrapper; these scheduling changes do not alter
the frozen simulation. `TMPDIR` is the job scratch directory; local verification
uses `$HOME/tmp`.

Twelve focused tests cover outcome-invariant planning, exact chemical identity,
bounded supplementation and withdrawal, no-added-source media, missing
ingredients, unknown solver outcomes, canonical GPR deletion with an isozyme
bypass, genuine WT controls, reversed exchange signs, independent GLPK agreement,
matched pairing, absent observations, changed/incomplete receipts and oracle
independence from the order of prior media and gene deletions. The full
benchmark gate passes 638 tests with three optional skips. No production JVM
source changes are part of this benchmark slice.

## Archived evidence

The complete evidence package is in
`benchmark/v3/results/nutrient-phenotypes-v1/`: frozen inputs and plans, both SBML
models, every production result, original and revised oracle outputs, solver
diagnoses, failed/cancelled attempts, source snapshots, test logs and scoring.
`storage-index.json` binds original and stored paths and uncompressed SHA-256
hashes; large payloads use deterministic gzip. `SHA256SUMS` covers stored files.
The archive verification concerns file integrity, not biological acceptance.
