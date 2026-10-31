# Germicidin: growth versus secondary-metabolism capacity

The real-model diagnostic shows why a growth-only experiment can miss a
secondary-metabolism function. In Sco-GEM, deleting the germicidin synthase gene
`G_SCO7221` removes synthesis capacity for all four represented germicidins while
leaving predicted growth unchanged. Closing phosphate uptake eliminates growth
but retains germicidin synthesis capacity under ATP-maintenance constraints.
These are model properties, not measurements of secretion or recovered unknown
functions.

## Primary observations

Thanapipatsiri et al. (2015), *A Streptomyces coelicolor host for the heterologous
expression of Type III polyketide synthase genes*,
[DOI 10.1186/s12934-015-0335-0](https://doi.org/10.1186/s12934-015-0335-0),
PMID 26376792, PMCID PMC4573997, provides gene-deletion and expression contrasts.
The complete Table 1 is retained: **24 cells across six constructs and four
reported product quantities**, comprising **11 detected, five not detected and
eight not determined**. Neither “None” nor “ND” becomes a numerical zero.
The reported uncertainty is retained without inventing an SD/SEM interpretation.

Table 2 and the methods identify M1317 as an M1152 derivative lacking `sco7221`,
the `srs` operon and the `rppA` operon. Reintroducing `sco7221` under the constitutive
promoter restores germicidin production in this background. The paper separately
reports similar growth/sporulation of the deletion derivatives on SFM agar, with
data not shown. The quantified production assays instead use five-day GYM
cultures; those two observations must not be combined into a matched growth and
production measurement. [Primary study](https://link.springer.com/article/10.1186/s12934-015-0335-0)

Full article XML, its download receipt, all raw table values and curated culture
metadata are archived. The culture metadata is explicitly marked as manual
curation with checked source phrases and paragraph hashes. The article is
available under CC BY 4.0; authorship and the source are preserved.

## Model experiment and complete results

The unchanged Sco-GEM v1.2.1 model and bundled minimal medium come from the
[pinned source intake](v3-secondary-source-intake.md), commit
`63003f8b01a692dc2b73f2c9e0bbde871a3fe6f0`. No product observation, fitted
product rate or experimental production label is passed to the solver.

The protocol includes all **18** combinations: three conditions (bundled minimal
medium, phosphate uptake closed, glucose uptake closed), three executions
(native, `G_SCO7221` deleted, native repeated after deletion), and two product
requirements (95% of the genotype's own maximum growth, or maintenance only).
The repeated native execution checks state isolation and determinism; it is not
an independent complementation experiment. ATP maintenance remains the original
2.64 mmol/gDW/h; glucose uptake in the bundled medium is capped at 0.8 mmol/gDW/h.
All other medium entries and source finite bounds are preserved by the declared
closed-uptake policy. The plan and input hashes are written before LP execution.

| Condition and requirement | Native growth, 1/h | Native maximum A/B/C/D demand flux, mmol/gDW/h | Deleted-gene maximum |
|---|---:|---|---|
| Minimal, 95% growth | 0.074194 | 0.018467 / 0.018390 / 0.016354 / 0.021105 | All zero |
| Minimal, maintenance | 0.074194 | 0.368327 / 0.366829 / 0.325776 / 0.421682 | All zero |
| Phosphate closed, 95% growth | 0 | Not assessed: growth requirement unavailable | Not assessed |
| Phosphate closed, maintenance | 0 | 0.368327 / 0.366829 / 0.325776 / 0.421682 | All zero |
| Carbon closed, either requirement | Infeasible with maintenance | No result | No result |

Every feasible product minimum is zero: product flux is optional, not required.
Native repetition gives identical results. The four maxima are **separate
optimizations**, not a jointly attainable production mixture or total titer.
Growth in the deletion strain is unchanged within numerical tolerance.
Product demands are intracellular sinks; none is relabelled as an external
exchange. The standard agent product menu still requires an external exchange,
so this diagnostic does not yet provide an agent route for intracellular demand
targets.

Final job **10218** independently compares all 18 scenarios, **90 LP results**,
with COBRApy/GLPK using its own SBML parser, GPR deletion and bound construction.
All agree with HiGHS at relative tolerance 1e-6 and absolute tolerance 1e-7,
including solver statuses and missing intervals. All 18 valid native-ID outputs
are exactly unchanged by the identifier repair. Two additional label tests
demonstrate the repaired failure on this same network. The complete run uses
216 production/before/label LP calls plus 90 independent GLPK solves.

The regression checks require positive native capacity, preserved growth after
deletion, zero capacity after deletion, maintenance-supported capacity without
growth, and carbon-starvation infeasibility. These controls were added after
inspecting the initial diagnostic, with unchanged medium, bounds and solver
settings. They are developmental mechanism tests, not preregistered biological
prediction performance.

## What this does not establish

The model represents M145, not the complete M1152/M1317 backgrounds; regulatory
mutations and plasmid copy number are absent. The minimal medium is not GYM,
which contains complex yeast and malt extracts. Demand flux is not supernatant
production, and relative HPLC areas are not steady-state fluxes. Model A/B/C/D
names and formulas have not been reconciled with the publication's A/B7 congeners.
Flaviolin and the experimental total-production quantity have no matched endpoint
in this diagnostic. Their observations remain in the complete table, unscored.

The source GPR already assigns all four germicidin synthase reactions to
`G_SCO7221`. Removing that known assignment tests the execution mechanism; it
cannot validate discovering the function of an unknown protein. No empirical
accuracy, EC improvement or recovered-function count is reported. The next
functional benchmark must independently establish product identity, transport,
medium and target assignments, then compare masked competing hypotheses without
feeding measured products back as bounds.

Archive and verifier: `benchmark/v3/results/knockout-identity-v1/`. All complete
outputs, source packages, primary observations, reviews and intermediate failures
are retained. All twelve whole acceptance items remain open.

A subsequent [product-identity audit](v3-germicidin-identity.md) confirms a concrete
name/formula conflict: measured germicidin A is C11H16O3, while the model's
`M_germicidinA_c` is C10H14O3. The model's C-labelled product matches that formula,
but structural identity is not established. Germicidin B has two formula-compatible
model candidates. No products, assays or existing simulations are relabelled;
name-only joins cannot support a biological accuracy claim.
