# Intracellular products in metabolic experiments

The experiment engine now accepts explicit intracellular product demands in
addition to extracellular exchanges. This addresses a concrete limitation in
the [germicidin diagnostic](v3-germicidin-diagnostic.md): Sco-GEM represents
germicidins with intracellular drains, and growth alone misses their loss.
The implementation reports conditional **synthesis capacity**, never secretion,
expression, a recovered function or a calibrated probability.

`baseline.product_demands` declares exact existing reaction identifiers.
`product_demand` names one of them on an evidence record or masked condition.
The reporter must drain one unit of one intracellular metabolite, have no gene
association, and differ from biomass and maintenance. An experiment cannot
delete its own reporter. The model's chemistry, transport and finite bounds
remain unchanged; uptake through the drain stays closed. Forward and reversed
drain conventions are supported. The existing `product_requirement` applies:
growth-fraction by default, or explicitly maintenance-only.
Demand results expose the signed reporter bounds. A blocked reporter is a
model-bound limitation, and a positive lower bound can force production;
neither is independent biological evidence. Standard results also flag when
maximum capacity reaches the reporter cap. Bounds are reported rather than
silently relaxed. Intracellular here means outside the explicitly configured
extracellular compartments; it does not mean exclusively cytosolic.

The standard menu combines a demand with interventions that share its evidence
or hypothesis. The masked builder combines it with the candidate's supported
conditions; every selected condition compares every assignment, including the
null assignment, with and without the target. Product LPs count against the
same budget. Demand results remain in the existing product category, which is
excluded from learned enzyme features. The standard planner prompt is version
3 and carries typed synthesis results. The masked planner also receives the
product requirement before choosing a condition.

For masked demands, `product_contrast` distinguishes lost, reduced, retained or
absent control capacity. Valid maintenance-only comparisons can complete at zero
growth while retaining the separate `no_growing_control` growth result. An LP
failure remains inconclusive. When the required growth cannot be achieved,
product capacity is unassessed rather than reported as zero.
`growth_status` and `product_status` report endpoint validity separately; the
cell is `completed` if either contrast is informative. Consumers must inspect
the status of the endpoint they use. Product failure therefore does not erase
a valid lethal-growth result. Under `growth_fraction`, each genotype uses its
own growth optimum; the recorded absolute floors may differ, so a capacity
contrast can reflect both resource allocation and reaction deletion. This is
not a controlled comparison at one shared growth rate.

The diagnostic script is
[`check_intracellular_products.py`](../../benchmark/v3/check_intracellular_products.py).
It checks both evidence-guided production CLI routes on an explicitly synthetic
network without product transport. On the pinned, unmodified Sco-GEM model it
tests all four germicidin demands, phosphate closure, target deletion and
carbon closure. Independent COBRApy/GLPK calculations check the numerical
results. Live Qwen runs retain their prompts, replies, request usage and selected
conditions. These are integration diagnostics, not a comparison of planner
accuracy or efficiency.

Validation on Unimatrix CPU jobs 10219–10222:

| Check | Result |
|---|---|
| Final Python experiment/community tests | 159 passed |
| Affected benchmark/oracle tests | 121 passed |
| Focused Java builder tests | 54 passed |
| Full Java suites | 1,436 cases, zero failures/errors, six optional skips |
| Legacy fixed/random/evidence runs | Six before/after outputs exactly equal |
| Production CLI | Standard and masked evidence-guided routes passed |
| Synthetic CLI masked oracle | 24 phenotype comparisons and six demand classifications agree |
| Real standard demands | All 12 capacity comparisons agree |
| Real masked demands | All 72 phenotype comparisons and 36 demand classifications agree |
| Final live Qwen | Six valid selections; zero endpoint failures/fallbacks |

The final fifteen-method demand regression suite produces one failure and twelve
errors against the pre-extension code; all fifteen pass after the changes. The
first Java run had one incorrect test assertion: reversing evidence records
also reversed the private audit list, although the experiment menu was identical.
The corrected test separately checks menu-order invariance and independence
from hidden gene rules. The original failure is retained. The full Java gate
preceded the final Python-only review fixes; both CLI routes were rerun with the
validated Java artifact and final Python sources.

Qwen's three standard choices were growth deletion, baseline germicidin-A
capacity and phosphate closure. It used eight LP calls and did not reach a
product-deletion test. Its rationale still described feasible synthesis as
showing an active pathway, which is unsupported. The masked agent chose the
three conditions for product A, using 52 LP calls; that product alone cannot
distinguish the four-reaction assignment from the A-only assignment. Final usage
was 19,716 prompt tokens and 221 completion tokens across both agents. These
results demonstrate execution and feedback, not better experiment selection.

The [retained archive](../../benchmark/v3/results/intracellular-products-v1/record.json)
contains sources, original failures, generated CLI manifests/results, independent
solver results, prompt transcripts, usage and source hashes. Its
[`verify_archive.py`](../../benchmark/v3/results/intracellular-products-v1/verify_archive.py)
checks the inventory, pinned inputs, test counts, exact legacy parity, numerical
comparisons and transcript bindings without running a model or endpoint.

The real-model hypotheses and targets are controlled inputs. Germicidin
reactions in this SBML lack EC annotations, so this run does not demonstrate
automatic evidence-guided discovery of their assignments. The original gene
rules are known; masking them and proposing controlled alternatives does not
make this an independent function benchmark. The source model's medium, strain,
intracellular product endpoints and congener identities still differ from the
published production experiment. No assay labels are read during execution,
and no experimental accuracy is calculated. All twelve whole-item acceptance
requirements remain open.
