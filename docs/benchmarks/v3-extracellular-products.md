# Secretion experiments under candidate-supported conditions

The standard planner combined intracellular demands with candidate-supported
media changes and knockouts, but extracellular exchanges had only baseline
product tests. The masked Python engine calculated secretion intervals but
classified only intracellular product endpoints. A maintenance-only secretion
contrast could therefore remain inconclusive solely because neither genotype
grew. The Java evidence-guided builder also generated environmental combinations
only for intracellular demands.

The repair admits the same linked interventions for either product type, keeping
intracellular synthesis separate from extracellular net secretion. Masked cells
now report `growth_status`, `product_status` and the appropriate
`secretion_capacity_*` or `synthesis_capacity_*` contrast. Valid product evidence
survives absent growth, while failed or unassessed product LPs remain unknown.
The same complete-hypothesis comparison and solver-call budgets apply. No
simulation creates an enzyme probability or promotes a function assignment.

Python job **10419** completes with exit zero. Nine new mechanism tests fail
against the original engine (one assertion failure and eight missing-field/menu
errors); all nine pass after the repair. The complete affected experiment gate
passes **168 tests**, and **70** related benchmark tests pass. These check blocked
transport, reversed exchange signs against an analytically derived capacity,
mixed intracellular/extracellular targets, unassessed growth-dependent products,
solver failure, budget accounting and the observations supplied to the agent.
The agent callback tests are controlled interface checks, not live Qwen efficacy
experiments.

The same job runs the production Python engine on the pinned, unmodified Sco-GEM
v1.2.1 model and its bundled minimal medium. Both model exchanges, actinorhodin
and gamma-actinorhodin, are tested at baseline, after removal of phosphate or
glucose uptake, and after deletion of `G_SCO5089`. The registered masked panel
compares an explicit null assignment with controlled `R_ACTS1` and `R_GERMA`
assignments over both products and three media. Registration precedes simulation;
no measured production labels are read or fitted as flux bounds.

| Check | Result |
|---|---|
| Standard product experiments | All 8 agree with independent COBRApy/GLPK calculations |
| Masked control/knockout phenotypes | All 36 agree |
| Masked endpoint classifications | All 18 agree |
| Native maximum net secretion, maintenance only | Actinorhodin 0.113509; gamma-actinorhodin 0.116439 mmol/gDW/h |
| `G_SCO5089` deletion | Both product maxima become zero; growth remains 0.074194/h |
| Phosphate uptake removed | Growth becomes zero; both product maxima remain positive |
| Glucose uptake removed | Infeasible with maintenance; no product result |

These numerical results establish execution on a real network, not measured
product accuracy or recovery of an unknown function. The gene rule and candidate
alternatives are controlled inputs. The assay-medium and chemical-identity
limitations in the [source intake](v3-secondary-source-intake.md) remain relevant;
none of those observations is reinterpreted as benchmark truth.

Initial job **10416** reproduced the omissions, but three independent-solver
tests could not import COBRA because its cache defaulted to a read-only home
directory inside the container. Job 10419 explicitly places its cache in writable
job scratch. The original failed logs are retained.

Java job **10422** completes with exit zero. The new builder regression fails
against the original code; all **55** focused checks pass after the change. The
full JVM gate runs **1,474 cases**, with six optional skips and no failures or
errors. Four actual CLI runs compare the old and new Java artifacts using the
same repaired Python engine. The old masked builder offers no product/environment
combination; the new builder offers the supported carbon-source combination.
All **42** control/knockout endpoint comparisons and **nine** product classifications
across the two masked CLI runs agree with the independent checker. The standard
CLI menus, complete experiment records and annotation TSV bytes are exactly
unchanged between Java versions; masked candidate hypotheses are also unchanged.
These fixture runs test evidence-guided route integration, not biological accuracy.

Sonnet reviewed the Python change. The review prompted reversed-sign, mixed-target
and unassessed-product checks, including explicit target identity on an unassessed
standard secretion result. Its suggestion to reject all conditions touching the
product exchange was not adopted: supplying a metabolite is a valid declared
condition, and signed net exchange is distinct from synthesis. No gross secretion
or expression claim follows from these intervals. Product-only results also
remain separate from growth-derived conditional evidence; that growth evidence
must not be interpreted as disproving a secondary-metabolism function.

Evidence is retained in [the run archive](../../benchmark/v3/results/extracellular-products-v1/).
All twelve whole acceptance requirements remain open.
