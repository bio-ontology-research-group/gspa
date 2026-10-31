# Evidence-guided metabolic experiments

`scripts/metabolic_experiments.py` runs a small, auditable set of in-silico
experiments (scipy/HiGHS LPs) on a real SBML-FBC model to test what supplied
evidence and competing hypotheses imply *under that model*. It sits next to
`scripts/stoichiometric_locus_refiner.py` (repair) and
`service/metabolic_conditions.py` (nutrient dependence) and does not replace
either, nor any scoring.

For competing **gene assignments with the original target GPR masked**, see
[masked assignment experiments](../benchmarks/v3-masked-assignment-experiments.md).
That format supports both explicit assignment clauses and automatic construction
from this run's raw/context candidates through `--experiment-mode evidence_guided`.
It compares every assignment, including null, under each selected condition. Its
conditional signatures do not become ordinary FBA confidence features. The
reaction-perturbation format documented below remains available separately.

For competing **absence explanations**, see [explicit absence comparisons](../benchmarks/v3-absence-comparison.md)
and [automatic comparisons from run gaps](../benchmarks/v3-guided-absence.md).
The guided route retains candidate/reaction alternatives, shares one solver and
condition budget, and preserves unavailable evidence. It does not classify
biological absence or promote annotations.

## What the output means

- Every result is **conditional** on the model, medium, bounds and maintenance
  constraint. A planner-chosen objective is never independent proof of gene
  function.
- Every record carries `calibrated_probability: null` and `promotion: "none"`.
  Nothing is scored and no gap-fill is promoted automatically.
- **Secretion checks are model properties, not production.** For a product
  exchange the record gives `max_secretion` at a growth floor (can this model
  secrete it?) and `min_secretion_at_growth_floor` (can it avoid secreting it
  while growing?). Classes: `no_secretion_possible_in_model_at_growth_floor`,
  `secretion_feasible_not_unavoidable`, `secretion_unavoidable_in_model`. A
  nonzero minimum only means the stoichiometry couples that exchange to growth;
  none of these are observed or expressed production. If either LP fails or is
  infeasible the class is `inconclusive` and never "no secretion possible".
- **Infeasible, failed and skipped runs are not function evidence.** Their
  `conditional_evidence.bearing` is `no_result`.
- **Optional vs. required flux.** For focus reactions (from the linked
  hypotheses) FVA at the growth floor gives `required_flux` / `optional_flux` /
  `carries_no_flux_at_growth_floor`, and a knockout re-solve gives
  `knockout_dispensable` (with `alternative_route_via` when a competing
  hypothesis' reaction takes over) / `knockout_reduces_growth` /
  `knockout_abolishes_growth` / `knockout_infeasible_with_maintenance`.
- Outcomes distinguish `growth_*`, `infeasible_with_maintenance` (no steady
  state meets the maintenance constraint) and `inconclusive` (solver error,
  non-finite or mass-balance-violating solution, budget exhausted mid-run).
  Skipped items are listed with a status (`skipped_budget`,
  `skipped_not_selected`). An empty menu gives
  `outcome: "no_candidate_experiments"` with zero solver calls.
- A non-growing or infeasible baseline does **not** stop condition selection.
  Alternative conditions can report `growth_recovered`; their growth ratio is
  undefined (`null`). A failed reference solve is unknown, not zero growth:
  successful conditions report `growth_without_baseline` in that case. Neither
  recovery nor lack of recovery alone establishes a candidate protein function.
- `config.product_requirement` explicitly selects `growth_fraction` (default)
  or `maintenance_only`. The latter tests product secretion even when biomass
  growth is zero, retaining the declared ATP maintenance and model bounds without
  an added biomass floor. Results record the requirement, actual biomass lower
  bound, and `min_secretion_at_required_state`. It can reveal feasible secondary
  metabolism in non-growing cells; it does not establish expression or observed
  production. The growth-based legacy field is omitted in maintenance-only mode.

Exchange signs are read from stoichiometry, so both `M_e ->` (v < 0 uptake) and
the refiner's reversed `-> M_e` (v > 0 uptake) are handled. Boundary species are
excluded from mass balance. Every solution is checked for finiteness,
`max|S v|` and bound violations.

**Objective.** `objectives.biomass` is preferred. Otherwise the SBML active
objective is used only if it is a single maximised term with a positive
coefficient; multi-term, minimised, negative or missing objectives are rejected
rather than guessed.

**Bounds.** Only infinite bounds are replaced (by `+-default_bound`, count
recorded); finite bounds are never clipped. NaN, `lb > ub` and `+inf` lower /
`-inf` upper bounds are rejected at load time.

**External boundary.** `baseline.exchanges` must list the trusted external
exchange reactions; each must be a single-metabolite reaction on an
extracellular species (`baseline.extracellular_compartments`, default `["e"]`).
Medium entries and the carbon/nitrogen/oxygen lists must be declared. Declared
exchanges not in the medium have uptake closed, and every *undeclared*
single-metabolite reaction (sinks, demands, biomass drains) has its uptake
direction closed, so nothing can be supplied from nowhere. Experiments can only
open or close declared exchanges.

## Manifest

```json
{
  "model": "model.xml",
  "config": {"strategy": "evidence", "seed": 0, "max_experiments": 8, "max_solver_calls": 60},
  "baseline": {
    "exchanges": ["EX_glc_e", "EX_ac_e", "EX_nh4_e", "EX_o2_e", "EX_vit_e"],
    "medium": {"EX_glc_e": 10, "EX_nh4_e": 10, "EX_o2_e": 10},
    "carbon_exchanges": ["EX_glc_e"], "nitrogen_exchanges": ["EX_nh4_e"],
    "oxygen_exchange": "EX_o2_e"
  },
  "objectives": {"biomass": "BIOMASS", "maintenance": {"reaction": "ATPM", "lower_bound": 1.0}},
  "context": {
    "hypotheses": [{"id": "H1", "gene": "b0001", "reactions": ["R1"], "group": "g",
                    "rescue_exchange": "EX_vit_e", "condition": {"carbon_exchange": "EX_ac_e"}}],
    "evidence": [{"id": "E1", "type": "transport", "hypothesis_id": "H1",
                  "substrate_exchange": "EX_ac_e", "role": "carbon"}]
  }
}
```

`maintenance` must be explicit (`"none"` is allowed but must be written).
Hypotheses sharing a `group` are competitors (declare it only when the candidates really are mutually
exclusive alternatives; the evidence-guided builder never does, because the same EC on several proteins can be
co-required subunits or isozymes). A hypothesis with `"perturbation": "reaction"` perturbs exactly its listed
`reactions`; a gene named beside them is an annotation and is never deleted, even when it is in the GPR. Genes resolve to reactions through
the SBML GPR (gene-product id or label).

The admissible menu is derived only from these records:

| Record | Experiment(s) |
|---|---|
| `transport` + `substrate_exchange` + `role` (carbon/nitrogen) | `carbon_source:<ex>` / `nitrogen_source:<ex>` (baseline sources swapped for the candidate) |
| `bgc` + `product_exchange` | `product:<ex>` |
| any record with `condition_hint: "oxygen_removal"` | `oxygen_removal` |
| `model_condition` + `condition_hint: "carbon_swap"`/`"nitrogen_swap"` + `exchange` | `carbon_source:<ex>` / `nitrogen_source:<ex>` (written by the evidence-guided builder from exact model identity; ignored when ordering by evidence) |
| `condition_hint: "supplement_removal"` + `exchange` | `supplement_removal:<ex>` |
| each hypothesis | `knockout:<id>` (optionally under its `condition`) |
| hypothesis with `rescue_exchange` | `ko_rescue:<id>:unrescued` / `:rescued` |

The **baseline** experiment also checks (FVA and a reaction knockout) the hypothesis reactions in the unchanged medium, up to
`config.max_baseline_focus_reactions` (default 6) and never at the expense of the planned experiments' minimum calls;
`measurements.focus_reactions_not_checked` lists the rest.

Items that do not resolve to the model go to `menu_rejections`. A hypothesis
may give candidate genes, model reactions, or both; reactions given directly
are used as-is and genes map through the model GPR. A candidate gene absent from
the GPR (for example when target assignments are masked) is listed in
`unmapped_candidate_genes` and does not remove a reaction-level hypothesis.
A gene-only hypothesis whose genes are unmapped has nothing to test and is
rejected.

## Selection strategies

| `strategy` | Behaviour |
|---|---|
| `fixed` | `config.panel` (list of menu ids) or the canonical kind order |
| `random` | seeded shuffle of the menu |
| `evidence` | deterministic: most distinct evidence types, then records, then competitors (an ordering, not a probability) |
| `qwen` | agent picks one id at a time given the menu and prior results; replies not in the menu are rejected and logged; after `qwen.max_invalid` invalid replies the run falls back to `evidence` order, logged |

Qwen only sees the menu and compact prior outcomes. It cannot edit the model or
run code. Endpoint conventions match `benchmark/agent/run_gspa_agent.py`: key
from `BORG_LLM_API_KEY` / `OPENAI_API_KEY` / `~/.config/borg-llm/api-key`
(never printed or written), base URL from `--api-base` or `GSPA_LLM_URL`.

## Budget

`max_solver_calls` is a hard cap across the baseline and all experiments
(inconsistent-bound checks cost no call). The policy is identical for all four
strategies and is written to `budget_policy` in every output; each menu item
lists its `min_calls` / `max_calls`:

- baseline: 1 LP, always first when the menu is non-empty;
- per experiment: 1 growth LP, plus 2 LPs (max/min secretion) for product items;
- per focus reaction: 3 LPs (FVA min, FVA max, knockout growth), at most
  `max_focus_reactions`, excluding reactions the item itself knocks out;
- an experiment is skipped (`skipped_budget`) if remaining calls < mandatory
  calls; focus checks are truncated to `floor((remaining - product LPs) / 3)`
  reactions and the shortfall is recorded;
- no calls are reserved for any planner, Qwen requests use no solver calls, and
  `max_experiments` counts every attempted item including budget-skipped ones.

## CLI

```bash
python3 scripts/metabolic_experiments.py --manifest frozen_manifest.json \
    --out experiments.json [--strategy fixed|random|evidence|qwen] [--seed N] \
    [--max-experiments N] [--max-solver-calls N] [--qwen-model M] [--api-base URL]
```

Output JSON: `hashes` (model SHA-256, config, context, manifest), `menu`,
`menu_rejections`, `selection_log`, `experiments` (each with mods,
`evidence_ids`, `hypothesis_ids`, measurements, status, outcome, solver-call
indices, `conditional_evidence`), `solver` (budget and per-call log),
`evidence` (used / not used), `assumptions`, `budget_policy`, `outcome`.

The top-level `conditional_evidence` list is the integration interface: one
record per experiment with `candidate_genes`, `unmapped_candidate_genes`,
`candidate_reactions`, `hypothesis_ids`, `evidence_ids`, `bearing`, `status`,
`outcome`, `model_sha256`, `config_sha256` and `solver_calls`. It carries no
confidence; connecting it to learned candidate attribution is left to the
integration step, which must supply its own calibration.

Python: `load_sbml_model(path)` then `run_experiments(manifest, model,
overrides=None, chat=None)`; pass `chat=callable(messages) -> str` to inject a
Qwen client.

## Production route

`gspa integrate --experiment-manifest manifest.json [--experiment-mode evidence_guided ...]` (or the
`experiment*` keys of the typed workflow config) runs this script as the `metabolic_experiments` stage under a hard
solver-call budget and wall-clock limit, with a fresh output, exit-status, schema, model-hash and "no probability, no
promotion" checks; a result that fails any of them fails the stage and is removed. Files go to
`<out>.experiments/` (`manifest.json`, `result.json`, `sidecar.log`).

* **`manifest` mode runs a user-supplied manifest as written.** The ledger says so
  (`context_source: user_manifest`) and makes no claim of evidence guidance. A hypothesis that names its
  `function_key` can become a candidate feature; others only appear in the raw record.
* **`evidence_guided` mode** treats the file as a template (the medium, biomass objective and explicit maintenance are
  its assumptions) and adds hypotheses and evidence from the run's candidates: an EC claim (retained or
  not) becomes a reaction-level hypothesis only through a validated EC to model-reaction mapping (refiner `reactions.tsv`,
  the model's EC annotation, or a reaction-graph id that is a model reaction id). A protein missing from the model's gene
  rules is still tested, as a reaction and never as a gene knockout; a gene label alone is not a mapping. Evidence is
  added where the run measured it (a raw claim, a context-positive window, a supporting cohort peer), where claim
  metadata names an exact model exchange id, or where `--experiment-evidence` supplies it; nothing is inferred from
  metabolite names. Features come from the per-reaction checks of a condition; a whole-gene knockout is never a flux class
  or a feature. Media are chosen from the template's `condition_menu` only where a hypothesis reaction handles the
  exchange's metabolite by the same species id or a shared compound cross-reference; the unchanged baseline is a
  matched control with its own checks. Proteins with the same EC share one hypothesis; a broad EC is reported as
  ambiguous, not as one reaction's activity. With nothing admissible the ledger records `no_admissible_experiment` and no solver
  runs.
* **Mapping to `ExternalContext.FbaConditional`**: conditions assessed, share with required flux, share with optional
  flux, mean relative growth change on knockout. A product-formation experiment, and a reaction the refiner added to
  make growth feasible, cannot confirm the candidate: those conditions are flagged, excluded from the evidence, and a
  candidate with only such conditions gets `objectiveRequiresCandidate = true`. The refiner builds its gene rules from the
  same EC evidence, so a gene in a rule is the candidate assignment, not an independent observation.
* `qwen` as strategy is forwarded as `--strategy qwen --qwen-model qwen3.8-27b --api-base http://unimatrix01.kaust.edu.sa:8000/v1`
  only; the sidecar uses its own key loader and the agent chooses only validated menu ids. The tests never call it.

Details and remaining gaps: [Learned, cohort and metabolic context](../concepts/context-workflow.md).

## Tests

```bash
python3 -m unittest scripts/tests/test_metabolic_experiments.py
```

Needs `libsbml`, `scipy`, `numpy`. The tests use a small stoichiometric network
to check mechanisms (carbon-source-dependent catabolism, auxotrophy rescue,
oxygen dependence, an unnecessary secondary product, alternative routes,
infeasibility, no-candidate, budget, selection, objective handling, bound
validation, closed internal sinks, masked-GPR hypotheses). They are not scientific
benchmarks.
