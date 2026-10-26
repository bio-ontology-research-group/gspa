# Evidence-guided metabolic experiments

`scripts/metabolic_experiments.py` runs a small, auditable set of in-silico
experiments (scipy/HiGHS LPs) on a real SBML-FBC model to test what supplied
evidence and competing hypotheses imply *under that model*. It sits next to
`scripts/stoichiometric_locus_refiner.py` (repair) and
`service/metabolic_conditions.py` (nutrient dependence) and does not replace
either, nor any scoring.

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
Hypotheses sharing a `group` are competitors. Genes resolve to reactions through
the SBML GPR (gene-product id or label).

The admissible menu is derived only from these records:

| Record | Experiment(s) |
|---|---|
| `transport` + `substrate_exchange` + `role` (carbon/nitrogen) | `carbon_source:<ex>` / `nitrogen_source:<ex>` (baseline sources swapped for the candidate) |
| `bgc` + `product_exchange` | `product:<ex>` |
| any record with `condition_hint: "oxygen_removal"` | `oxygen_removal` |
| `condition_hint: "supplement_removal"` + `exchange` | `supplement_removal:<ex>` |
| each hypothesis | `knockout:<id>` (optionally under its `condition`) |
| hypothesis with `rescue_exchange` | `ko_rescue:<id>:unrescued` / `:rescued` |

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
