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
- **Conditional feasibility vs. actual production.** For a product exchange the
  record gives `max_secretion` at a growth floor (feasible?) and
  `min_secretion_forced_by_growth` (actually required?). Classes:
  `not_producible_at_growth_floor`, `feasible_but_not_required`,
  `production_required_by_growth`.
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
  `skipped_not_selected`, `skipped_baseline_not_growing`). An empty menu gives
  `outcome: "no_candidate_experiments"` with zero solver calls.

Exchange signs are read from stoichiometry, so both `M_e ->` (v < 0 uptake) and
the refiner's reversed `-> M_e` (v > 0 uptake) are handled. Boundary species are
excluded from mass balance. Every solution is checked for finiteness,
`max|S v|` and bound violations.

## Manifest

```json
{
  "model": "model.xml",
  "config": {"strategy": "evidence", "seed": 0, "max_experiments": 8, "max_solver_calls": 60},
  "baseline": {
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

Items that do not resolve to the model go to `menu_rejections`.

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
(inconsistent-bound checks cost no call). An experiment is skipped when its
mandatory calls do not fit, and FVA/knockout focus checks (3 calls per reaction,
at most `max_focus_reactions`) are truncated to what remains.

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
`evidence` (used / not used), `outcome`.

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
infeasibility, no-candidate, budget, selection). They are not scientific
benchmarks.
