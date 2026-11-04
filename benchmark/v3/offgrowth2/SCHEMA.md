# Output Schema (`genome_results2.json.gz`)

```json
{
 "schema": "gspa.offgrowth2.genome/1",
 "genome_id": "str",
 "candidate_table_kind": "six|clean",
 "inputs": {"path": "sha256"},
 "panel": [{"panel_id": "str", "element": "str", "seed_id": "str", "exchange_rxn": "str"}],
 "minimal_medium": {"exchange_id": 0.0},
 "default_medium": {"exchange_id": 0.0},
 "model_genes": ["str"],
 "draft_reactions": ["str"],
 "known_labels": {"protein_id": ["ec"]},
 "skipped_cases": [
   {
     "case": {"case_id": "str", "protein_id": "str", "truth_ecs": ["str"], "case_source": "str", "hedged": 0},
     "reason": "masking_audit_failed",
     "audit": {"case_id": "str", "protein_id": "str", "truth_ecs": ["str"], "equivalent_ecs": ["str"], "target_reactions": ["str"], "passed": false, "violations": ["str"]}
   }
 ],
 "runs": [
   {
    "run_kind": "masked|unmasked",
    "case": {"case_id": "str", "protein_id": "str", "truth_ecs": ["str"], "case_source": "str", "hedged": 0} | null,
    "audit": {"case_id": "str", "protein_id": "str", "truth_ecs": ["str"], "equivalent_ecs": ["str"], "target_reactions": ["str"], "passed": true, "violations": []} | null,
    "stratum": {"stratum": "off_growth|on_growth", "wt_biomass": 0.0, "ko_biomass": 0.0, "ko_fraction": 0.0} | null,
    "target_reactions": ["str"],
    "draft_reactions": ["str"],
    "baseline_fill": {
      "fill": ["str"],
      "replay_fill": ["str"] | null,
      "stable_fill": ["str"],
      "n_lps": 0,
      "cpu_seconds": 0.0,
      "grows_without_fill": false
    },
    "tests": [
      {
        "test_id": "str",
        "test_type": "T1|T2|T3|A1old|A0",
        "target_id": "str",
        "panel_id": "str" | null,
        "element": "str" | null,
        "gapfill_ran": false,
        "feasible_in_universal": false,
        "already_satisfied": false,
        "fill": ["str"],
        "replay_fill": ["str"] | null,
        "stable_fill": ["str"],
        "n_lps": 0,
        "cpu_seconds": 0.0
      }
    ],
    "n_lps": 0,
    "cpu_seconds": 0.0,
    "wall_seconds": 0.0,
    "peak_rss_kb": 0
   }
 ]
}
```

## Notes and documented additions

- All reaction identifiers (`draft_reactions`, `target_reactions`, `fill`, `replay_fill`, `stable_fill`) are normalized base ModelSEED reaction IDs (`base_seed_rxn_id`, e.g. `rxn00001`).
- All lists of reaction IDs, gene IDs, and EC numbers are deterministically sorted.
- `draft_reactions` is recorded both at the top level (the unmasked draft model's gene-evidenced base reaction IDs) and inside each run dictionary in `runs` (that run's draft model base reaction IDs, i.e. after case masking for `"masked"` runs).
- `skipped_cases` (top-level list, empty when all cases pass the masking audit) records any case whose mechanical masking audit failed (`MaskingAuditError`) so that failed-audit cases are never executed and never dropped silently.
- Test IDs follow the deterministic conventions:
  - `T1`: `T1_<panel_id>` (`test_type = "T1"`, `target_id = <biomass_id>`, `panel_id = <panel_id>`, `element = "C"|"N"|"P"|"S"`)
  - `T2`: `T2_unblock_<rxn_id>` (`test_type = "T2"`, `target_id = <rxn_id>`, `panel_id = null`, `element = null`)
  - `T3`: `T3_deadend_<met_id>_<unproduced|unconsumed>` (`test_type = "T3"`, `target_id = "<met_id>:<unproduced|unconsumed>"`, `panel_id = null`, `element = null`)
  - `A1old`: `A1old_<panel_id>` (`test_type = "A1old"`, `target_id = <biomass_id>`, `panel_id = <panel_id>`, `element = "C"|"N"|"P"|"S"`)
  - `A0`: `A0_biomass_default` (`test_type = "A0"`, `target_id = <biomass_id>`, `panel_id = null`, `element = null`)
