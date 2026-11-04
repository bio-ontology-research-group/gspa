# Report: `offgrowth2-measured-panel`

## 1. Cause of the defect

Tracing `--org-map-tsv` and `--matches-tsv` from `benchmark/v3/offgrowth2/score_driver.py` through `report.generate_report` and `report._score_genome_phase1_worker` into `score.build_genome_context` and `panels.load_org_id_map` / `panels.measured_panel` identified the root cause:

1. **Wrong genome-identifier column name in `panels.load_org_id_map`:** `load_org_id_map` read `row.get("genome_id")`, whereas `benchmark/v3/results/stage2b-followups-v1/price22_fitness_browser_organisms.tsv` stores the RefSeq/GenBank assembly accession under the column header `assembly` (`assembly, org_id, matched_proteins, assembly_proteins, fitness_browser_proteins, fraction_of_fitness_browser_proteome`). Because `row.get("genome_id")` was always `None`, `load_org_id_map` returned an empty dictionary `{}`.
2. **Unmapped assembly accession passed to `measured_panel`:** With `org_id_map` empty and `org_id` absent from `genome_results2.json.gz`, `report._score_genome_phase1_worker` fell back to the assembly accession (`GCA_...` / `GCF_...`) as `org_id`. Looking up assembly accessions against `experiment_panel_matches.tsv` (whose `org_id` column holds Fitness Browser identifiers such as `Smeli`, `Phaeo`, `BFirm`) matched no rows and returned `None` for every genome, causing every `N2` run to fall back to `transporter_panel` with `measured_panel_fallback = 1`.

## 2. Fix applied

- **`benchmark/v3/offgrowth2/panels.py`:**
  - Updated `load_org_id_map` to read the genome identifier from `assembly` (falling back to `genome_id` or `host_assembly`).
  - Updated `measured_panel` to validate TSV headers (`org_id` and `panel_id`), accepting both the 7-column `experiment_panel_matches.tsv` (`org_id, exp_name, exp_group, condition, exp_desc, media, panel_id`) and the 2-column Tier C table (`org_id, panel_id`, `benchmark/v3/offgrowth2/data/tierc_measured_panel.tsv`) via `--matches-tsv`, and to avoid restricting `allowed_ids` to an empty set when `panel` is empty.
- **`benchmark/v3/offgrowth2/score.py`:**
  - Updated `build_genome_context` to accept an optional `org_map_tsv` argument so callers invoking `build_genome_context` directly can resolve `org_id` from the organism map, and to treat an empty `measured_panel` list as `None`.
  - Updated `evaluate_arm_on_run` on `N2` to require a non-empty `ctx.measured_panel_ids`, setting `measured_panel_fallback = False` (`0`) whenever a measured panel is used and `True` (`1`) when falling back to `ctx.transporter_panel_ids`.
- **`benchmark/v3/offgrowth2/tests/test_panels.py` and `benchmark/v3/offgrowth2/tests/test_score.py`:**
  - Added unit tests covering `load_org_id_map` on the `assembly, org_id, ...` schema, `measured_panel` on both the 7-column Tier A and 2-column Tier C TSV formats, and end-to-end `build_genome_context` / `score_genome` / `generate_report` verification that a genome with rows in the matches table receives a non-empty measured panel distinct from its transporter panel (`measured_panel_fallback = 0`) while a genome without rows falls back to its transporter panel (`measured_panel_fallback = 1`).

## 3. Unit-test verification

Command:

```bash
PYTHONPATH=benchmark/v3 /home/leechuck/tmp/gspa-v3-testenv/bin/python -m unittest discover -s benchmark/v3/offgrowth2/tests -p 'test_*.py' -v
```

Output:

```text
test_clean_sidecar_conversion_and_helpers (test_candidates.TestCandidates) ... ok
test_score_bounds_enforced (test_candidates.TestCandidates) ... ok
test_baseline_fill_and_base_model_and_differential_t1 (test_engine.TestOffGrowth2Engine) ... ok
test_determinism_across_runs (test_engine.TestOffGrowth2Engine) ... ok
test_stability_intersection_across_equivalent_alternatives (test_engine.TestOffGrowth2Engine) ... ok
test_minimal_medium_filters_organics_and_preserves_bounds (test_medium.TestMinimalMedium) ... ok
test_sole_source_replacement_from_minimal_for_c_n_p_s (test_medium.TestMinimalMedium) ... ok
test_measured_panel_with_and_without_rows (test_panels.TestPanels) ... ok
test_random_panels_reproducible_by_seed (test_panels.TestPanels) ... ok
test_transporter_panel (test_panels.TestPanels) ... ok
test_report_consumes_runner_output_directly (test_report.TestReport) ... ok
test_report_end_to_end_files_and_markdown_tables (test_report.TestReport) ... ok
test_two_phase_byte_identical_to_single_process_and_resume (test_report.TestReport) ... ok
test_assign_gene_margin_floor_free_mutual_and_conaming (test_rules.TestRules) ... ok
test_constants_defined (test_rules.TestRules) ... ok
test_old_rule_and_budget_cut (test_rules.TestRules) ... ok
test_pathway_calls_adjacent_vs_distant_and_currency (test_rules.TestRules) ... ok
test_reaction_gate_fill_size_stability_draft_and_baseline (test_rules.TestRules) ... ok
test_masking_audit_failure_skips_case_and_records_it (test_runner.TestOffGrowth2Runner) ... ok
test_resume_skips_finished_cases (test_runner.TestOffGrowth2Runner) ... ok
test_runner_end_to_end_schema_and_determinism (test_runner.TestOffGrowth2Runner) ... ok
test_genome_seed_caches_reuse_across_calls (test_score.TestScore) ... ok
test_measured_panel_n2_differs_from_transporter_panel_and_flags_fallback (test_score.TestScore) ... ok
test_score_genome_recovery_decisions_q_controls_and_stats (test_score.TestScore) ... ok

----------------------------------------------------------------------
Ran 24 tests in 3.818s

OK
```

## 4. Full Tier A scoring re-run

`~/tmp/gspa-v3-offgrowth2/score-tierA-six/out/per-genome/` was deleted before launching the full re-run with `OFFGROWTH2_SCORE_JOBS=8`.

| Run | Wall time | User CPU (s) | System CPU (s) | Peak RSS (kB) | Exit status |
| --- | --- | ---: | ---: | ---: | ---: |
| Full Tier A re-run (`score-tierA-six/out`) | 7:03.84 | 1568.01 | 29.23 | 2194808 | 0 |

## 5. Per-genome measured-panel sizes and `measured_panel_fallback` counts

| Genome ID | Fitness Browser `org_id` | Matched experiment rows | Measured panel wells ($P_M$) | Transporter panel wells ($P_T$) | Unmasked admitted tests (`N2`) | Unmasked admitted tests (`N1`) | Total `N2` runs | `N2` runs with `measured_panel_fallback = 0` | `N2` runs with `measured_panel_fallback = 1` |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `GCA_000006965.1` | `Smeli` | 40 | 39 | 119 | 1641 | 1721 | 83 | 83 | 0 |
| `GCA_003443575.1` | `Phaeo` | 110 | 41 | 105 | 1467 | 1531 | 48 | 48 | 0 |
| `GCA_004402975.1` | `BFirm` | 122 | 49 | 124 | 1708 | 1783 | 2 | 2 | 0 |
| `GCF_000015245.1` | `SB2B` | 80 | 31 | 111 | 1262 | 1342 | 36 | 36 | 0 |
| `GCF_000016065.1` | `PV4` | 18 | 9 | 113 | 1353 | 1457 | 7 | 7 | 0 |
| `GCF_000018145.1` | `Dino` | 70 | 33 | 108 | 1490 | 1565 | 25 | 25 | 0 |
| `GCF_000020125.1` | `BFirm` | 122 | 49 | 124 | 1593 | 1668 | 55 | 55 | 0 |
| `GCF_000146165.2` | `MR1` | 56 | 27 | 113 | 1349 | 1435 | 41 | 41 | 0 |
| `GCF_000166295.1` | `Marino` | 83 | 30 | 108 | 1318 | 1396 | 17 | 17 | 0 |
| `GCF_000203935.1` | `ANA3` | 57 | 27 | 114 | 1388 | 1475 | 30 | 30 | 0 |
| `GCF_000236665.1` | `PS` | 20 | 4 | 108 | 1091 | 1195 | 5 | 5 | 0 |
| `GCF_000237365.1` | `azobra` | 54 | 28 | 117 | 1449 | 1538 | 8 | 8 | 0 |
| `GCF_000327065.1` | `psRCH2` | 107 | 43 | 118 | 1478 | 1553 | 45 | 45 | 0 |
| `GCF_000698265.1` | `WCS417` | 96 | 49 | 119 | 1512 | 1582 | 69 | 69 | 0 |
| `GCF_000832305.1` | `Cup4G11` | 61 | 35 | 120 | 1837 | 1922 | 32 | 32 | 0 |
| `GCF_001307155.1` | `pseudo3_N2E3` | 102 | 64 | 123 | 1487 | 1546 | 70 | 70 | 0 |
| `GCF_001307275.1` | `pseudo5_N2C3_1` | 92 | 55 | 121 | 1485 | 1551 | 70 | 70 | 0 |
| `GCF_001440325.1` | `pseudo13_GW456_L13` | 68 | 49 | 121 | 1561 | 1633 | 60 | 60 | 0 |
| `GCF_001623525.1` | `pseudo6_N2E2` | 110 | 52 | 121 | 1495 | 1564 | 60 | 60 | 0 |
| `GCF_001625455.1` | `pseudo1_N1B4` | 72 | 52 | 121 | 1490 | 1559 | 30 | 30 | 0 |
| `GCF_001633105.1` | `acidovorax_3H11` | 50 | 32 | 111 | 1445 | 1524 | 28 | 28 | 0 |
| `GCF_002797435.1` | `Korea` | 27 | 12 | 115 | 1322 | 1425 | 24 | 24 | 0 |
| **Total** | — | **1617** | — | — | — | — | **845** | **845** | **0** |

## 6. Open issues

None.
