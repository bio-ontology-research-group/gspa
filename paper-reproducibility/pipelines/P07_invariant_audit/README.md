# P07: 108-genome stoichiometric invariant audit (Fig. 6a, Suppl. S12)

Compares four ways of turning the DeepProZyme (DPZ) tau = 0.5 draft into a network: the raw
threshold draft, METEOR v8 (negative-cost MILP), the Grow-Support threshold repair (positive-cost LP;
the PSB-revision Section 3.2 comparator) and GSPA's StoichiometricLocusRefiner (positive-cost LP plus a
pFBA support certificate). The entry point is `run_108_genome_invariant_audit()` in
`benchmark/unified/pillar2_fba_lro_price149/run_pillar2_all.py`, and `run.sh` drives both the ibex
and the local steps.

## Archive provenance

| Item | Source (ibex, read-only) | How it enters the audit |
| --- | --- | --- |
| Cohort | `psb_revision/code_snapshot/panel108_gram.tsv` (108 GCA/GCF, Gram) | array index |
| Baseline variant | **`dpz_vanilla`** (DeepProZyme, no homology filter). Table 1 `baseline_dpz`/`meteor_dpz`, `meteor108_solutions_compact.pkl.gz` and `thresh_gapfill_deadends.py` all use it | |
| tau = 0.5 draft | DPZ-vanilla prediction pkl resolved by `code_snapshot/eval/baseline_io.resolve_baseline_pkl` | `draft_idx`: reactions with any mapped EC scoring >= 0.5 (same code as `thresh_gapfill_deadends.py`) |
| METEOR solution | `meteor_v8_evw_p2mu3_run/meteor_out/dpz_vanilla/meteor_sol_<GCA>.pkl` (`y_vals`, `v_vals`, 47,880 reactions) | `meteor_active_idx` = y > 0.5; `meteor_v_idx/val` = MILP flux |
| Grow-Support repair | no per-genome network was archived; re-derived with `meteor_v8.repair.grow_support` (same LP, PuLP 3.3.1 / CBC) | `gs_repaired_idx`, `gs_added_idx`, `gs_support_idx`, `gs_v_val` |
| MEMOTE profiles | `psb_revision/results/thresh_gapfill_deadends/tgd_<GCA>.json` | stored as `tgd`; dead ends and mass-imbalance counts of the Grow-Support network |

Checks: the re-derived draft, core and Grow-Support additions reproduce `tgd_<GCA>.json`
(`n_draft`, `n_core`, `n_added_by_repair`) for **108/108** genomes. METEOR active sets and flux-carrying
sets are identical to the earlier `meteor108_solutions_compact.pkl.gz` for 108/108 genomes. The ibex
`code_snapshot/data` files (`fba_matrices_v6RE.pkl`, tight bounds, `medium.pkl`, `seedr2ec.pkl`) have the
same sha256 as the local copies under `benchmark/manifests/reference_data/psb_revision/code_snapshot/data`,
and `allrxns` in the archive equals the refiner's reaction order.

Other files in the METEOR tree (inspected but not used):
- `meteor_df_<GCA>.pkl` holds the per-protein EC posterior (proteins x 7,433 ECs, float64, 100-330 MB).
  It has no reaction-level content and does not give the draft.
- `psb_revision/meteor_out/repair_rerun/<baseline_variant>/` contains METEOR verify-and-repair reruns
  (Task 6, 11 configurations, no `dpz_vanilla`), which are not the Grow-Support repair.

Jobs (ibex): 53986089 (venv with the archived versions numpy 1.23.4, pandas 1.5.3, scipy 1.13.1, cobra 0.31.1,
pulp 3.3.1, plus a schema probe), 53988457 (108-task extraction array, about 40 s per genome),
53988464 (merge and sha256 manifest). Write area: `/ibex/scratch/projects/c2014/hohndor/meteor108_repro_2026-10-19/`.

## Archive schema: `benchmark/manifests/reference_data/meteor108_drafts_compact.pkl.gz` (4.4 MB)

Top level:
- `allrxns`: 47,880 reaction ids. Every `*_idx` array indexes this list.
- `media[gram]`: `lb`, `ub`, `biomass_id`, `biomass_idx`, `excluded_idx` and `media_exchange_idx`, with
  gram = negative or positive.
- `genomes[GCA]`, `panel`, `missing` (empty), `S_sha256`, `description`.

Per genome:
- `draft_idx`, `draft_maxscore`, `core_idx`
- `meteor_active_idx`, `meteor_v_idx`, `meteor_v_val`, `meteor_biomass_flux`, `meteor_status`
- `gs_repaired_idx`, `gs_added_idx`, `gs_support_idx`, `gs_v_val`, `gs_*_growth_milp_medium`
- `tgd` (full T9 JSON), `tgd_match`, `pred_path`, `n_proteins`

## Definitions

- **added** = internal reactions of the network (not `EX_`, not `biomass_*`) that are not in the
  tau = 0.5 draft. The archived Grow-Support `n_added_by_repair` (81.3) also counts the exchange
  reactions the repair opens; the internal-only count is 52.7.
- **flux-certified added** = added reactions with |v| > 1e-6 in a pFBA solution: minimise sum |v|
  subject to S v = 0 and v_bio >= gamma_min = 0.1, with the network restricted to its own reactions.
- **zero-flux padding** = added minus flux-certified added.

Media:
- **METEOR's own medium** (`media[gram]`): the FVA-tight bounds plus `apply_media(["default"])`.
  `fba_matrices_v6RE` reverses the sign of exchange columns, so v > 0 is uptake. Every non-medium
  exchange therefore has (lb, ub) = (0, 1000): uptake is open and secretion is closed. Medium exchanges
  have (-100, 100). Only the exchanges inside the network carry flux.
- **Audit medium**: `StoichiometricLocusRefiner.get_media_bounds(gram, "minimal")`. This is the named
  minimal medium plus the 36-compound ModelSEED supplement. Uptake is allowed only for those compounds,
  and every exchange may secrete; the same bounds are used for Stage 2B.
  - `lp_biomass_feasible_pct` is measured on this medium as given (strict).
  - Flux certification on this medium adds the greedy-minimal set of the network's non-medium uptakes
    needed to reach gamma_min (`*_required_non_medium_uptakes`), because no METEOR or Grow-Support
    network grows on the strict medium.

## Caveats

- **Shared draft core.** GSPA's draft core is seeded from METEOR's archived zero-flux active set
  (p = 0.75, see `_solve_single_genome_slr`), not from the tau = 0.5 draft now in the archive. The GSPA
  row therefore isolates the repair policy rather than an independent reconstruction. The GSPA row is
  unchanged by this pipeline.
- **Grow-Support certification.** All Grow-Support additions carry flux in Grow-Support's own LP by
  construction (support |v| > 1e-9). The pFBA columns test whether they are needed at minimum total flux.
- **Mixed column.** `mean_imbalanced_added_or_active` holds MEMOTE mass-imbalance counts over the
  active network for the raw threshold, METEOR and Grow-Support rows (from Table 1 and the T9 JSON),
  and the count over added reactions for GSPA. `mean_imbalanced_added` uses one definition for all rows.
- **Raw threshold feasibility.** For the raw threshold, `lp_biomass_feasible_pct` is the archived
  Table 1 `fba_growth > 0`. `lp_biomass_feasible_pct_own_medium` is the T9 core-draft growth >= 0.1.
