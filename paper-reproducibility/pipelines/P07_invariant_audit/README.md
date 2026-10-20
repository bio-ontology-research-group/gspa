# P07: 108-genome stoichiometric invariant audit (Fig. 6a, Suppl. S12)

Compares four ways of turning the DeepProZyme (DPZ) tau = 0.5 draft into a network: the raw
threshold draft, METEOR v8 (negative-cost MILP, archived solutions), the Grow-Support threshold repair
(positive-cost LP; the PSB-revision Section 3.2 comparator) and GSPA's StoichiometricLocusRefiner
(positive-cost LP plus a pFBA support certificate). The entry point is `run_108_genome_invariant_audit()`
in `benchmark/unified/pillar2_fba_lro_price149/run_pillar2_all.py`, and `run.sh` drives the ibex and the
local steps.

Every computation is GSPA code. The comparator (threshold draft, structurally excluded reactions,
METEOR's default medium, the Grow-Support LP, the maximal-biomass LP and the structural network profile)
is implemented in `scripts/meteor_threshold_repair.py` (numpy, scipy with HiGHS, pandas) with unit tests in
`scripts/tests/test_meteor_threshold_repair.py`. The archived METEOR tree is read as data only: the
per-protein EC score matrices, the METEOR MILP solutions and the archived counts used for comparison.
Until 2026-10-20 the ibex extraction imported METEOR's `code_snapshot/src` and `eval` modules; that
dependency is gone.

## Inputs (ibex, read only)

| Item | Source | How it enters the audit |
| --- | --- | --- |
| Cohort | `panel108_gram.tsv` in this folder (108 GCA/GCF, Gram; copy of `psb_revision/code_snapshot/panel108_gram.tsv`) | array index |
| Universal model | `fba_matrices_v6RE.pkl` (S, lb, ub), `tight_bounds_v6_{neg,pos}.pkl`, `medium.pkl`, `seedr2ec.pkl`, `all_ancestors.txt`, `enzymeobsolete.pkl` (plain pickles and text, copies under `benchmark/manifests/reference_data/psb_revision/code_snapshot/data`) and `universal_metadata.json` (reaction and metabolite tables written by GSPA) | `UniversalData` |
| Baseline variant | **`dpz_vanilla`** (DeepProZyme, no homology filter) | |
| tau = 0.5 draft | DPZ-vanilla prediction pkl: `meteor_v7_run/reinfer/dpz_vanilla/<GCA>_DPZ.pkl` when present, else `paperA_2026/baseline_preds/dpz_vanilla_genome_collection/<GCA>_DPZ.pkl` | `draft_idx`: reactions with any mapped EC scoring >= 0.5 |
| METEOR solution | `meteor_v8_evw_p2mu3_run/meteor_out/dpz_vanilla/meteor_sol_<GCA>.pkl` (`y_vals`, `v_vals`, 47,880 reactions) | `meteor_active_idx` = y > 0.5; `meteor_v_idx/val` = MILP flux |
| Grow-Support repair | re-derived: `grow_support(avail = not excluded, gamma_min = 0.1, core = core)`, HiGHS | `gs_repaired_idx`, `gs_added_idx`, `gs_support_idx`, `gs_v_val`, `gs_lp_objective` |
| Archived counts | `psb_revision/results/thresh_gapfill_deadends/tgd_<GCA>.json` | stored as `tgd` for comparison only |

The METEOR root is `/ibex/scratch/projects/c2014/kexin/funcarve` (kept in place); the work area is
`/encrypted/e3008/gspa/meteor108_repro_2026-10-20` (`scripts/`, `data/`, `per_genome/`, `logs/`).

## Agreement with the archived values

Two comparisons, both written by this pipeline:

1. `ibex/extract_meteor108_drafts.py` recomputes draft, core and Grow-Support repair per genome and
   stores the archived `tgd_<GCA>.json` counts next to them. Result (jobs 54048078, 54048495): draft and core
   sizes equal the archived values for **108/108** genomes; the number of reactions added by the repair
   equals the archived value for 48/108 and differs by at most five otherwise.
2. `smoke_test_reimpl.py` (`run.sh smoke`) compares the implementation with the archived
   `meteor108_drafts_compact.pkl.gz` on the local copies of the data files
   (`results/smoke_test_reimpl.{tsv,summary.json,media.json}`). Over 108 genomes: medium bounds, excluded
   reactions and open exchanges identical for both Grams; Grow-Support LP objective identical for 108/108
   (relative 1e-4) and core growth identical for 108/108; network sizes and mass-unbalanced reaction
   counts of the repaired, draft and METEOR networks identical for 108/108; FBA growth of the profile
   protocol within 0.01 for 108/108 (draft), 99/108 (repaired) and 105/108 (METEOR).

Two quantities are comparable rather than identical, for reasons in the method and the archive:

- **Grow-Support support.** The LP `min sum_{j not in core} |v_j|` has several optimal supports in most
  genomes. PuLP/CBC (archive) and HiGHS (GSPA) return different ones with the same objective value, so the
  added sets overlap (mean Jaccard 0.66) and the added count differs by at most five (mean -0.7). The
  maximal biomass of the repaired network follows the support and differs accordingly.
- **Dead ends.** The archived MEMOTE profiles were computed on a cobra sub-model built from the universal
  model's own reaction and metabolite objects, so every metabolite's reaction list still spanned the whole
  universal model and the count depended on reactions outside the network. GSPA counts the dead-end
  metabolites of the network itself (every reaction of the metabolite in the network is irreversible, not
  an exchange, and produces it), which gives 9.6 for METEOR and 16.6 for the raw draft per genome instead of
  the archived 2.4 and 9.8. The smoke test records a replica of the shared-object count
  (`*_deadends_shared_objects`) only to document the discrepancy.

## Uniform structural profile in the audit table

`run_108_genome_invariant_audit()` computes one structural profile for every method with
`NetworkProfiler` (raw draft, METEOR active set, Grow-Support network, GSPA network): dead-end metabolites
of the network, elementally unbalanced internal reactions of the network, unbalanced reactions among the
added ones, and FBA growth under the profile protocol. Before 2026-10-20 the first three rows of Table S12
carried archived MEMOTE counts over the whole network (exchanges included) while the GSPA row counted
unbalanced and dead-end reactions among its additions only; the columns `mean_imbalanced_added_or_active`
and `mean_deadend_active` now hold the same quantity for all rows. The per-genome table keeps the
archive-derived columns (`base_*`, `meteor_*_active`, `gs_imbalanced_active`, `gs_deadend_active`,
`slr_imbalanced_added`, `slr_deadend_added`) next to the uniform ones (`*_deadends_network`,
`*_imbalanced_internal`, `*_added_imbalanced_formula`).

Jobs (ibex): 54048078 (108-task extraction array, 20-100 s per genome, venv of 2026-10-19 with numpy 1.23.4,
scipy 1.13.1, pandas 1.5.3), 54048495 (merge and sha256 manifest). Earlier runs with the METEOR snapshot
code: 53986089, 53988457, 53988464 (work area `/ibex/scratch/projects/c2014/hohndor/meteor108_repro_2026-10-19`).

## Archive schema: `benchmark/manifests/reference_data/meteor108_drafts_compact.pkl.gz` (4.6 MB, schema 2)

Top level:
- `allrxns`: 47,880 reaction ids. Every `*_idx` array indexes this list.
- `media[gram]`: `lb`, `ub`, `biomass_id`, `biomass_idx`, `excluded_idx` and `media_exchange_idx`, with
  gram = negative or positive.
- `genomes[GCA]`, `panel`, `missing` (empty), `S_sha256`, `description`.

Per genome:
- `draft_idx`, `draft_maxscore`, `core_idx`
- `meteor_active_idx`, `meteor_v_idx`, `meteor_v_val`, `meteor_biomass_flux`, `meteor_status`
- `gs_repaired_idx`, `gs_added_idx`, `gs_support_idx`, `gs_v_val`, `gs_lp_objective`, `gs_*_growth_milp_medium`
- `profiles[repaired | raw_threshold | meteor]`: `n_rxn`, `n_internal`, `deadends`, `mass_imbal`,
  `mass_imbal_internal`, `fba_growth`, `deadends_shared_objects`
- `tgd` (archived T9 JSON), `tgd_match` (draft, core and added count equal), `tgd_match_draft_core`,
  `pred_path`, `n_proteins`, `implementation`

## Definitions

- **added** = internal reactions of the network (not `EX_`, not `biomass_*`) that are not in the
  tau = 0.5 draft. The archived Grow-Support `n_added_by_repair` (81.3) also counts the exchange
  reactions the repair opens; the internal-only count is about 53.
- **flux-certified added** = added reactions with |v| > 1e-6 in a pFBA solution: minimise sum |v|
  subject to S v = 0 and v_bio >= gamma_min = 0.1, with the network restricted to its own reactions.
- **zero-flux padding** = added minus flux-certified added.
- **network profile** = the selected reactions plus the universal model's single-metabolite boundary
  reactions of their metabolites plus the Gram biomass; internal reactions clipped to +-100 of their MILP
  bounds; exchanges first unbounded to find the uptakes an optimum uses, then those uptakes limited to 10
  and every other exchange closed for uptake.

Media:
- **METEOR's own medium** (`media[gram]`): the FVA-tight bounds plus the default medium.
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
  (p = 0.75, see `_solve_single_genome_slr`), not from the tau = 0.5 draft in the archive. The GSPA
  row therefore isolates the repair policy rather than an independent reconstruction.
- **Grow-Support certification.** All Grow-Support additions carry flux in Grow-Support's own LP by
  construction (support |v| > 1e-9). The pFBA columns test whether they are needed at minimum total flux.
- **Raw threshold feasibility.** For the raw threshold, `lp_biomass_feasible_pct` is the profile FBA growth
  > 0 and `lp_biomass_feasible_pct_own_medium` is the core-draft growth >= 0.1 on METEOR's medium.
- **Heavy compute kept as data.** The DeepProZyme score matrices and the METEOR MILP solutions are not
  recomputed; they enter as checksummed inputs (`manifest_rows.tsv`).
