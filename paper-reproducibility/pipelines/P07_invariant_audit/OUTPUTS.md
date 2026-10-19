# P07 outputs

| File | Written by | sha256 (2026-10-19) |
| --- | --- | --- |
| `results/pillar2_invariant_audit_108.tsv` | `run.sh audit` (108 rows, per-genome measured columns; `meteor_n_net_count_diff` keeps the old count difference) | c24ded4e21f54816ddee6e7e68655af1efd4f6e68bec31b1c7d474982706de7c |
| `results/pillar2_invariant_audit_summary.tsv` | `run.sh audit` (4 method rows; own-medium and audit-medium columns) | f20eed81d517df16924835dfad8058e7f3468eb19f00dfebf2de14925340fb5f |
| `results/pillar2_curated6_gems_real.tsv` | `run.sh audit` (byproduct, unchanged by P07) | |
| `benchmark/manifests/reference_data/meteor108_drafts_compact.pkl.gz` | ibex `merge.sbatch` then `run.sh fetch` | eaf6bf7ff7c7c702ceffae282b181b77fd2d2a8915d72d07b2fd800f8d54ca71 |
| `manifest_rows.tsv` | sha256 of every ibex input, computed inside job 53988464 (450 rows: 108 x {meteor_sol, dpz_pred, tgd_json, table1_v2_json} + code/data + output) | |

The `slr_solve_time_sec` column is wall-clock time and differs from run to run. Every other column is deterministic for fixed inputs (HiGHS LP).
