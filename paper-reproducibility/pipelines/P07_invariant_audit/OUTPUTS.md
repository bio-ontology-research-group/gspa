# P07 outputs

| File | Written by | sha256 |
| --- | --- | --- |
| `results/pillar2_invariant_audit_108.tsv` | `run.sh audit` (108 rows, per-genome measured columns; `meteor_n_net_count_diff` keeps the old count difference) | a40b803eb9b28a8218d93fdfb8c01c55e4cae10b39827e44dd7967df6d29cca6 (2026-10-20) |
| `results/pillar2_invariant_audit_summary.tsv` | `run.sh audit` (4 method rows; own-medium and audit-medium columns) | ff2491cd0df8cb984799eaf2546b9db17ecdbfef9c02027738d52c223ac28484 (2026-10-20) |
| `results/pillar2_curated6_gems_real.tsv` | `run.sh audit` (byproduct, unchanged by P07) | |
| `benchmark/manifests/reference_data/meteor108_drafts_compact.pkl.gz` | ibex `merge.sbatch` (job 54048495, GSPA code) then `run.sh fetch` | 12604d59f0c2f9b89f93603d31cac99bc0ec8c6fb27d59fa921eec5c6dfc0809 (2026-10-20) |
| `manifest_rows.tsv` | sha256 of every ibex input, computed inside job 54048495 (442 rows: 108 x {meteor_sol, dpz_pred, tgd_json, table1_v2_json} + data files + panel + output) | |
| `results/smoke_test_reimpl.{tsv,summary.json,media.json}` | `run.sh smoke` (comparison of `scripts/meteor_threshold_repair.py` with the archived values, 108 genomes) | |

The `slr_solve_time_sec` column is wall-clock time and differs from run to run. Every other column is deterministic for fixed inputs (HiGHS LP).
