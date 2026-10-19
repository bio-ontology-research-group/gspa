# P15 outputs

Input: `results/multikingdom_per_genome_eval.tsv` (88 proteomes x 3 configurations, written on
unimatrix01 by `scripts/run_multikingdom_batch_unimatrix.py`; not re-run here).

| File | Written by | Content |
| --- | --- | --- |
| `results/multikingdom_bootstrap_ci.tsv` | `python3 bootstrap_ci.py --repo .` (deterministic: 10,000 resamples, seed 20261019) | Per domain (and `ALL`), GSPA configuration (`GSPA_Iter1`, `GSPA_Iter3_Relax`) and metric (`fmax_all`, `fmax_lt30`, `fmax_mf`, `fmax_bp`, `fmax_cc`): mean of the per-proteome difference from `Noisy-OR`, 95% paired percentile bootstrap interval over proteomes, and counts of proteomes with a positive / zero / negative difference |

Used by Suppl. Table S9 (`tab:supp_multikingdom_domains`, `paper-reproducibility/tables/make_tables.py`).
