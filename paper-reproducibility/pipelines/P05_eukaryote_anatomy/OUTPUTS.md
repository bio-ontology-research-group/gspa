# P05 outputs

All three tables are written by `benchmark/unified/pillar3/run_pillar3_eukaryotic_anatomy.py`
(`run.sh run both`) to `$WORK/results/` on unimatrix01
(`/data/hohndor/gspa/proteomes/pillar3_unified/repro_2026-10-19/P05_eukaryote_anatomy/work/results/`)
and copied to `results/` by `run.sh fetch`. Every table has a leading `input_source` column:

- `real`: stress test built on the real blind DIAMOND + Pfam claims (the default; use these rows).
- `synthetic`: the earlier design with claims generated from Swiss-Prot truth. Kept so the earlier
  numbers stay reproducible and the two designs can be compared.

All other legacy column names are unchanged. New columns are appended after them.

| File | Rows | Content |
| --- | --- | --- |
| `results/pillar3_eukaryote_anatomy.tsv` | 6 species x 2 input sources | Violations per configuration (B1 raw DIAMOND, E1..E5) on the perturbed cohort, the share eliminated, developmental rescue (attenuated true pairs reaching posterior >= 0.50 under E1/E3/E5), injection diagnostics (`n_injected_forbidden_for_species`, `n_injected_already_true`), how many attenuations were possible (`n_dev_true_pairs`, `n_attenuation_eligible`), violations among the injected pairs (`injected_pair_violations_e*`), and which constraint tables were used (`constraint_tables`) |
| `results/pillar3_relax_vs_remove.tsv` | species x mode (B1, E1..E5) x scope | Micro F-max (`cafa_fmax`, the legacy definition: pooled TP/FP/FN, thresholds 0.1..0.9), propagated precision/recall, mean IC of raw terms with posterior >= 0.50, plus CAFA protein-centric F-max (`fmax_protein_centric`, thresholds 0.05..0.95). `scope` is one of `injected_cohort_all_proteins`, `injected_cohort_injected_proteins` (the only scope in the earlier table), or `endogenous_cohort_all_proteins` (real only) |
| `results/pillar3_eukaryote_endogenous_vs_injected.tsv` | 6 species + `MEAN_6_EUKARYOTES`, per input source | Endogenous (unperturbed real claims) vs injected summary. In `real` rows, `endog_e4_remove_fmax` comes from an actual E4 run (`endog_e4_method = measured_run`). In `synthetic` rows it is the legacy formula `max(E1, E5 - 0.0008)`, flagged as such. `endog_*_fmax` uses CAFA protein-centric F-max over all Swiss-Prot truth proteins of the panel (n = 1,200; 155 for P. patens), the same definition as `multikingdom_per_genome_eval.tsv`. `inj_e*_fmax/recall/mean_ic` cover the injected proteins only (legacy scope). `inj_all_e*` cover every protein in the perturbed cohort |

Deleted: `results/pillar3_eukaryotic_anatomy.tsv`, a stale duplicate from an older run of the synthetic
design with different constraint tables.

Diagnostics (remote only, in `$WORK`):

- `run_info.json`: argv, host, SLURM job, sha256 of the jar, go.obo and go.owl, sha256 of every
  input claims file and every perturbed claims file, and per-species perturbation counts.
- `diagnostics.json`: a cross-check that the generic evaluator matches the legacy evaluator on the
  synthetic rows, and the multikingdom noisy-OR outputs re-evaluated with this evaluator.
- `real/{species}_injections.tsv`: each injected pair, with flags for whether it is forbidden for the
  species and whether it is already true.
- `real/{species}_attenuated.tsv`: each attenuated true developmental pair.
- `real/{species}_{inj,endog}_e{1..5}.tsv`: the integrated outputs, each with a `.stamp.json`
  (claims sha256 plus the exact command) and a `.log`.
- `synthetic/{species}_e{1..5}.tsv`: the same for the synthetic design.
