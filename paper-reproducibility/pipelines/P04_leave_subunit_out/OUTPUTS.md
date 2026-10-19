# P04 outputs

Remote root: `unimatrix01:/data/hohndor/gspa/proteomes/pillar3_unified/repro_2026-10-19/`
(`$REPRO`). Nothing outside `$REPRO` is written; the Oct-15/17 files in `pillar3_unified/` are
read only.

## Copied into the repository (`bash run.sh fetch`)

| File | Content |
| --- | --- |
| `results/pillar3_complex_lso.tsv` | One row per decoy_design × genome/cohort × method: `decoy_design, genome, method, method_id, n_cases, hit_at_1, hit_at_1_strict, recall_at_05, precision_at_05, f1_at_05, complex_coherence`. Cohorts: six genomes, `ALL (Weighted Mean)`, `Out-of-Sample 5 (Weighted Mean)` (all but ecoli). Superset of the Oct-17 table except that the synthetic "DIAMOND-BLASTp" / "InterProScan / Pfam" rows are replaced by `input_attenuated`, `real_diamond_channel`, `real_interproscan_channel`; `coherence_operon_stage2` is new. |
| `results/pillar3_complex_lso_mcnemar.tsv` | Exact McNemar (two-sided binomial on discordant pairs, strict Hit@1) per decoy design and cohort (ALL, Out-of-Sample 5, synechocystis): coherence_operon vs input_attenuated, coherence_operon vs coherence, full_stage3 vs coherence_operon. Columns `b_a_only` (method_a right, method_b wrong), `c_b_only`, `p_exact_two_sided`. |

## Remote only (`$REPRO`)

| Path | Content |
| --- | --- |
| `regen_check.tsv` | Whether the original generator (seed 42) reproduces the frozen Oct-15 cases and the complex term/module files (all `True` on 2026-10-19) |
| `complex_terms_extended.tsv`, `complex_modules.tsv` | Regenerated; byte-identical to `pillar3_unified/` copies |
| `within_operon_decoy_stats.tsv` | Per genome: decoys from operon / directon / fill, cases needing fill, synthetic claims removed, pre-existing real claims on new decoys |
| `{design}/{g}_m1_baseline.tsv`, `_m2_coherence.tsv`, `_m3_coherence_operon.tsv`, `_m3_stage2.tsv`, `_b3_operon_vote.tsv` | Jar outputs (unrelated_loci reruns are byte-identical to the Oct-15/17 files) |
| `{design}/{g}_lso_anomalies.jsonl`, `{g}_lso_resolutions.jsonl`, `{g}_m4_agent.tsv` | Stage 3 input, agent resolutions, Stage 2 output with resolutions applied |
| `within_operon/{g}_lso_cases.json`, `within_operon/{g}_claims_lso.jsonl` | within_operon design inputs |
| `pillar3_complex_lso.tsv`, `pillar3_complex_lso_mcnemar.tsv` | Same as the repository copies |
| `pillar3_complex_lso_per_case.tsv` | Per case × method scores and indicators (input to McNemar) |
| `real_channel_coverage.tsv` | Held-out subunits / decoys with a real DIAMOND or InterProScan claim for the complex term |
| `stage3_status.tsv` | PROMOTE / DEMOTE / ABSTAIN / other counts and guardrail rejections per design × genome |
| `P04_inputs.sha256.tsv` | path, size, sha256, mtime of every input (`sha_inputs.sh`) |
| `slurm/p04_lso.{out,err}` | Job log |

## Copies of the small diagnostics in this folder (`aux/`)

`pillar3_complex_lso_per_case.tsv`, `real_channel_coverage.tsv`, `stage3_status.tsv`,
`within_operon_decoy_stats.tsv`, `regen_check.tsv`, `P04_inputs.sha256.tsv` (from job 9245 and
the `sha_inputs.sh` job, 2026-10-19).
