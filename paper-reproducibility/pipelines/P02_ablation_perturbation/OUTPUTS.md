# P02 outputs

## Results files (local repo)

| File | Content | Produced by |
| --- | --- | --- |
| `results/pillar1_2_perturbation_designs.tsv` | All five perturbation runs, 9 genomes × 5 configs each (225 rows). Columns: `design jar claims_ids genome config injected_pairs injected_claims surviving_violations surviving_anatomical surviving_taxon suppression_rate_pct fmax_micro_all fmax_micro_exp fmax_cafa_all fmax_cafa_exp` | `run_perturbation.py` ×5 → `collect_results.py` |
| `results/pillar1_2_perturbation_summary.tsv` | **Current numbers**: design v2, current jar, `claims_ids=normalised` (45 rows). Historical column layout with a leading `design` column (`injected_false_positives` = `injected_pairs`) | same |
| `results/pillar1_2_ablation_summary.tsv`, `results/pillar1_2_ablation_per_genome.tsv` | Ablation grid (not re-run; shown below to be unaffected by the claims-ID bug) | `ablation/01`–`03`, `ablation/05_compile_results.py` |
| `ablation/ablation_idfix_compare_2026-10-19.tsv` (this folder) | Ablation A1–A9 × 9 genomes re-integrated on normalised claims vs the published outputs | `ablation/ablation_idfix_check.py` |

Runs in `pillar1_2_perturbation_designs.tsv`:

| design | jar | claims_ids | run dir | meaning |
| --- | --- | --- | --- | --- |
| v2 | gspa-unified.jar | normalised | `v2_current_jar_idfix` | **current injection, current integrator, fixed IDs (report this)** |
| v1 | gspa-unified.jar | normalised | `v1_current_jar_idfix` | **old injection, current integrator, fixed IDs (report this)** |
| v2 | gspa-unified.jar | raw | `v2_current_jar` | reproduces the pre-fix paper numbers 612/109/109/21/21 |
| v1 | gspa-unified.jar | raw | `v1_current_jar` | old injection, current integrator, pre-fix |
| v1 | gspa-unified-bench.jar | raw | `v1_old_jar` | reproduction check of the 326c112 numbers (old jar, one JVM per job) |

`claims_ids=raw` rows are affected by the phantom-ID bug (below) for hpylori,
mgenitalium and mjannaschii and should not be reported as results.

F-max values are rounded to 4 decimals by `collect_results.py`; the per-run
`perturbation_summary.tsv` on unimatrix01 keeps full precision.

## Raw outputs (unimatrix01)

`/data/hohndor/gspa/proteomes/bench9_unified/repro_2026-10-19/<run dir>/`:

* `{tag}_injected_pairs.tsv` — 100 injected (protein, GO term) pairs with category and primary raw score
* `{tag}_claims_perturbed.jsonl` — base claims (preprocessed when `claims_ids=normalised`) + injected claims (100 for v1, 160 for v2)
* `claims_preprocessing.tsv` — per genome: lines in, IDs normalised, duplicates removed, lines out, proteins (idfix runs only)
* `{tag}_{config}_perturbed_integrated.tsv` — integrator output
* `{tag}_{config}_cafa_eval.json` — `benchmark_pgap_v2.py` output (exp and all truth)
* `{tag}_A8/A9_perturbed_actions.tsv` — Stage 2A enforcement actions (per-job runner only)
* `perturb_batch_manifest.tsv` (batch runner), `logs/`
* `perturbation_summary.{tsv,json}`, `run_info.json` (claims_ids, jar sha256, go.owl sha256, host, SLURM job id, timings)

`repro_2026-10-19/ablation_idfix_check/`: normalised claims, batch manifest,
81 integrated TSVs + CAFA JSONs, `ablation_idfix_compare.tsv`.
SLURM stdout/stderr: `repro_2026-10-19/slurm/p02_<run>.{out,err}`.

## Regenerating

```bash
# on unimatrix01, after copying this folder to $PIPE (see run.sh)
bash $PIPE/run.sh submit_idfix    # v1 + v2 (normalised claims) + ablation ID check, ~10 min each
bash $PIPE/run.sh submit          # optional: the three historical raw-ID runs (old-jar run ~36 min)
# locally
R=unimatrix01:/data/hohndor/gspa/proteomes/bench9_unified/repro_2026-10-19
for r in v1_old_jar v1_current_jar v2_current_jar v1_current_jar_idfix v2_current_jar_idfix; do
  mkdir -p runs/$r; scp $R/$r/perturbation_summary.tsv runs/$r/
done
python3 paper-reproducibility/pipelines/P02_ablation_perturbation/collect_results.py \
  --runs runs/v1_old_jar runs/v1_current_jar runs/v2_current_jar \
         runs/v1_current_jar_idfix runs/v2_current_jar_idfix --results-dir results
```

## Claims-ID bug and fix (2026-10-19)

In `bench9/with_iprscan/{hpylori,mgenitalium,mjannaschii}_claims_merged.jsonl`
every pfam claim uses the raw FASTA ID (`sp|P21762|AHPC_HELPY`; 11,325 / 5,739
/ 13,773 claims). The other six genomes have no such IDs. The injection step
sampled proteins from the raw ID set, so 47 / 42 / 55 of the 100 injected pairs
were placed on phantom IDs. `gspa integrate` normalises IDs, so those pairs
were integrated under the accession but never matched by the survivor count,
which undercounted survivors for these three genomes.

Fix (`run_perturbation.py --claims-ids normalised`, the default): before sampling,
all nine claims files are ID-normalised (`^(sp|tr)|ACC|NAME` → `ACC`) and exact
duplicate records are removed (6–96 per genome; normalisation creates no new
duplicates). For the six unaffected genomes the injected pairs are unchanged and
every survivor and F-max value is identical to the raw run, so de-duplication has
no effect on results.

### Perturbation, normalised claims, current jar (jobs 9305 v1, 9306 v2)

Survivors summed over nine genomes (out of 900), anatomical / taxon in brackets:

| config | v1 (job 9305, 629 s) | v2 (job 9306, 596 s) | v2 raw (pre-fix) |
| --- | --- | --- | --- |
| A1 | 149 (71/78) | 735 (371/364) | 612 |
| A6 | 2 (2/0) | 137 (43/94) | 109 |
| A7 | 10 (2/8) | 137 (43/94) | 109 |
| A8 | 2 (2/0) | 24 (24/0) | 21 |
| A9 | 2 (2/0) | 24 (24/0) | 21 |

Affected genomes, survivors raw → normalised:

| genome | design | A1 | A6 | A7 | A8 | A9 |
| --- | --- | --- | --- | --- | --- | --- |
| hpylori | v2 | 41 → 77 | 6 → 14 | 6 → 14 | 3 → 2 | 3 → 2 |
| mgenitalium | v2 | 43 → 88 | 5 → 16 | 5 → 16 | 0 → 3 | 0 → 3 |
| mjannaschii | v2 | 36 → 78 | 6 → 15 | 6 → 15 | 2 → 3 | 2 → 3 |
| hpylori | v1 | 8 → 14 | 0 → 0 | 0 → 1 | 0 → 0 | 0 → 0 |
| mgenitalium | v1 | 6 → 22 | 0 → 0 | 0 → 1 | 0 → 0 | 0 → 0 |
| mjannaschii | v1 | 13 → 12 | 0 → 0 | 1 → 1 | 0 → 0 | 0 → 0 |

### Ablation grid (job 9307, 651 s)

`ablation/ablation_idfix_check.py` re-ran the full ablation manifest (A1–A9 × 9
genomes, same jar and references) on normalised, de-duplicated claims. All 81
integrated outputs match the published `bench9_unified/integrated/` files
exactly: same rows, max |Δposterior| = 0, no call flips, ΔF-max = 0 (micro and
CAFA, exp and all). GAEF metrics are computed from these same integrated files,
so they are unchanged too. The ablation does not inject claims, and integrate
normalises IDs internally, so the bug does not affect it.

## Validation of the historical raw-ID runs (jobs 9231, 9234, 9235)

| run | A1 | A6 | A7 | A8 | A9 | expected |
| --- | --- | --- | --- | --- | --- | --- |
| v2 / gspa-unified.jar / raw (job 9234, 446 s) | 612 | 109 | 109 | 21 | 21 | 612/109/109/21/21. Reproduced exactly; all 45 rows identical to the pre-fix summary after 4-dp rounding |
| v1 / gspa-unified-bench.jar / raw (job 9235, 2174 s) | 128 | 2 | 2 | 2 | 2 | 125/3/–/3/3 (326c112). **Not reproduced exactly** |
| v1 / gspa-unified.jar / raw (job 9231, 466 s) | 128 | 2 | 8 | 2 | 2 | n/a (first run of this combination) |

The v1 / old-jar deviation comes from the injected set, not the integrator. A1
(no priors) gives the same per-genome survivors under both jars (128), while the
326c112 table has 125 with a different per-genome split (e.g. mgenitalium
14 vs 6, mtb 12 vs 18). Claim and (protein, function) group counts in the
Oct 14 integrate logs (`bench9_unified/logs/*_perturbed.log`) match this run
for all nine genomes, and mean F-max per config matches the 326c112 table to
within 1e-4. Only the drawn proteins and scores differ. The 326c112 script
reuses an existing `{tag}_injected_pairs.tsv` if one is present, so its numbers
probably came from an injected set written by an earlier version of the script.
Those files were overwritten on 2026-10-17 14:24 and cannot be checked.
