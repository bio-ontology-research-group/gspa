# P12 outputs

## Results files (local repo)

| File | Content | Produced by |
| --- | --- | --- |
| `results/pillar1_suite2_1000bacteria_sota.tsv` | Suppl. Table S8, 8 rows (six tools, GSPA Stage 0, GSPA Stage 1+2). The historical columns are kept. Added columns: `source` (computed / unavailable), `n_genomes_evaluated`, `essential_full_n_terms`, the four taxon-violation variants, `smin`, the `gaef_*` columns (GSPA evaluator), `runtime_sec_per_genome_median`, `peak_ram_mb_max`, `unavailable_columns` and `notes` | `05_aggregate.py` → `run.sh fetch` |
| `paper-reproducibility/pipelines/P12_refseq1000/results/per_genome_suite.tsv` | One row per genome × method from the unified structural suite | `p12_structural_suite.py` via `03_eval_chunk.sh` |
| `…/results/per_genome_gaef.tsv` | One row per genome × method from `gspa evaluate --lite` | `03_eval_chunk.sh` |
| `…/results/per_genome_runtime.tsv` | GSPA wall time and max RSS for each genome and step | `02_gspa_genome.sh` (`/usr/bin/time -v`) |
| `…/results/matrix_reproduction.tsv` | The six tool rows recomputed from the audit's derived matrices | `p12_matrix_rows.py` |
| `…/results/progress.json` | Genomes finished for each method at aggregation time | `05_aggregate.py` |
| `…/results/{input_sha256_ibex.txt,tool_versions.txt,taxon_map.tsv}` | Checksums, DIAMOND/HMMER versions, organism taxon passed to GSPA | `01_prep_ibex.sbatch` |
| `…/manifest_rows.tsv` | Rows to append to `paper-reproducibility/manifest/data_manifest.tsv` (`used_by_pipeline = P12`) | by hand |

### Column semantics in the regenerated table

* `essential_full29_pct` keeps its historical name. It is the 32-term list
  (`essential_full_n_terms = 32`), propagated with GO 2025-03-16. The old
  value came from matrices propagated with an older GO release.
* `taxon_constrained_rows_per_genome` keeps its historical name, but now
  holds **violating** propagated rows per genome, which is what the paper's
  text says it reports. The old numbers (for example 14,198.84 for DeepGOMeta)
  were the count of rows whose GO term carries *any* taxon constraint,
  violated or not. `p12_matrix_rows.py` recomputes the old numbers exactly
  under that definition (`taxon_constrained_rows_per_genome` in
  `matrix_reproduction.tsv`).
* `cafa_macro_fmax` and `smin` are `NA` for every row (see README, "CAFA truth").
* `runtime_sec_per_genome` and `peak_ram_mb` are `NA` for the six tools,
  because the audit logged neither.

## Raw outputs (IBEX)

`/ibex/scratch/projects/c2014/hohndor/refseq1000_repro_2026-10-19/`:

* `work/<genome>/`:
  * `stage0_integrated.tsv`, `stage12_integrated.tsv`: integrator output
  * `stage12_actions.tsv`: Stage 2A enforcement actions
  * `claims_{dp,ipr}.jsonl.gz`: claims
  * `preds/{diamond_results.tsv.gz,pfam_results.domtbl.gz}`: DIAMOND and HMMER results
  * `times/*.time`: timing logs
  * `logs_*.txt`: step logs
  * `DONE`: completion marker
* `eval/{tools,gspa}/`:
  * `suite/chunk_K.jsonl`: unified suite records
  * `leaf/<g>.<method>.leaf.tsv`: leaf GO sets used by the GSPA evaluator
  * `gaef/<g>.<method>_quality.json` and `.gaf`: GSPA evaluator reports
  * `manifests/`: evaluator batch manifests
* `results/matrix_rows.json`: matrix reproduction
* `results/aggregate/`: aggregator outputs
* `gff/`, `operons/`: RefSeq GFF3 and derived operons
* `db/`, `ref/`, `bin/`: inputs (see `manifest_rows.tsv`)

unimatrix01 `/data/hohndor/gspa/refseq1000_repro_2026-10-19/`:

* `goa_reference_loo9_subset.gaf`: GOA subset
* `truth/exp_truth_rows.tsv`, `truth/exp_truth_summary.tsv`: experimental-truth check
* `input_sha256.txt`: input checksums

## Job log (2026-10-19)

| Job | Host | Purpose | State |
| --- | --- | --- | --- |
| 9254 | unimatrix01 | GOA subset (first try; failed instantly on node007) | FAILED |
| 9258 | unimatrix01 | GOA subset, 857,017,321 lines scanned, 6,936,109 kept | COMPLETED |
| 9265 | unimatrix01 | truth check, landed on the head node, cancelled | CANCELLED |
| 9267 | unimatrix01 | truth check (`04_goa_exp_truth.py`) | COMPLETED |
| 53986257 | IBEX | prep (16 cores, cancelled and resubmitted with 4) | CANCELLED |
| 53986264 | IBEX | prep: makedb, hmmpress, 999/1,000 GFF3, operons, taxon map (20 min 53 s) | COMPLETED |
| 53986354 | IBEX | matrix reproduction (Python 3.7 lacks `statistics.fmean`) | FAILED |
| 53987063 | IBEX | matrix reproduction | COMPLETED |
| 53986356 | IBEX debug | evaluation test, tools, 2 genomes | COMPLETED |
| 53986642 | IBEX | evaluation, tools, 40 × 25 genomes (DeepFRI GAEF empty, see README) | COMPLETED |
| 53997270 | IBEX | GSPA-evaluator re-run for DeepFRI only (`GAEF_METHODS=deepfri`) | COMPLETED (40/40 chunks) |
| 53987583 | IBEX debug | GSPA test, genomes 1–5 (tasks 0–3 completed; task 4 never started and was cancelled; its genome is covered by the array) | 4/5 COMPLETED |
| 53989989 | IBEX debug | evaluation test, GSPA rows, 4 genomes | COMPLETED |
| 53989601 | IBEX | GSPA array, 100 tasks × 10 genomes, at most 60 concurrent | COMPLETED (1,000/1,000 genomes with `DONE`, no `[FAIL]` line in any task log) |
| 54026849 | IBEX | evaluation, GSPA rows, 40 × 25 genomes (submitted after the array finished) | COMPLETED (40/40 chunks, 2,000 genome × configuration rows) |
| 54027265 | IBEX | `05_aggregate.py`, final (all eight rows `source=computed`; fetched into `results/`) | COMPLETED |
| 53991948, 53997214, 54002223, 54002429 | IBEX debug | `05_aggregate.py` snapshots (the last one is in `results/`) | COMPLETED |

## Finishing the run

```bash
P=paper-reproducibility/pipelines/P12_refseq1000
ssh ibex 'ls /ibex/scratch/projects/c2014/hohndor/refseq1000_repro_2026-10-19/work/*/DONE | wc -l'   # want 1000
ssh ibex 'squeue -u hohndor -r | grep -c p12_'                                                          # want 0
bash $P/run.sh gspa 60        # only if genomes are missing after the array ends (DONE genomes are skipped)
bash $P/run.sh eval gspa      # structural suite + gspa evaluate for the two GSPA rows
bash $P/run.sh aggregate      # then, when the job has finished:
bash $P/run.sh fetch          # copies the table into results/pillar1_suite2_1000bacteria_sota.tsv
```

`05_aggregate.py` writes a row (and each `gaef_*` / runtime group) only when
all 1,000 genomes are evaluated (`--min-genomes 1000`). Until then the row
has `source = unavailable`, and the partial values are in the per-genome TSVs.
