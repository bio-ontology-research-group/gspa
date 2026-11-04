# Off-growth v2 Scoring Memory and Runtime Fix

## Changes Implemented

The scoring and reporting pipeline in `benchmark/v3/offgrowth2` was updated across `report.py`, `score.py`, `rules.py`, and `score_driver.py`, along with unit tests in `tests/test_report.py` and `tests/test_score.py`:

- Two-phase reporting in `benchmark/v3/offgrowth2/report.py`: `generate_report` no longer loads all genome payloads and `GenomeScoringContext` objects into the parent process before forking. In the worker phase, each worker receives only a genome result file path, loads the payload, builds its context, scores all arms, writes `<out-dir>/per-genome/<genome>.json.gz`, clears per-genome caches, and exits after finishing its task. In the aggregation phase, the parent loads the compact per-genome files, resolves the `Q_fixed` threshold if needed, and writes `summary.json`, `arms.tsv`, `cases.tsv`, `decisions.tsv`, `pathway_calls.tsv`, `orphan_steps.tsv`, `cost.tsv`, and `report.md`.
- Resume support in `benchmark/v3/offgrowth2/report.py` and `benchmark/v3/offgrowth2/score_driver.py`: added `--resume` to both CLIs and `resume: bool = False` to `generate_report` so the worker phase skips any genome whose `<out-dir>/per-genome/<genome>.json.gz` is already present.
- Per-genome threshold support in `benchmark/v3/offgrowth2/score.py`: `GenomeScoreResult` stores compact `q_threshold_items` and `q_best_per_protein` records in each per-genome `.json.gz` file so `find_q_fixed_precision_threshold` and `apply_q_fixed_threshold_to_result` operate directly on the serialized per-genome results without keeping `GenomeScoringContext` instances in memory. When `--q-fixed-threshold` is explicitly passed, the threshold search is skipped.
- Per-genome and per-seed candidate caching in `benchmark/v3/offgrowth2/score.py` and `benchmark/v3/offgrowth2/rules.py`: `permute_candidates_map_within_genome`, `candidate_index`, `cached_best_by_protein`, and `_normalized_candidate_rows` are cached once per `(genome_id, seed)` across all runs of a genome while preserving fallback `id(candidates)` caching for unkeyed callers.
- Elimination of redundant per-run candidate scans in `benchmark/v3/offgrowth2/score.py` and `benchmark/v3/offgrowth2/rules.py`: `old_rule_calls`, `_evaluate_masked_recovery`, and `find_orphan_steps` reuse `candidate_index(candidates, genome_id, seed).ec_to_prot_scores` instead of scanning the full candidate dictionary on every arm and run, and `score_genome` caches reaction-gate calls and gene assignments across arms that share gate parameters.

## Verification Metrics

Table 1. Unit test and real-data verification runs (`OFFGROWTH2_SCORE_JOBS=8`).

| Verification run | Genomes | Tests | Output directory | Wall time (s) | User time (s) | System time (s) | Peak RSS (KB) | Exit code |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Unit test suite (`unittest discover`) | n/a | 23 | n/a | 3.09 | n/a | n/a | n/a | 0 |
| Subset Tier A `six` (`GCF_000236665.1`, `GCF_001633105.1`) | 2 | n/a | `/home/leechuck/tmp/gspa-v3-offgrowth2/verify-two/out` | 50.33 | 63.68 | 1.84 | 1141428 | 0 |
| Full Tier A `six` | 22 | n/a | `/home/leechuck/tmp/gspa-v3-offgrowth2/score-tierA-six/out` | 394.72 | 1367.51 | 31.62 | 2196116 | 0 |

Table 2. Output TSV line counts (`wc -l`, including header line) and data row counts for the two-genome and full Tier A `six` verification runs.

| Output table | Two-genome `wc -l` | Two-genome data rows | Full Tier A `wc -l` | Full Tier A data rows |
| --- | --- | --- | --- | --- |
| `arms.tsv` | 611 | 610 | 916 | 915 |
| `cases.tsv` | 652 | 651 | 17284 | 17283 |
| `cost.tsv` | 694 | 693 | 17746 | 17745 |
| `decisions.tsv` | 13078 | 13077 | 150997 | 150996 |
| `orphan_steps.tsv` | 507 | 506 | 5862 | 5861 |
| `pathway_calls.tsv` | 52269 | 52268 | 1306802 | 1306801 |
