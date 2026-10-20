# P11 outputs

Everything below is computed by `benchmark/allthebacteria/aggregate_atb_java.py` from GSPA
integrator (gspa-unified.jar v2.1.0, sha256 `e415d9c6…72a7`) output only: `<gid>.noisy_or.tsv`
(Stage 0 noisy-OR, `max_iter 0`, no priors, no Stage 2) and `<gid>.gspa_full.tsv` (`max_iter 3`,
six priors, Stage 2A relax). No table here reads `delta_events.tsv`.

The older files in `results/allthebacteria/` without `_java` in the name (`insight1_…` to
`insight6_…`, `atb_insights_summary.json`) come from the Python heuristic in Stage C and are kept
only for comparison.

## Files in `results/allthebacteria/`

| File | Run | Genomes | Content |
|---|---|---|---|
| `insight1_dark_matter_rescue_by_phylum_java.tsv` | (a) | Java output from the current Stage C config | one row per domain × phylum × homology regime: protein counts, annotated before/after, rescued proteins and pairs by prior, suppressed pairs by category |
| `insight2_complex_coherence_vs_singleton_java.tsv` | (a) | same | one row per phylum × complex GO term (jar `complex_terms.tsv`, class h/n/a): coherence-rescued pairs, Stage 2A singleton-pruned pairs, number of genomes with each |
| `insight3_per_genome_java.tsv` | (a) | same | one row per genome: metadata, ORF/dark/twilight/zero-alignment counts, rescue/suppression counts, gaps.jsonl size, rows on which each prior fired, JVM runtime |
| `insight6_genome_size_strata_java.tsv` | (a) | same | six ORF-count strata: means/medians of per-genome rates; trailing `#` lines give the Spearman rho of each rate against n_orfs (all genomes, and below/above the median), with p and n |
| `atb_java_summary.json` | (a) | same | totals, regime fractions, rescue by prior, suppression categories, Spearman results, list of batches, date |
| `*_java_legacy.tsv`, `atb_java_summary_legacy.json` | (a-legacy) | the 3,312 genomes in w2_0013/0018/0019/0020/0021/0024 that an older Stage C integrated with `max_iter 5` and pairwise `operons.tsv` | same tables. Not pooled with (a), because the configuration differs |
| `*_java_stageD_snapshot.tsv`, `atb_java_summary_stageD_snapshot.json` | (b) snapshot | all 7,009 genomes of the 15 wave-1 batches Stage D had finished by 2026-10-19 14:04 | same tables |
| `heuristic_vs_java_same_genomes_java.tsv` | — | the 1,138 genomes present in both the heuristic insight3 table and run (a) | heuristic event totals against JVM totals (`compare_heuristic_vs_java.py`) |
| `*_java_stageD.tsv`, `atb_java_summary_stageD.json` | (b) | every genome of the batches Stage D has finished | same tables. Written when the Stage D arrays finish (see README for the command) |

The legacy batches w2_0001, w2_0009 and w2_0017 (6,854 genomes) have JVM outputs with
`n_claims = 0`: the older Stage C passed empty claims files. They are excluded everywhere.

## Column glossary (insight1 / insight3)

- `n_proteins` / `n_orfs`: pyrodigal ORFs (`dedup/genome_protein_to_rep.tsv`).
- homology regime: max % identity of the protein's linclust representative (90% id / 80% cov) to
  the CAFA training Swiss-Prot subset (`ref/train_db.dmnd`), DIAMOND blastp `--sensitive -k 25`, default
  E ≤ 1e-3. Regimes are `dark_matter_lt30` (<30, which includes zero-alignment), `twilight_30_50`,
  `moderate_50_70` and `close_homolog_ge70`. `zero_alignment` means no DIAMOND hit at all.
- `n_annot_noisy_or` / `n_annot_gspa_full`: proteins with at least one GO row with posterior ≥ 0.50.
- `n_rescued_proteins`: annotated in gspa_full but not in noisy_or. `n_lost_proteins` is the reverse.
- `rescued_pairs__X`: protein–GO pairs with p0 < 0.50 ≤ p1 where prior X has a positive entry in
  `priors_fired`. A pair counts under every prior that fired on it. `unattributed` means no positive
  entry; the shift is left over from a prior that fired in an earlier, damped iteration.
  `rescued_proteins__X` counts proteins with at least one such pair.
- `suppressed__*`: pairs with p0 ≥ 0.50 and p1 < 0.50, or removed from gspa_full. Each pair gets the
  first category that matches: `complex_singleton_prune`, `taxon_forbidden`, `anatomy_forbidden`,
  `sat_coannotation_conflict`, `stage2_removed_coannotation_conflict`, `anatomical_prior_other`,
  `other`. The aggregator's docstring defines each one.
- `taxid_in_jar_hierarchy`: 0 when the manifest taxid is not among the 1,302 taxa in the jar's
  bundled NCBI hierarchy. Such a genome has no organism-level forbidden terms, so only genome-wide
  co-annotation conflicts can suppress its terms.
- `n_gaps_jsonl`: number of KEGG-coverage gap claims that Stage C wrote and passed to the integrator.
  They are not gapseq output. The rule: KEGG pathways with between 40% and 100% (exclusive) of their
  GO terms covered at ≥ 0.45 by the Stage B scores; the missing reactions are listed, at most 350
  per genome. `n_gap_filling_fired_rows` and `n_gap_filling_rescued_pairs` count the gspa_full rows
  where the gap_filling prior fired, and the subset where it crossed 0.50.

## Intermediate files on Ibex (`/encrypted/e3008/allthebacteria/repro_2026-10-19/`)

- `survey/batch_inventory.tsv`, `survey/java_runtime.tsv`: batch inventory and every JVM runtime
  row joined with n_proteins (job 53986011).
- `test_bytecompat/batch_0001/out/*.tsv`: the 20-genome debug rerun (job 53986478). 40 of 40 outputs
  and 40 of 40 regenerated inputs are byte-identical to Stage C.
- `<batch>/out/<gid>.{noisy_or,gspa_full}.tsv.gz`, `<batch>/java_batch_stats.tsv`, `<batch>/DONE`:
  Stage D outputs.
- `partials_stageC/`, `partials_stageC_legacy/`, `partials_stageD/`: per-batch
  `<batch>.{genomes,regimes,complex}.tsv` written by `aggregate_atb_java.py summarize`.
- `ref/forbidden_by_taxon{,.sets}.tsv`, `ref/complex_terms.tsv`: constraint tables taken from the jar.
- `checksums.tsv`: sizes and sha256 of the inputs (job 53986494), copied into `manifest_rows.tsv`.
