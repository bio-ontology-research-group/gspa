# Many genomes in batch

`gspa integrate --batch-manifest` runs the same integration stages as a
single-genome `gspa integrate` invocation. Enabled reaction-context, dark-matter,
recurrence, symbolic-enforcement, FBA and agent stages run for each genome.
GO, anatomical constraints, complex definitions, pathway databases and the
reaction graph are shared; per-genome state and taxon checkers are separate.

```bash
gspa integrate --batch-manifest manifest.tsv --batch-threads 8 \
  --batch-stats-out batch_performance.tsv --go-owl reference_data/go.owl --lite
```

The manifest is a positional TSV. The first three columns are required; trailing
columns are optional. Existing 11-column manifests remain valid. Blank optional
fields or `-` inherit the command-line value. Paths are relative to the working
directory of the command. Genome IDs and output paths must be unique.

| Column | Name | Meaning |
|---|---|---|
| 0 | `genome_id` | Unique genome or sample ID |
| 1 | `claims_path` | Evidence claims JSONL |
| 2 | `out_path` | Integrated TSV output |
| 3 | `taxon` | Taxon ID or domain alias |
| 4 | `essential_profile` | Essential-function profile |
| 5 | `operons_path` | Operon assignments TSV |
| 6 | `gaps_path` | Metabolic gaps JSONL |
| 7 | `max_iter` | Refinement iteration limit |
| 8 | `enforce_stage2` | `true` or `false` |
| 9 | `enforce_mode` | Symbolic enforcement mode |
| 10 | `enable_priors` | Comma-separated prior names; `none` disables priors |
| 11 | `genome_layout_path` | Gene coordinates for reaction context |
| 12 | `orthogroups_path` | Protein-to-orthogroup assignments |
| 13 | `cluster_consensus_path` | Cluster functional consensus |
| 14 | `fba_medium` | FBA medium name |
| 15 | `fba_gram` | FBA biomass template (`neg` or `pos`) |

A header beginning `genome_id`, `tag` or `id` and followed by `claims_path` or
`claims` is recognized. Lines beginning with `#` are ignored. Files used for
orthogroups and consensus must be prepared without held-out annotation labels
when evaluating prediction accuracy.

Each integrated TSV has a companion `<out_path>.stages.json` ledger. It records
executed, disabled, failed and unreached stages, plus score increases, decreases,
additions and removals for individual gene/function pairs. These are prediction
changes, not measurements of biological accuracy. An explicitly enabled stage
with missing required inputs fails rather than silently claiming to have run.

With `--cohort` a batch is a **two-pass** run: the stage-0 raw baseline of every row (and of the genomes in
`--cohort-manifest`) is built first, one joint cohort inference follows, and each genome's workflow then consumes only
its own held-out features (`<out_path>.cohort.json` reports the snapshot, peers and gaps). Rows may share protein ids,
because every candidate is namespaced by genome; claims, layout and orthogroups must use the FAA sequence ids. If the
cohort pass fails, no genome runs and every ledger records `cohort_context: failed`. `--feature-export-dir` writes one
sub-directory per genome; `--community-manifest` runs the community once and every ledger points at that result. See
[Learned, cohort and metabolic context](../concepts/context-workflow.md).

Optional artifact flags select separate files for each batch row:
`<out_path>.provenance.json`, `.suggestions.tsv`, `.features.tsv`,
`.enforcement.tsv`, `.anomalies.jsonl`, and `.agent_resolutions.jsonl`.
The filename supplied to such a flag is used for single-genome runs; batch mode
uses these per-output suffixes to prevent workers overwriting each other. FBA
artifacts go in `<out_path>.metabolism/` in both modes. Agent execution uses
`--agent-python-bin` (default `python3`); install its dependencies in that
interpreter's environment before running it.

`--batch-stats-out` records claim and candidate counts, iteration settings,
per-genome elapsed time, and JVM heap usage at completion. Heap measurements are
for the shared JVM, not the memory attributable to one genome. Set thread counts
and JVM heap limits for the available hardware and enabled stages; v3 throughput
claims require the production benchmarks and are not yet established.
