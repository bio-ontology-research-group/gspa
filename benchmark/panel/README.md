# Panel expansion — KAUST culture genomes

Phase 12 (cross-genome conditional-LR dark-matter) scales from the 29-
genome GTDB panel + 4-culture pilot to the ~340 KAUST-isolate/MAG
genomes under `/data/emptyquarter/sequencing-results/` on unimatrix01.

## Strategy (approved 2026-04-20)

- **Strategy C** (tiered panel): 29 GTDB reps + ~100 KAUST isolates/MAGs
  stratified by phylum as the **panel**; remaining genomes as
  **query**.
- Dereplication: strict **ANI-95%** (skani).
- Scope: full — expect ~150-200 genomes post-derep.

## Phase 1 — inventory + QC + dereplication → manifest

Pipeline (run on unimatrix01):

1. `enumerate_genomes.sh` — walks `isolates/`, `enrichment/`,
   `rh_sequencing/`; emits `genome_inventory.tsv` (one row per
   candidate assembly FASTA). Applies per-source "canonical FASTA"
   rules so a single sample never contributes multiple redundant rows.
2. `stage_fastas.sh` — turns inventory into symlinks
   `staged/<genome_id>.fna`, writes `genome_list.tsv`.
3. `run_checkm2.sh` (SLURM) — CheckM2 v1.1.0 on staged dir;
   `quality_report.tsv` with completeness / contamination.
4. `run_gtdbtk.sh` (SLURM) — GTDB-Tk `classify_wf --skip_ani_screen`
   on staged dir; phylum-level taxonomy.
5. `run_skani_derep.sh` (SLURM) — skani all-vs-all ANI, then
   greedy-cluster at ANI ≥ 95, AF ≥ 50; cluster reps preferred by
   contig count.
6. `build_genome_manifest.py` (SLURM, depends on 3+4+5) — joins
   everything into `genome_manifest.tsv` with per-genome completeness,
   contamination, classification, cluster rep, quality tier
   (high/medium/low/excluded).

Orchestrated by `phase1_orchestrate.sh` which submits all SLURM jobs
with the right dependencies.

### Output

`/data/hohndor/gspa/proteomes/culture_panel/phase1/`:

```
genome_inventory.tsv       # raw enumeration (source_dir, sample_id, path, size)
genome_list.tsv            # genome_id → original path map
staged/<genome_id>.fna     # symlinked flat FASTA dir
checkm2/quality_report.tsv
gtdbtk/{bac120,ar53}.summary.tsv
skani/{ani_pairs.tsv,clusters.tsv}
genome_manifest.tsv        # final per-genome manifest
```

## Phase 2+ (downstream — not in this dir)

Phase 2 runs `annotate_culture.sh` as a SLURM array over the dereplicated
representative set, then re-clusters proteins with MMseqs2. Phase 3
rebuilds the non-anchor catalog. See `benchmark/cultures/` and
`benchmark/cross_genome/` for the predict + augment scripts.
