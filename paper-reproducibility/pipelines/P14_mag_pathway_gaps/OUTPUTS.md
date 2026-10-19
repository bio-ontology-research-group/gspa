# P14 outputs

Everything is written under `OUT=/ibex/scratch/projects/c2014/hohndor/mags_repro_2026-10-19/`
(the Empty Quarter MAG tree `rob/metagenomics/`, the canonical `EmptyQuarter_Data/` tree, the
P13 reference directory and the gapseq conda env are only read).

| Path (ibex, under `$OUT`) | Produced by | Content |
| --- | --- | --- |
| `code/` | `run.sh sync` | this pipeline, `mag_pathway_gaps.py`, `bgc_benchmark.py` (claims builder only) |
| `ref/` | `00_select` | snapshot of all references (Swiss-Prot experimental-GO DIAMOND DB and its GO table, GO obo/owl/closure, pfam2go, ec2go, KEGG pathways/modules, gapseq SEED reaction graph files, essential-function profiles extracted from the jar) |
| `ref/sha256.txt`, `ref/tool_versions.txt` | `00_select` | checksums of jar, Pfam-A, MAG metadata tables, scripts and every reference; tool versions |
| `eligible_mags.tsv` | `00_select` | all 3,794 HQ MAGs passing the filters (CheckM2, GTDB phylum, completed gapseq, Bakta) |
| `mag_list.tsv`, `select.log` | `00_select` | the 100 sampled MAGs with taxonomy and input paths; allocation log |
| `work/<mag>/input_sha256.txt` | `01_mag` | sha256 of that MAG's FASTA, Bakta faa/gff3/fna and gapseq xml/tbl inputs (computed inside the job) |
| `work/<mag>/<mag>_layout.tsv`, `_operons.tsv` | `01_mag` (prep) | Bakta CDS coordinates (locus_tag ids) and same-strand ≤300 bp operons |
| `work/<mag>/<mag>_gap_reactions.tsv`, `_gaps.jsonl` | `01_mag` (prep) | gapseq gap reactions (final − draft model) with ECs, GO, MetaCyc pathways, gapseq GPR genes and sub-threshold hit loci; the `--gaps` input |
| `work/<mag>/<mag>_pathway_reactions.tsv` | `01_mag` (prep) | gapseq-predicted MetaCyc pathways and their EC-bearing reactions |
| `work/<mag>/<mag>_{diamond.tsv,pfam.domtbl,claims.jsonl}.gz` | `01_mag` | raw channel outputs and Stage 0 claims |
| `work/<mag>/<mag>_stage0.tsv`, `_stage12.tsv` | `01_mag` | GSPA posteriors: Stage 0 (noisy-OR) and Stage 1+2A (+gaps, operons, RLGC) |
| `work/<mag>/<mag>_rlgc_suggestions.tsv`, `_stage2_actions.tsv` | `01_mag` | RLGC suggestions per gap; Stage 2A actions |
| `work/<mag>/<mag>_score.tsv`, `_fills.tsv` | `01_mag` (score) | per-MAG metrics; every new fill with source and control class |
| `work/<mag>/timing.tsv`, `times/*.time` | `01_mag` | wall clock / max RSS per step (`/usr/bin/time -v`) |
| `results/pillar3_dark_matter_mags_v2.tsv` | `02_aggregate` | per-phylum table (+ `ALL` row) |
| `results/pillar3_dark_matter_mags_v2_per_mag.tsv` | `02_aggregate` | per-MAG table |
| `results/mag_gap_fills_v2.tsv.gz`, `timing_v2.tsv`, `input_sha256_per_mag_v2.tsv` | `02_aggregate` | all fills; timing; input checksums |
| `test/…`, `results/*_test*` | `run.sh test` | 3-MAG debug run |

`run.sh fetch` copies the two result tables to `results/` in the repository and the audit files
(fills, timing, checksums, versions, selection log, eligible list) into `outputs/` here.

## Columns of `pillar3_dark_matter_mags_v2.tsv`

`phylum_or_group, domain, n_mags`; `candidate_fills_q030`, `hq_fills_q050` (sum over MAGs of new
(gap reaction, locus) fills at posterior ≥0.30 / ≥0.50) and their per-MAG means;
`pathway_completeness_before/after` (Stage 0 vs Stage 1+2A+RLGC, q = 0.5, mean over MAGs of the
mean over gapseq-predicted MetaCyc pathways); `essential_completeness_before/after`;
`gapseq_pathway_completeness` (gapseq's own Completeness column, for reference);
`mean_gap_rxns_per_mag`, `mean_fillable_gaps_per_mag`; `gaps_resolved_q030/q050` (fillable gaps with
a new fill and no Stage 0 fill); `gaps_covered_stage0_q050`; `fills_prior_*` / `fills_rlgc_*`
(fill source; a fill can have both); and the control columns
`gapseq_gene_same_locus_*`, `gapseq_gene_other_locus_*`, `gapseq_subthr_hit_same_locus_*`,
`gapseq_subthr_hit_other_locus_*`, `gapseq_no_hit_*` (see README).
