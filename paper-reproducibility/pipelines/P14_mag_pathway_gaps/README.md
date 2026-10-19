# P14 — Pathway gap resolution across uncultured MAGs (Empty Quarter HQ MAGs)

Regenerates, as a real experiment, the analysis behind Suppl. Table S11 ("pathway gap
resolution across uncultured metagenome-assembled genomes"). The old table
`results/pillar3_dark_matter_97mags.tsv` is a static transcription of a `unimatrix01` run that
cannot be traced to inputs or code. It is left untouched. The new outputs are
`results/pillar3_dark_matter_mags_v2.tsv` (per phylum) and
`results/pillar3_dark_matter_mags_v2_per_mag.tsv` (per MAG).

Code: `benchmark/unified/pillar3/mag_pathway_gaps.py` (stdlib Python; subcommands `select`,
`prep`, `score`, `aggregate`). The module docstring holds the exact definitions. Job scripts
live in `jobs/`, and `run.sh` is the driver. Ibex output root:
`/ibex/scratch/projects/c2014/hohndor/mags_repro_2026-10-19/`.

## How to run

```bash
bash paper-reproducibility/pipelines/P14_mag_pathway_gaps/run.sh sync        # laptop -> ibex code/
ssh ibex 'bash /ibex/scratch/projects/c2014/hohndor/mags_repro_2026-10-19/code/run.sh select'   # 00 (debug)
ssh ibex 'bash .../code/run.sh test'       # 3 MAGs on debug (Actinomycetota, Pseudomonadota, 1 archaeon)
ssh ibex 'CONC=10 bash .../code/run.sh all'  # 100-MAG array, 10 concurrent, then aggregate (afterany)
ssh ibex 'bash .../code/run.sh resume'     # resubmit only MAGs without work/<mag>/DONE, then aggregate
ssh ibex 'bash .../code/run.sh rescore'    # scoring changes only: re-score finished MAGs, re-aggregate
bash paper-reproducibility/pipelines/P14_mag_pathway_gaps/run.sh fetch       # results/ + outputs/
```

All computation runs under sbatch (account `pi-hohndor`). The login node only submits jobs.

## MAG selection (`00_select`, seed 20261019)

Source list: `rob/metagenomics/05_metabolism/hq_mag_list.tsv` (3,835 single-sample SemiBin MAGs
labelled HQ; these are the MAGs gapseq was run on). The selection script keeps a MAG if it meets
all of the following:

* CheckM2 completeness ≥ 90 % and contamination ≤ 5 %
  (`rob/metagenomics/results/checkm2_single_sample_full.tsv`). All 3,835 MAGs pass.
* It has a GTDB-Tk r226 classification with a phylum (`02_taxonomy/gtdbtk_out/gtdbtk.{bac120,ar53}.summary.tsv`).
  All 3,835 MAGs pass.
* It has a completed gapseq run: `<mag>.xml`, `<mag>-draft.xml`, `-all-Pathways.tbl` and
  `-all-Reactions.tbl` exist and are non-empty in `05_metabolism/gapseq_out/<mag>/`.
  41 MAGs fail this check.
* It has a Bakta 1.x annotation (`.faa`, `.gff3`, `.fna`) under
  `EmptyQuarter_Data/.../<sample>/wgs/annotation/single_sample/bakta/<bin>/`.
  Bakta's CDS calls come from Pyrodigal, and these are the "Prodigal" proteome and GFF used here.

This leaves **3,794 eligible MAGs in 35 phyla** (`eligible_mags.tsv`). That is more than 150, so
the script draws a phylum-stratified sample of **100 MAGs**:

* The 17 phyla with ≥ 10 eligible MAGs are the strata.
* Each stratum gets 4 slots. The remaining 32 slots go to strata in proportion to their
  remaining size (largest remainder).
* MAGs are drawn uniformly within each phylum with `random.Random(20261019)` from the
  `mag_id`-sorted pool.
* The 18 phyla with < 10 eligible MAGs (34 MAGs in total) are not sampled.

| Phylum | n | Phylum | n |
| --- | --- | --- | --- |
| Actinomycetota | 17 | Bacteroidota_A | 5 |
| Pseudomonadota | 12 | Acidobacteriota | 5 |
| Bacteroidota | 6 | Verrucomicrobiota | 4 |
| Bacillota | 6 | Deinococcota | 4 |
| Planctomycetota | 5 | Bdellovibrionota | 4 |
| Patescibacteriota | 5 | Thermoproteota (Archaea) | 4 |
| Myxococcota | 5 | Thermoplasmatota (Archaea) | 4 |
| Gemmatimonadota | 5 | Halobacteriota (Archaea) | 4 |
| Chloroflexota | 5 | **total** | **100 (17 phyla; 88 bacterial, 12 archaeal)** |

`mag_list.tsv` in this directory is the exact list. It has the full GTDB string, CheckM2
values, the asserted taxon (`NCBITaxon_2` for Bacteria, `NCBITaxon_2157` for Archaea), the
essential profile, and every input path.

## Per-MAG pipeline (`01_mag.sbatch`; 4 CPUs, 24 GB, ≤ 6 h)

1. **Inputs are checksummed inside the job** (`work/<mag>/input_sha256.txt`).
2. **Prep** (`mag_pathway_gaps.py prep`):
   * **Layout and operons.** The genome layout (`protein_id contig start end strand`) and the
     operons (same strand, intergenic gap ≤ 300 bp) come from the Bakta GFF3. Both use the
     `ID=`/locus_tag key, which equals the `.faa` id.
     `benchmark/make_operons.py` is not used because it keys on `Name=`, which in Bakta is the
     product name.
   * **Gap reactions.** These are the SEED reactions in gapseq's gap-filled model
     (`<mag>.xml`) that are absent from the draft (`<mag>-draft.xml`), i.e. the reactions
     gapseq had to add to reach growth on its predicted medium.
     Each gap reaction is annotated with:
     * its EC numbers (SBML `ec-code`),
     * GO terms via `ec2go`,
     * MetaCyc pathways (from `Reactions.tbl` `dbhit`),
     * gapseq's GPR genes (`G_<contig>_<start>_<end>`), mapped onto Bakta loci by contig
       sequence identity and coordinate overlap,
     * gapseq's own sub-threshold tblastn hits for the reaction (`Reactions.tbl`
       `stitle/sstart/send/status`), mapped the same way.
   * **Gaps file.** Gaps with ≥ 1 complete EC ("fillable") become `--gaps` JSONL records,
     one per (reaction, EC, GO), with `gapseq_guessed` set to whether gapseq attached a gene.
3. **Stage 0 claims (two channels).**
   * DIAMOND 2.1.16 `blastp --evalue 1e-5 --max-target-seqs 5` against the experimental-GO
     Swiss-Prot DB built by P13 (`bgc_repro_2026-10-19/ref/sprot_exp.dmnd`: 83,592 Swiss-Prot
     entries with ≥ 1 experimental GO annotation, GO release 2026-03-25). It is copied into
     `ref/`.
   * HMMER 3.3.2 `hmmsearch --cut_ga` against Pfam-A (27,481 models) with `pfam2go`.
   * Claims are built with `bgc_benchmark.py build-claims`, the same code and calibration as P13.
   * **InterProScan was not run.** Ibex has an install (`module interproscan/5.61-93.0`; also
     5.71/5.77 under `rob/gspa/databases/interproscan/`). However, the per-user CPU quota was
     saturated (`QOSMaxCpuPerUserLimit` on every queued job) when this was set up. Stage 0 is
     therefore DIAMOND + Pfam only.
4. **GSPA v2.1.0** (`allthebacteria/bin/gspa-unified.jar`, sha256 `e415d9c6…f72a7`; Java 21.0.5):
   * **Stage 0:** `integrate --lite --enable-priors none` (noisy-OR; ablation A1).
   * **Stage 1+2A + RLGC:** `integrate --lite --enable-priors
     essentiality,coherence,consistency,anatomical,gap_filling,genomic_context
     --essential-profile {bacteria|archaea} --taxon {NCBITaxon_2|NCBITaxon_2157} --pathways
     kegg_pathways.tsv --modules kegg_modules.tsv --ec2go ec2go.txt --operons <mag>_operons.tsv
     --gaps <mag>_gaps.jsonl --rlc-suggester --reaction-graph seed_reactions_corrected.tsv
     --diffusion-mets diffusion_mets.tsv --reaction-ec-aliases
     seed_Enzyme_Class_Reactions_Aliases_unique.tsv --genome-layout <mag>_layout.tsv
     --suggestions-out … --enforce-stage2 --enforce-consistency-mode relax` (A9 + gaps + RLGC).
   * The reaction-graph files are gapseq 1.4.0's own `dat/` tables, from the same conda env
     that built the GEMs. RLGC parameters are the defaults: k = 2, bandwidth 5 kb, anchor ≥ 0.3.
   * `--lite` (no ELK) matches P12/P13 and the ablation runs.
5. **Score** (`mag_pathway_gaps.py score`; the docstring has the full definitions):
   * **`candidate_fills_q030` / `hq_fills_q050`.** These count (gap reaction, Bakta locus) pairs
     where the locus carries a gap function after Stage 1+2A (posterior ≥ q), or where RLGC
     suggests the locus with q_value ≥ q.
     * A gap function is a complete EC, its ec2go GO term, or a GO descendant of that term.
       Generic catalytic roots never count.
     * A pair only counts if its Stage 0 posterior was < 0.30 for every gap function, so the
       0.50 fills are a subset of the 0.30 fills.
   * **Promotions.** Pairs that move from [0.30, 0.50) to ≥ 0.50 are counted separately as
     `promoted_fills_q050`.
   * **Pathway completeness.** For each gapseq-predicted MetaCyc pathway (`Prediction == true`),
     this is the fraction of its EC-bearing reactions that have a locus with posterior ≥ 0.5 for
     one of the reaction's functions; the MAG value is the mean over pathways.
     "Before" uses Stage 0. "After" uses Stage 1+2A plus RLGC suggestions with q ≥ 0.5.
     gapseq's own mean pathway completeness is reported for reference.
   * **Essential completeness.** The fraction of the jar's bacteria/archaea essential GO terms
     covered at posterior ≥ 0.5, including via descendants.
   * **Control** (what gapseq itself had for each fill's reaction):
     * a GPR gene in the gap-filled model at the same locus (`gapseq_gene_same_locus`) or at
       another locus (`gapseq_gene_other_locus`);
     * no GPR (gapseq only added the reaction) but gapseq's own sub-threshold tblastn hit at the
       same locus (`gapseq_subthr_hit_same_locus`) or at another locus
       (`gapseq_subthr_hit_other_locus`);
     * no gapseq hit at all (`gapseq_no_hit`).
   * **RLGC diagnostics.** `rlgc_max_q` is the highest RLGC q over the MAG's suggestions.
     `rlgc_top1_is_gapseq_hit_locus` counts gaps where RLGC's top locus equals gapseq's hit locus.
6. **Aggregate** (`02_aggregate.sbatch`): per-phylum sums and means plus an `ALL` row; all fills
   go to `results/mag_gap_fills_v2.tsv.gz`, timings to `results/timing_v2.tsv`, and per-MAG
   input checksums to `results/input_sha256_per_mag_v2.tsv`.

## Tool versions and references (`ref/tool_versions.txt`, `ref/sha256.txt`)

| Item | Version / source |
| --- | --- |
| GSPA | 2.1.0 unified jar |
| Java | OpenJDK 21.0.5 |
| DIAMOND | 2.1.16 |
| HMMER | 3.3.2 |
| Python | 3.11.0 |
| gapseq | 1.4.0 (GEMs and `dat/`) |
| Bakta | 1.x (DB 6.0) |
| GTDB-Tk | r226 |
| CheckM2 | as in `checkm2_single_sample_full.tsv` |
| GO | releases/2026-03-25 (`go.obo`, `go.owl`, `go-dag.tsv`) |
| Pfam-A | 27,481 models |

All reference files are snapshotted into `$OUT/ref/` and checksummed there.

## Job IDs

| Step | Job | State |
| --- | --- | --- |
| 00 select | 53996209 | COMPLETED (19 s) |
| test (3 MAGs, debug) | 53996226 (array 0-2), aggregate 53996227 | COMPLETED (~6 min per MAG) |
| test re-score after the scoring change | 53996792, aggregate 53996834 | COMPLETED |
| full array on batch | 53996870, aggregate 53996873 | cancelled: never started after ~45 min (`QOSMaxCpuPerUserLimit`) |
| full array on debug (`PART=debug`, 100 MAGs, %10) | 54001102 | COMPLETED: 100/100 DONE, 13:26–14:31 |
| full aggregate (debug) | 54003848 (54001103 cancelled) | COMPLETED |

## Result (100 MAGs, 17 phyla)

| | ALL |
| --- | --- |
| fillable gap reactions per MAG | 92.6 |
| new fills at posterior ≥ 0.30 / ≥ 0.50 | 1,272 / 1,258 (12.7 / 12.6 per MAG) |
| additional promotions [0.30, 0.50) → ≥ 0.50 | 842 |
| fills from the gap_filling prior / from RLGC (q ≥ 0.30) | 1,254 / 18 |
| MetaCyc pathway completeness, Stage 0 → Stage 1+2A | 0.140 → 0.274 (gapseq's own: 0.945) |
| essential completeness, Stage 0 → Stage 1+2A | 0.782 → 0.932 |
| control at q ≥ 0.50: gapseq GPR gene (same / other locus) | 0 / 0 |
| control at q ≥ 0.50: gapseq sub-threshold hit at the same locus / other locus / no hit | 742 / 491 / 25 |

Of the fills at q ≥ 0.50, 59 % land on the locus that gapseq's own sub-threshold tblastn hit
already pointed to. 39 % name a different locus, and 2 % fall on reactions where gapseq had no
hit at all. gapseq never gives a gene-level assignment for its gap-filled reactions in these
models. The per-phylum table is `results/pillar3_dark_matter_mags_v2.tsv`.

## Caveats

* Only two Stage 0 channels were used (no InterProScan); see step 3.
* None of the gapseq gap-filled reactions carries a GPR gene in these models. The meaningful
  gapseq-side control is therefore its sub-threshold homology hit, which is present for most
  fillable gaps.
* RLGC softmax q values are flat: the maximum over all 100 MAGs is 0.61. RLGC therefore
  contributes only 18 of the 1,272 fills at q ≥ 0.30, and its top-ranked locus equals gapseq's
  hit locus for 329 of 8,449 gaps (3.9 %). Almost every fill comes from the gap_filling prior
  acting on existing Stage 0 claims.
* The old 97-MAG table lists phyla (Spirochaetota, Campylobacterota, Asgardarchaeota,
  Thermodesulfobacteriota) and counts (8 Cyanobacteriota HQ MAGs) that do not occur among the
  3,794 eligible Empty Quarter HQ MAGs, which have 5 Cyanobacteriota. It cannot have been
  computed on this MAG set.
