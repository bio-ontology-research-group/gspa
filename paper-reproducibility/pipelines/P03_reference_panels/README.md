# P03: reference panels (Fig. 2d, Suppl. Table S2, former Table 2 external block)

This pipeline replaces the hand-transcribed `results/pillar1_suite2_panel13_21_sota.tsv` with values computed from the raw outputs that exist on unimatrix01. All methods are scored by one evaluator under one protocol.

```bash
./run.sh            # sync code, sbatch both panels on unimatrix01 (node007), wait, fetch, aggregate
./run.sh aggregate  # fetch existing remote outputs and rebuild results/*.tsv only
```

## Panels

**13-genome panel.** This is the panel scored for metagenomic-deepFRI (`/data/hohndor/mdf-runs/scoring/`). It contains:
- the 10 bench10 genomes: vcholerae, saureus, spneumoniae, ccrescentus, rprowazekii, tpallidum, tthermophilus, dradiodurans, scoelicolor and pfuriosus (`/data/hohndor/gspa/proteomes/bench10/`);
- 3 bench9 genomes: hpylori, mgenitalium and mjannaschii (`/data/hohndor/gspa/proteomes/{bench9,truth_dual}/`).

All 13 genomes are prokaryotes, 11 bacteria and 2 archaea, and every one of them has an NCBI PGAP annotation: 10 in `bench10/pgap/` and 3 in `proteomes/{tag}_pgap.tsv`. The manuscript calls this "11 prokaryotes = 13 minus 2 eukaryotes", but no eukaryote and no genome without PGAP exists in the inputs, so that split could not be reproduced. All tables report the 13 genomes (`Panel13_Prok13`).

**21-genome panel.** This is `/data/hohndor/gspa-neural/panel/` (21 bacteria, RefSeq-keyed truth). Only the methods with raw outputs are evaluated. See the `notes` column of the result table.

## Methods (13-genome panel)

| id | Source of predictions | Version |
| --- | --- | --- |
| pgap | GO qualifiers in the NCBI PGAP GFF, extracted with `benchmark/extract_pgap_annotations.py` and mapped from RefSeq to UniProt; binary score of 1.0 | RefSeq GFFs as downloaded in `bench10/step1_download.sh` |
| pfam | `source=pfam` claims (`hmmsearch --domtblout -E 1e-5` against Pfam-A, then pfam2go) | HMMER 3.4 (Aug 2023); `/storage/software/databases/hmmer/Pfam-A.hmm` (27,481 families; HMMER3/f 3.3 header; release number not recorded in the file); pfam2go dated 2026/03/02 |
| interproscan | `source=interproscan` claims (`parse_interproscan_claims.py`, InterPro2GO, evalue-derived score) | InterProScan 5.72-103.0 (`/data/hohndor/envs/interproscan-5.72-103.0`, per the job logs) |
| diamond | `source=diamond` claims: best-hit GO transfer from Swiss-Prot with the panel proteins left out (`reference_loo19.dmnd` for bench10, `reference_loo9.dmnd` for bench9) | as in the bench10/bench9 predictor runs |
| mdf | `/data/hohndor/mdf-runs/{tag}/results.tsv`, maximum score per (protein, GO term) | metagenomic-deepFRI, model `pretrained_deepfri_wg20697.pkl` |
| gspa3 | `gspa integrate` on the 3-channel claims (`bench10/claims/{tag}_dpi_merged.jsonl`, `bench9/with_iprscan/{tag}_claims_merged.jsonl`) with the A9_stage2_relax flags of `benchmark/unified/pillar1_2/01_run_ablation_grid.py` | `/data/hohndor/gspa/bin/gspa-unified.jar` (gspa 2.1.0); GO `releases/2026-03-25` |

Notes on the gspa3 inputs:
- **Taxon.** The species-level taxon is used where the bundled hierarchy (`gspa-core/.../taxon-constraints/ncbi-taxon-hierarchy.tsv`) contains it. Otherwise the deepest ancestor it contains is used: Alphaproteobacteria for ccrescentus and rprowazekii, Spirochaetia for tpallidum, Deinococci for tthermophilus, and Thermococcaceae for pfuriosus. See `p03_common.py`.
- **Operons.** These are UniProt-keyed operon files: `bench10/ops/*_operons_up.tsv` and `proteomes/operons/`.
- **Gaps.** Only real gapseq gaps are used: `bench9/gapseq/` and `bench10/gaps/*_real.jsonl` for 4 genomes. The `bench10/gaps/*_gaps.jsonl` files are synthesised from integrated output and are therefore not used.

## Protocol (`p03_eval.py`)

The metric definitions are those of `benchmark/benchmark_pgap_v2.py`, as used by `benchmark/unified/pillar1_2/02_run_cafa_eval.py`: micro F-max, protein-centric CAFA F-max, coverage and IC recall at 0.5, and S-min. That evaluator applies no ancestor propagation, but the paper's protocol requires it. The P03 evaluator therefore adds the following:
- Truth and predictions are propagated over `is_a` and `part_of` in `go.obo`.
- Alternative IDs are mapped to their primary IDs, and obsolete terms are dropped.
- A propagated term takes the maximum score of its descendants.
- The three root terms are removed.
- The aspect of each term comes from its GO namespace.
- Thresholds run from 0.01 to 0.99 in steps of 0.01, and prediction rows below 0.01 are dropped.
- IC is -log2 of a term's frequency in the genome's propagated all-evidence truth. The same IC values are used for experimental-only truth.

**Caveat on all-evidence truth.** Roughly 85 % of the GOA annotations behind the truth files are IEA: 668,604 of the rows in `bench10/goa_subset.gaf`, 219,897 of them assigned by InterPro. All-evidence F-max therefore rewards methods that reproduce GOA's electronic pipelines, which includes InterProScan itself and every GSPA configuration that has an InterProScan channel. Experimental-only truth has no IEA, but it is small: 323 proteins across the 13 genomes and 228 across the 21.

## Structural audit (`p03_audit.py`)

This step uses `make_synth_gff_gaf.py --threshold 0.5` and then `gspa evaluate --lite` with the KEGG pathway, module and ec2go files, run as a batch manifest. The settings are identical to `03_run_gaef_eval.py`. This is variant `t05`.

A fixed cut-off of 0.5 is arbitrary for uncalibrated external scores. Variant `fmaxt` therefore applies each method's per-genome micro F-max threshold, taken from the all-evidence evaluation. Both variants are reported.

## Perturbation (`p03_build.py`, `p03_perturb_eval.py`)

This is design v2 of `benchmark/unified/pillar1_2/04_run_perturbation.py`:
- **Injected terms.** 100 cross-clade terms are chosen, with the same terms and seed: 50 anatomical (5 mandatory) and 50 restricted to eukaryotes by taxon constraints.
- **Proteins.** Each term goes to one protein, sampled with `seed = md5(tag)`.
- **Scores.** Raw scores are drawn from U(0.62, 0.90).

**Protein-ID fix for the bench9 genomes.** In `bench9/with_iprscan/{tag}_claims_merged.jsonl`, every pfam claim carries the raw FASTA id (`sp|P21762|AHPC_HELPY`); there are 11,325, 5,739 and 13,773 such claims for hpylori, mgenitalium and mjannaschii. `p03_build.normalize_claims()` strips these to the accession and drops 27 and 45 exact duplicates (hpylori, mjannaschii).

Without this fix the code reproduces the original injected pairs byte for byte, but 47, 42 and 55 of the 100 injected pairs then sit on these phantom ids. As a result, `04_run_perturbation.py` undercounts survivors for these genomes whenever `gspa integrate` normalises the ids. For hpylori A1 it counts 41 survivors where 77 is correct; A9 is unaffected (3 in both cases). With normalised ids, the injected pairs for those three genomes differ from the P02 files.

The GSPA-integrated F-max values did not change with the fix, which indicates that `gspa integrate` already normalises these ids internally.

Each method's output is used as Stage 0 claims (`source` = pgap / pfam / interproscan / diamond / mdf). Unknown sources get gspa's default calibration, sigmoid(4r-2). Each injected pair is added as one claim from that source. Survivors are counted in two ways:
- **(a) No GSPA:** the maximum raw score for the pair is at least 0.5.
- **(b) With GSPA:** `gspa integrate` A9 runs on the method's claims, and the pair survives if its posterior is at least 0.5.

The GSPA row injects claims exactly as `04_run_perturbation.py` does: a diamond or deepgo claim, plus a corroborating claim for 60 % of the pairs.

## 21-genome panel (`p03_panel21.py`)

Methods are evaluated against `gspa-neural/panel/truth/{tag}_truth_{all,exp}_refseq.tsv`, propagated by the evaluator.

| Method | Raw output used | Note |
| --- | --- | --- |
| esm2-deepgoplus, proteinfer, clean | `preds/{tag}/{tag}.<m>.tsv` | clean is EC only |
| diamond | `bench_gtdb30/claims/{tag}_dp_claims.jsonl` | UniProt IDs mapped to RefSeq via `panel/maps/{tag}_map.tsv` |
| gspa2 | the same D+P claims, run through A9 | The Pfam part of these claims is truncated (15 pairs per genome), so this is effectively DIAMOND-only |
| esm2t30-knn | `bench_gtdb30/plm/*_esm2t30.npy` | This is ESM-2 t30 (150M), not 650M. Cosine k = 10, similarity-weighted labels, leave-one-genome-out over the 30 bench_gtdb30 genomes with GOA truth. |

## Compute

- **Host.** unimatrix01 Slurm `debug` partition, on node007. node007 exports `/data` and sees it as `/srv/data`, so the jobs set `P03_DATA=/srv/data`.
- **Logs.** `~hohndor/p03logs/`.
- **Software.** Python `/storage/miniforge3/envs/metagenomics/bin/python3` (3.11, numpy 1.26.4) and Java `/data/hohndor/envs/java21`.

## Files

`p03_common.py` (configuration), `p03_build.py`, `p03_integrate.py`, `p03_audit.py`, `p03_eval.py` (evaluator), `p03_eval13.py`, `p03_perturb_eval.py`, `p03_panel21.py`, `p03_aggregate.py`, `p03_manifest.py` (writes the data_manifest rows, P03), `p03_panel13.sbatch`, `p03_panel21.sbatch`, `run.sh`, `OUTPUTS.md`.
