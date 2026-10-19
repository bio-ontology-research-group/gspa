# P13: GSPA biosynthetic-gene-cluster benchmark (antiSMASH 8 regions, 108 genomes)

Manuscript element: the BGC specificity table (Suppl. Table S10). It replaces
`results/pillar2_3_antismash_106bgc_specificity.tsv` ("106 BGCs from 12
species"). No code, raw data or job record for that table could be found on
any host, so it cannot be reproduced. The new result is
`results/pillar2_3_antismash_bgc_specificity_v2.tsv`. The old file is left as
it was.

Compute host: ibex (SLURM, account `pi-hohndor`). Nothing runs on the login node.
Output root: `OUT=/ibex/scratch/projects/c2014/hohndor/bgc_repro_2026-10-19/`.
The METEOR antiSMASH tree (`/ibex/scratch/projects/c2014/kexin/funcarve/…`) is only read.

## Contents

| File | Role |
| --- | --- |
| `run.sh` | `sync` (laptop → ibex), `ref`, `test`, `all`, `integrate`, `fetch` |
| `jobs/common.sh` | all paths, tool versions, InterProScan application list |
| `jobs/00_ref.sbatch` | genome list, experimental-GO Swiss-Prot DIAMOND DB, class-term check, sha256 |
| `jobs/01_genome.sbatch` | per genome: antiSMASH GenBank → proteins / GFF / labels → operons → DIAMOND, HMMER, InterProScan → claims |
| `jobs/02_integrate.sbatch` | `gspa-unified.jar integrate --batch-manifest` for one configuration |
| `jobs/03_evaluate.sbatch` | result table with bootstrap CIs, per-genome audit counts, region inventory |
| `../../../benchmark/unified/pillar2_fba_lro_price149/bgc_benchmark.py` | all parsing, claim building, manifests and evaluation (Python standard library only) |
| `../../../benchmark/make_operons.py` | operons from the GFF (same strand, intergenic gap ≤ 300 bp) |

## Inputs

| Input | Path (ibex) | Version / note |
| --- | --- | --- |
| antiSMASH output, 108 genomes | `kexin/funcarve/paperA_2026/results/antismash_gc/<acc>/<acc>.gbk` | antiSMASH **8.0.4** (JSON `version`), schema 4, taxon bacteria; genes called by prodigal inside antiSMASH (`ctgN_M` locus tags); run without smCOG gene functions, so `gene_kind` is only `biosynthetic` or `biosynthetic-additional` |
| Genome list + taxonomy | `kexin/funcarve/paperA_2026/results/gcf_taxonomy_gc.tsv` | 108 = all antiSMASH directories; same set as METEOR `cohorts/panel108_gram.tsv` |
| Swiss-Prot | `rob/gspa-neural-deploy/reference/uniprot_sprot.dat.gz` (2026-04-22) | 574,627 entries; 83,592 have ≥1 experimental GO annotation (EXP, IDA, IPI, IMP, IGI, IEP, TAS, IC, HTP, HDA, HMP, HGI, HEP); 501,980 (accession, GO) pairs |
| Pfam-A | `rob/gspa/databases/pfam/Pfam-A.hmm` | 27,481 models; `hmmsearch --cut_ga` (HMMER 3.3.2) |
| pfam2go, KEGG pathways/modules, ec2go, go.owl | `hohndor/allthebacteria/ref/` | pfam2go 2026/03/02; go.owl releases/2026-03-25 |
| go.obo | `rob/gspa-neural-deploy/reference/go.obo` | releases/2026-03-25, the same release as go.owl |
| InterProScan | ibex module `interproscan/5.61-93.0` | apps Pfam-35.0, TIGRFAM-15.0, Hamap-2021_04, PIRSF-3.10, SMART-9.0, CDD-3.20, ProSiteProfiles/Patterns-2022_05, PRINTS-42.0, SFLD-4; `-goterms -iprlookup -dp`. Gene3D, SUPERFAMILY, PANTHER and FunFam were left out to save compute: in the test they added few matches that carry InterPro2GO terms |
| GSPA jar | `hohndor/allthebacteria/bin/gspa-unified.jar` | v2.1.0, sha256 `e415d9c6e13fbf96382ba7d798074b2ca4ec2991499de6714a5e4da2b19f72a7` (the "current jar" of P02) |
| Java | `rob/jdk-21/bin/java` | |
| DIAMOND | module `diamond/2.1.16` | `blastp --evalue 1e-5 --max-target-seqs 5` (default sensitivity) |

The sha256 sums of every reference are in `outputs/sha256.txt` (copied from `$OUT/ref/sha256.txt`).

## Benchmark design

**Loci.** Each antiSMASH region is one BGC. CDS features (with their
translations) are read from the full antiSMASH GenBank, so protein IDs,
coordinates and labels all come from one file. Each CDS gets one label:

* `core`: in a region, `gene_kind=biosynthetic`
* `tailoring`: in a region, `gene_kind` ∈ {`biosynthetic-additional`, `transport`, `regulatory`}. Only `biosynthetic-additional` occurs in these runs.
* `region_other`: in a region with no gene kind. These are excluded from both positives and negatives.
* `nonbgc`: outside every region. These are the negatives.

**Classes and labels.** Each region's class comes from its antiSMASH
protocluster `category`. Siderophores are the exception and are identified by
protocluster `product`. A core gene gets the class of the antiSMASH rule that
fired on it (`gene_functions`, e.g. `biosynthetic (rule-based-clusters) T1PKS: PKS_AT`).
A tailoring gene gets the classes of the protoclusters that span it. In hybrid
regions, this keeps a PKS core gene from being scored against the NRPS term.

| Class | antiSMASH source | GO term (checked against go.obo 2026-03-25) |
| --- | --- | --- |
| NRPS | category `NRPS` (NRPS, NRPS-like, NRP-metallophore, thioamide-NRP, NAPAA, isocyanide-nrp …) | GO:0019184 nonribosomal peptide biosynthetic process |
| PKS | category `PKS` (T1PKS, T2PKS, T3PKS, transAT-PKS, PKS-like, HR-T2PKS, arylpolyene …) | GO:0030639 polyketide biosynthetic process |
| terpene | category `terpene` (terpene, terpene-precursor) | GO:0016114 terpenoid biosynthetic process |
| siderophore | product `NI-siderophore` or `NRP-metallophore` | GO:0019290 siderophore biosynthetic process |
| RiPP | category `RiPP` (lanthipeptides, lassopeptide, RiPP-like, RRE-containing, thiopeptide …) | GO:0030152 bacteriocin biosynthetic process |
| ALL_BGC | every region, including `saccharide` and `other` | GO:0044550 secondary metabolite biosynthetic process **or** any of the five class terms. GO places GO:0019184 outside GO:0044550, so the generic term alone would miss NRPS genes |

A protein receives a term when the highest posterior over that term and its
`is_a`/`part_of` descendants is ≥ 0.50. Recall is pooled over the target genes
of all genomes. Specificity is 1 − FP/N over the non-BGC genes of all 108
genomes. Precision is TP/(TP+FP), where TP covers core and tailoring genes and
FP covers non-BGC genes. F1 combines precision with the pooled core+tailoring
recall. The 95 % CIs are percentile bootstrap intervals over genomes (all
genomes resampled, 1000 replicates, seed 20261019). ECs are not evaluated: the
channels emit GO only.

**Stage 0 claims.** Three classical channels are combined into one JSONL per
genome. Each channel keeps the maximum score per (protein, GO) pair. Aspects
come from go.obo, and obsolete terms are dropped.
1. `diamond`: DIAMOND against Swiss-Prot entries with experimental GO. raw = pident/100.
2. `pfam`: hmmsearch Pfam-A `--cut_ga` + pfam2go. raw = domain bitscore/100.
3. `interproscan`: InterPro2GO. raw = 0.95/(1+exp(−0.1(−log10 E − 10))), clipped to [0.1, 0.95]. This is the same calibration as `benchmark/02b_parse_predictors_to_claims.py`.

Operons come from the CDS coordinates via `benchmark/make_operons.py` (gap ≤ 300 bp).
No gapseq gaps are supplied, so the `gap_filling` prior in A9 has no input.

**Configurations** (one `--batch-manifest` run each; shared flags
`--go-owl ref/go.owl --pathways kegg_pathways.tsv --modules kegg_modules.tsv --ec2go ec2go.txt --essential-profile bacteria --taxon bacteria --lite`):

| Name | Manifest `enable_priors` | Stage 2 | Operons |
| --- | --- | --- | --- |
| `stage0_noisy_or` | `none` | no | no |
| `operon_only` | `genomic_context` | no | yes |
| `gspa_A9` | `essentiality,coherence,consistency,anatomical,gap_filling,genomic_context` | `--enforce-stage2`, mode `relax` | yes |

The asserted taxon is `bacteria` (NCBITaxon_2) for every genome. The strain-level
taxids in `gcf_taxonomy_gc.tsv` are not in GSPA's bundled taxon hierarchy, so the
SAT checker could not place them.

## Running

```bash
bash paper-reproducibility/pipelines/P13_bgc_antismash/run.sh sync
ssh ibex 'bash /ibex/scratch/projects/c2014/hohndor/bgc_repro_2026-10-19/code/run.sh ref'
ssh ibex 'bash /ibex/scratch/projects/c2014/hohndor/bgc_repro_2026-10-19/code/run.sh test'   # optional
ssh ibex 'bash /ibex/scratch/projects/c2014/hohndor/bgc_repro_2026-10-19/code/run.sh all'
bash paper-reproducibility/pipelines/P13_bgc_antismash/run.sh fetch
```

`all` submits 27 array tasks with 4 genomes each (`PER_TASK=4`). This keeps the
job count under ibex's per-user limit of 2000 queued jobs. Integration then
runs with `afterany`, and evaluation with `afterok`.

## Jobs of the run reported here (2026-10-19)

| Step | Job ID(s) | Partition | Wall | Note |
| --- | --- | --- | --- | --- |
| `00_ref` | 53986260 | debug | 1 min 15 s | first submission 53986254 (batch) cancelled while pending, resubmitted to debug |
| `test`: genomes / integrate ×3 / eval | 53986551 (array 0-2) / 53986552, 53986553, 53986554 / 53986555 | debug | 7–9 min / 17 s–65 s / 3 s | 3 genomes, InterProScan on the first 200 proteins; outputs in `$OUT/test/` and `results/*_test.tsv`. The test's InterProScan/claims files were moved to `$OUT/test/work_ips200/` before the full run |
| `01_genome` | 53987114 (array 0-26, 4 genomes/task) | batch | 1 h 20 min – 2 h 42 min per task | submitted as 16 CPU/48 GB; changed while pending to 8 CPU/24 GB (`scontrol update`) to get past the per-user CPU cap. Task 0 failed on its three test genomes (shell precedence bug `[[ -s f ]] \|\| cmd && mv`, now fixed) |
| `01_genome` rerun | 54002128 (array 0) | batch | 1 h 04 min | the fixed script on task 0 (GCA_003598085.1, GCA_014873095.1, GCA_022305115.1, GCA_046479585.1) |
| `02_integrate` | 54025653 (stage0_noisy_or), 54025654 (operon_only), 54025656 (gspa_A9) | batch | 31 s / 1 min 40 s / 8 min 12 s | 108/108 genomes each. The first submission (53987115-7, 160 GB) was cancelled while pending and resubmitted at 64 GB, `-Xmx56g` |
| `03_evaluate` | 54025658 | batch | 18 s | `results/pillar2_3_antismash_bgc_specificity_v2.tsv` |
| sensitivity (≥ 0.30) | 54026488 | debug | 32 s | `outputs/pillar2_3_antismash_bgc_specificity_v2_thr030.tsv` |
| posterior diagnostic | 54026509 | debug | 23 s | `jobs/diag_posteriors.py` → `outputs/diag_posteriors_54026509.txt` |

Compute: 518,063 proteins. The channel steps took 416 CPU-h at 8 threads per
genome: InterProScan 160,742 s wall, HMMER 26,102 s and DIAMOND 259 s, summed
over genomes. SLURM allocated 16 logical CPUs per 8-CPU task, so the reserved
allocation was about 860 core-h. Integration and evaluation are negligible.

## Result (posterior ≥ 0.50; 95 % bootstrap CI over the 108 genomes)

| Class | Genomes | Clusters | Core / tailoring loci | Stage 0 core recall | GSPA A9 core recall | Operon-only core recall | Stage 0 tail recall | GSPA A9 tail recall | GSPA specificity | Stage 0 F1 | GSPA F1 | Any-claim core / tail |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| NRPS (GO:0019184) | 72 | 359 | 778 / 1877 | 0.069 (0.045–0.104) | 0.069 (0.045–0.104) | 0.069 (0.045–0.104) | 0.019 (0.011–0.029) | 0.019 (0.011–0.029) | 0.9996 (0.9995–0.9996) | 0.061 (0.040–0.090) | 0.061 (0.040–0.090) | 0.713 / 0.093 |
| PKS (GO:0030639) | 75 | 270 | 649 / 1645 | 0.003 (0.000–0.008) | 0.003 (0.000–0.008) | 0.003 (0.000–0.008) | 0.002 (0.000–0.004) | 0.002 (0.000–0.004) | 0.9999 (0.9999–1.0000) | 0.004 (0.001–0.009) | 0.004 (0.001–0.009) | 0.542 / 0.089 |
| terpene (GO:0016114) | 106 | 303 | 426 / 649 | 0.157 (0.123–0.192) | 0.157 (0.123–0.192) | 0.157 (0.123–0.192) | 0.133 (0.107–0.168) | 0.133 (0.107–0.168) | 0.9993 (0.9992–0.9994) | 0.196 (0.169–0.222) | 0.196 (0.170–0.222) | 0.575 / 0.185 |
| siderophore (GO:0019290) | 63 | 117 | 305 / 462 | 0.528 (0.448–0.627) | 0.528 (0.448–0.627) | 0.528 (0.448–0.627) | 0.078 (0.050–0.111) | 0.078 (0.050–0.111) | 0.9998 (0.9998–0.9999) | 0.377 (0.330–0.433) | 0.377 (0.330–0.433) | 0.764 / 0.162 |
| RiPP (GO:0030152) | 93 | 314 | 570 / 554 | 0.000 (0.000–0.000) | 0.000 (0.000–0.000) | 0.000 (0.000–0.000) | 0.000 (0.000–0.000) | 0.000 (0.000–0.000) | 1.0000 (1.0000–1.0000) | 0.000 (0.000–0.000) | 0.000 (0.000–0.000) | 0.000 / 0.000 |
| ALL_BGC (GO:0044550 or any class term) | 106 | 1292 | 2778 / 4683 | 0.095 (0.079–0.113) | 0.093 (0.077–0.111) | 0.095 (0.079–0.113) | 0.028 (0.022–0.036) | 0.028 (0.022–0.036) | 0.9986 (0.9985–0.9988) | 0.093 (0.080–0.109) | 0.092 (0.079–0.108) | 0.503 / 0.188 |

* GSPA Stage 1+2A and the operon-only prior leave the result unchanged. The
  only differences in the third decimal are terpene precision and the ALL_BGC
  core recall (0.0947 → 0.0932). In the posterior diagnostic, 1,385 (core
  gene, class term) pairs have at least one claim, and A9 moves 19 of them by
  more than 0.01. The `genomic_context` prior only boosts GO terms that a KEGG
  pathway in `kegg_pathways.tsv` requires, and none of the five class terms
  occurs there. Synteny therefore cannot carry a class term from a core gene to
  its tailoring neighbours.
* The bottleneck is claim strength, not claim coverage. 50–76 % of core genes
  have a claim for their class term (none for RiPP), but the median Stage 0
  posterior of those claims is 0.24. Of the 2,568 class-term claims on core
  genes, 2,255 come from DIAMOND hits to distant experimentally annotated
  Swiss-Prot proteins, and pident/100 calibrates them below 0.5.
* RiPP gets zero recall: no channel produces GO:0030152 (or a descendant) for
  any RiPP core or tailoring gene.
* Specificity is ≥ 0.9986 for every class and configuration (FP ≤ 0.14 % of 480,236 non-BGC genes).
* At posterior ≥ 0.30 (`outputs/…_thr030.tsv`) recall rises, e.g. core recall
  for PKS 0.18, terpene 0.29, siderophore 0.57 and ALL_BGC 0.17. Stage 0 and
  GSPA remain within 0.003 of each other.

## Comparison with the old table

`results/pillar2_3_antismash_106bgc_specificity.tsv` reported "106 genomes,
106 BGCs" in five classes. It gave Stage 0 core recall of 0.72–0.83, rising to
0.89–0.97 under GSPA, tailoring recall of 0.48–0.60 rising to 0.80–0.89,
specificity of 0.9955–0.9972 and F1 of 0.68–0.79 rising to 0.86–0.94. No
generating code, inputs or job record exist for it. Its counts do not match
any antiSMASH output on ibex: these 108 genomes have 1,292 regions, not 106.
This real experiment does not reproduce its numbers. Classical
homology/domain claims recover the class term for 0–53 % of core genes
(0 % RiPP, 0.3 % PKS, 7 % NRPS, 16 % terpene, 53 % siderophore). GSPA's
Stage 1/2A priors do not raise that, and the old table's large GSPA gains are
not supported. The specificity claim (> 0.995) is supported, and is in fact
higher here (≥ 0.9986).


See `OUTPUTS.md` for every output file.
