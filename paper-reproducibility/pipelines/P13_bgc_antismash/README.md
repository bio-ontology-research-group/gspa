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

JOB_TABLE_PLACEHOLDER

See `OUTPUTS.md` for every output file.
