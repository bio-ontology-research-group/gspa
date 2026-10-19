# P13 outputs

All outputs are written under `OUT=/ibex/scratch/projects/c2014/hohndor/bgc_repro_2026-10-19/`
(the METEOR antiSMASH tree under `kexin/funcarve/` is only read).

| Path (ibex, under `$OUT`) | Produced by | Content |
| --- | --- | --- |
| `genomes.tsv` | `00_ref` | genome accession, NCBI taxid, organism (antiSMASH dirs ∩ `gcf_taxonomy_gc.tsv`) |
| `ref/sprot_exp.fasta`, `ref/sprot_exp_go.tsv`, `ref/sprot_exp.dmnd` | `00_ref` | Swiss-Prot entries with ≥1 experimental GO annotation, and their (accession, GO) pairs |
| `ref/class_terms_verified.tsv` | `00_ref` | class→GO map checked against `ref/go.obo` (releases/2026-03-25) |
| `ref/sha256.txt`, `ref/build_ref.log` | `00_ref` | checksums of jar and references; antiSMASH version; Pfam model count |
| `work/<genome>/<genome>.faa`, `.gff` | `01_genome` | proteins and CDS coordinates taken from the antiSMASH GenBank (prodigal calls, `ctgN_M` ids) |
| `work/<genome>/<genome>_genes.tsv` | `01_genome` | per-CDS label: `core` / `tailoring` / `region_other` / `nonbgc`, region id, benchmark classes |
| `work/<genome>/<genome>_regions.tsv` | `01_genome` | antiSMASH regions with protocluster categories, products, benchmark classes |
| `work/<genome>/<genome>_operons.tsv` | `01_genome` | `benchmark/make_operons.py` (same strand, gap ≤ 300 bp) |
| `work/<genome>/<genome>_{diamond.tsv,pfam.domtbl,ips.tsv}` | `01_genome` | raw channel outputs |
| `work/<genome>/<genome>_claims.jsonl` | `01_genome` | Stage 0 claims (sources `diamond`, `pfam`, `interproscan`) |
| `work/<genome>/<genome>_timing.tsv` | `01_genome` | wall-clock seconds per step |
| `manifests/manifest_<cfg>.tsv`, `manifests/batch_stats_<cfg>.tsv` | `02_integrate` | `gspa integrate --batch-manifest` input and per-genome runtime |
| `integrated/<cfg>/<genome>.tsv` | `02_integrate` | GSPA posteriors, `cfg` ∈ {`stage0_noisy_or`, `operon_only`, `gspa_A9`} |
| `results/pillar2_3_antismash_bgc_specificity_v2.tsv` | `03_evaluate` | the result table (one row per class) |
| `results/bgc_per_genome_counts_v2.tsv` | `03_evaluate` | per genome × config × class TP / N / FP counts (audit) |
| `results/bgc_regions_all_v2.tsv`, `results/bgc_timing_v2.tsv` | `03_evaluate` | region inventory, timing |
| `test/…`, `results/*_test.tsv` | `run.sh test` | 3-genome debug run (InterProScan on the first 200 proteins only) |

`run.sh fetch` copies the result table to `results/pillar2_3_antismash_bgc_specificity_v2.tsv`
in the repository and the audit files into `outputs/` of this pipeline directory.

## Result table columns

`bgc_class, go_id, go_name, n_genomes_total, n_genomes_with_class, n_clusters,
n_core_loci, n_tailoring_loci, n_nonbgc_genes, threshold`, then for each of
`stage0` (noisy-OR, no priors), `operon_only` (genomic_context prior only) and
`gspa` (A9: all Stage 1 priors + Stage 2A relax): `*_core_recall`,
`*_tail_recall`, `*_specificity`, `*_precision`, `*_f1`, each with `_lo`/`_hi`
(95 % percentile bootstrap over genomes, 1000 resamples, seed 20261019).
`claim_any_core_recall` / `claim_any_tail_recall` give the fraction of target
genes that have *any* Stage 0 claim for the class term (or a descendant) at any
score — the ceiling reachable by re-weighting existing claims.

## Additional outputs

| Path | Produced by | Content |
| --- | --- | --- |
| `results/pillar2_3_antismash_bgc_specificity_v2_thr030.tsv` (+ `bgc_per_genome_counts_v2_thr030.tsv`) | `run.sh sensitivity` (job 54026488) | same table at posterior ≥ 0.30 |
| `logs/diag-<job>.out` | `run.sh sensitivity` → `jobs/diag_posteriors.py` (job 54026509) | posterior quantiles of class-term claims on core genes per configuration; number of pairs changed by A9 |

Local copies in `outputs/` (from `run.sh fetch`): `genomes.tsv`, `sha256.txt`, `build_ref.log`,
`class_terms_verified.tsv`, `bgc_per_genome_counts_v2.tsv`, `bgc_regions_all_v2.tsv`, `bgc_timing_v2.tsv`,
`pillar2_3_antismash_bgc_specificity_v2_thr030.tsv`, `bgc_per_genome_counts_v2_thr030.tsv`, `diag_posteriors_54026509.txt`.
