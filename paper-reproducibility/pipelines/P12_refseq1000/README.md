# P12 — 1,000-genome RefSeq bacterial structural audit (Suppl. Table S8)

Manuscript element: Suppl. Table S8 (`results/pillar1_suite2_1000bacteria_sota.tsv`).
Compute: IBEX (Slurm, account `pi-hohndor`), plus two GOA scans on unimatrix01
(Slurm `debug`, head node and node007 excluded). No step runs on a login node.

The audit compares six published GO predictors and GSPA on 1,000 RefSeq
bacterial assemblies (`annot_evaluation/1000files`). The six tool rows come
from the predictions and derived matrices of the `annot_evaluation` project
(`/ibex/scratch/projects/c2014/annot_evaluation`, read-only, a colleague's
project). Before this package existed there were no GSPA outputs for these
genomes. The old Stage 0 and Stage 1+2 rows, and every CAFA, Smin, runtime and
RAM value in the old table, had no producer.

## What is computed, and by which code

| Table column | Six tools | GSPA rows | Code |
| --- | --- | --- | --- |
| essential core-22 %, full % (32 terms), complex coherence %, MetaCyc completeness %, has_part coherence % | recomputed from the per-genome prediction files | computed from the GSPA output at posterior ≥ 0.50 | `p12_structural_suite.py` (one evaluator for all 8 rows) |
| the same five metrics, from the audit's derived matrices | matrix reproduction (check only) | — | `p12_matrix_rows.py` |
| taxon-violating rows per genome (propagated; leaf; leaf implied; proteins) | suite | suite | `p12_structural_suite.py` |
| `gaef_*` columns: completeness, complex and KEGG-pathway coherence, SAT consistency violations | `gspa evaluate --lite` on leaf sets | same | `03_eval_chunk.sh` |
| runtime_sec_per_genome, peak_ram_mb | unavailable (the audit logged neither) | `/usr/bin/time -v` for each step | `02_gspa_genome.sh`, `05_aggregate.py` |
| cafa_macro_fmax, smin | unavailable | unavailable | see "CAFA truth" below |

### Metric definitions (`p12_structural_suite.py`)

These follow the audit's own scripts exactly. Each was checked against the
audit's per-genome matrices for all six tools on all 1,000 genomes.

* **Essential:** presence of the 22 core and 32 full essential GO terms
  anywhere in the genome's propagated annotation. The old table header said
  "full-29". The audit's list and matrices contain 32 terms.
* **MetaCyc:** a pathway is *annotated* when its GO term is present. It is
  *complete* when, in addition, every EC-derived GO term of at least one EC
  combination is present. The genome score is complete / annotated × 100,
  or 0 when no pathway is annotated. Inputs: `metacyc_GO_v2025-03-16_with_EC.tsv`
  and `ec2go_v2025-03-16`.
* **Complex coherence:** every is_a/part_of descendant of GO:0032991 carried
  by a protein that also carries GO:0032991 is a complex term. A term is
  coherent when two or more proteins carry it. No homodimer list is used,
  matching the audit's `J.sh`. The score is coherent / max(1, coherent + incoherent).
* **has_part coherence:** 100 − the percentage of has_part-required terms
  (`has-part-latest.txt`) that are missing from the genome.
* **Taxon violation:** a (protein, GO) row violates when a `never_in_taxon`
  taxon lies in the genome lineage, or when none of its `only_in_taxon` taxa
  does (Union classes are ignored). Inputs: `go_taxon_constraints_updated.tsv`
  and `propagated_taxon_ids.tsv`. Four counts are reported:
  * **Propagated rows:** rows of the `*_GO_ancestors.tsv` sets.
  * **Leaf rows that violate directly:** rows of the `*_GO_specific.tsv` sets.
  * **Leaf rows with a violating term in their ancestor closure.**
  * **Proteins with at least one violating row.**
* **Propagation:** GO 2025-03-16 `go-basic.obo`, re-implementing the audit's
  `utils.Ontology(with_rels=True)`. Every `relationship:` line counts as a
  parent edge.
* **Thresholds:** tools use the audit's thresholded files, which apply
  tool-specific cut-offs. Examples are DeepGOMeta (mf/cc 0.27, bp 0.11) and
  TALE (mf 0.28, cc 0.56, bp 0.15). InterProScan is unscored. Most tools'
  raw scores are not available (DeepGOMeta and InterProScan have none), so a
  common ≥ 0.50 cut-off can only be applied to GSPA (`posterior_prob ≥ 0.50`).

## GSPA run (`02_gspa_genome.sh`)

Each genome runs the same channels, flags and jar as ablation A1/A9 of
pillar1_2 (`benchmark/unified/pillar1_2/01_run_ablation_grid.py`):

| Item | Value |
| --- | --- |
| Proteins | `annot_evaluation/done_proteins/<g>_protein.faa` (RefSeq WP_ ids) |
| DIAMOND | v2.1.8 (IBEX module). pillar1_2 used v2.0.15, so the DB was rebuilt from the same FASTA. `blastp --evalue 1e-5 --max-target-seqs 50 --query-cover 50 --subject-cover 50 --id 30`, outfmt 6 as in `run_predictors_9.sh` |
| DIAMOND reference | `reference_loo9.fasta`: Swiss-Prot with the nine bench genomes removed (556,196 sequences), the DB behind the pillar1_2 claims. It is all of Swiss-Prot, not only experimentally annotated entries. |
| GO transfer | `02b_parse_predictors_to_claims.py` (byte-identical to the unimatrix01 copy that built bench9) with `goa_uniprot_all.gaf.gz` (GOA 2026-03-04, all evidence codes, NOT rows dropped), pre-filtered by `00_build_goa_subset.py` to the reference accessions |
| HMMER | 3.3.2, `hmmsearch --noali -E 1e-5`, Pfam-A.hmm from unimatrix01 (27,481 models, file dated 2026-01-13), `pfam2go.txt` (2026-03-02) |
| InterProScan | reused: `annot_evaluation/preds_IPscan/<g>_IPscan.tsv` (5.x run of 2024-06-30), parsed with `parse_interproscan_claims.py` |
| Operons | RefSeq GFF3 from NCBI FTP → `make_operons.py` (300 bp). 999 of 1,000: the GFF of GCF_020097475.1 returns 404, so that genome runs without operons |
| Gaps | **not used.** gapseq was not run for 1,000 genomes, so the gap_filling prior is enabled but has no gap claims |
| Organism taxon | deepest lineage taxon present in the jar's `ncbi-taxon-hierarchy.tsv` (`01_taxon_map.py`, `taxon_map.tsv`), mostly class or order |
| Stage 0 | `integrate --lite --enable-priors none` (= A1) |
| Stage 1+2A | `integrate --lite --enable-priors essentiality,coherence,consistency,anatomical,gap_filling,genomic_context --essential-profile bacteria --taxon T --pathways kegg_pathways.tsv --modules kegg_modules.tsv --ec2go ec2go.txt --enforce-stage2 --enforce-consistency-mode relax --operons …` (= A9) |
| Jar | `gspa-unified.jar` v2.1.0, sha256 `e415d9c6e13fbf96382ba7d798074b2ca4ec2991499de6714a5e4da2b19f72a7` (identical to the pillar1_2 jar) |
| Java | `/ibex/scratch/projects/c2014/rob/jdk-21/bin/java`, `-Xmx12g`, 4 cores per genome |

**Cost estimate**, from the four debug test genomes (job 53987583), with
4 cores per genome. Wall time, in seconds:

| Proteins | DIAMOND | hmmsearch | Claim parsing | Stage 0 integrate | Stage 1+2A integrate |
| --- | --- | --- | --- | --- | --- |
| 887 | 17 | 158 | 10 | 27 | 56 |
| 1,572 | 18 | 174 | 10 | 28 | 67 |
| 2,527 | 18 | 227 | 12 | 31 | 82 |
| 5,422 | 20 | 406 | 12 | 33 | 109 |

hmmsearch has a fixed cost of about 2 min for loading Pfam-A. Peak RSS is
3.1–4.3 GB, all of it in the `integrate` JVM. At the mean genome size of
3,197 proteins this comes to about 7 min of wall time on 4 cores, or about
0.5 CPU-h per genome: roughly 500 CPU-h, or about 120 task-hours, for all
1,000 genomes. With 60 concurrent tasks that is about 2 h of wall time once
the tasks are scheduled.

Runtime per genome is the wall clock of DIAMOND + HMMER + claim parsing +
`integrate` on 4 cores. The InterProScan time is not included because those
results were reused. Peak RAM is the largest max-RSS over these steps.
Input checksums are in `results/input_sha256_ibex.txt` and
`results/input_sha256_unimatrix.txt`.

## CAFA truth (`04_goa_exp_truth.py`)

RefSeq WP_ proteins were mapped to UniProtKB with NCBI's
`gene_refseq_uniprotkb_collab` (2,347,248 of 3,185,556 WP mapped). GOA rows
with experimental codes (EXP IDA IMP IPI IGI IEP HTP HDA HMP HGI HEP TAS IC)
were then collected for the mapped entries:

| Accepted when the annotated entry's taxon is … | genomes with ≥ 1 protein | proteins |
| --- | --- | --- |
| the assembly taxid | 15 | 125 |
| in the lineage at genus level or below | 17 | 130 |
| any taxon (an identical WP in another organism) | 66 | 741 |

That is 0.13 experimentally annotated proteins per genome on average, and
985 genomes have none. Genome-level F-max is not defensible on this truth.
The few covered proteins are well-studied Swiss-Prot entries, which sit in
GSPA's DIAMOND reference and in the tools' training data. `cafa_macro_fmax`
and `smin` are therefore `NA` for all eight rows. The old values for the five
deep-learning tools match `essential_terms/new_GO/average_completeness.tsv`
in the audit, a table without a producer script that was not computed on
these genomes. The InterProScan and GSPA values had no source at all.

## Findings relative to the old table

* **Exact matches.** Essential core-22, complex coherence, MetaCyc
  completeness and has_part coherence reproduce the old values exactly for
  all six tools. Both the matrix reproduction and the unified suite agree:
  the suite matches the audit's per-genome values for all 6,000
  genome × tool pairs.
* **Essential full.** The old column holds 32 terms, not 29. Its values came
  from matrices propagated with an older GO release; the GO 2025-03-16 values
  differ by up to 0.9 points. The old InterProScan matrix also contained one
  genome outside the 1,000 (GCF_000010245.2), and the old SPROF-GO matrix
  contained E. coli K-12 (GCF_000005845).
* **Taxon column.** The old "taxon-violating rows per genome" were the
  number of protein–GO rows (propagated) whose term carries any taxon
  constraint, violated or not. For DeepGOMeta that is 14,198.84; its real
  violations are 109.35 propagated rows per genome. All six old values
  reproduce exactly as constrained-row counts.
* **CAFA, Smin, runtime, RAM.** The tool values come from a summary table
  with no producer (CAFA) or have no source at all (runtime, RAM, the
  InterProScan CAFA value). All are now `NA`.
* **DeepFRI identifiers.** DeepFRI predictions are keyed by UniProt accession
  (structure models), not by RefSeq WP_ id. The structural suite works at
  genome level and is unaffected. The GSPA-evaluator step lists those
  accessions explicitly (`03_eval_chunk.sh`), because the first pass dropped
  all DeepFRI rows (completeness 0).

## Running

```bash
P=paper-reproducibility/pipelines/P12_refseq1000
bash $P/run.sh sync && bash $P/run.sh goa_subset     # then, when finished:
bash $P/run.sh stage && bash $P/run.sh prep
bash $P/run.sh matrix
bash $P/run.sh eval tools
bash $P/run.sh test                                  # 5 genomes, debug
bash $P/run.sh gspa 60                               # 100 tasks x 10 genomes
bash $P/run.sh eval gspa                             # after the GSPA array
bash $P/run.sh truth                                 # unimatrix01
bash $P/run.sh aggregate && bash $P/run.sh fetch
```

## Job IDs (2026-10-19)

See `OUTPUTS.md` § Job log.
