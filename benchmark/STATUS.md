# GSPA Benchmark Status

Current snapshot of where the 9-genome benchmark stands, plus the work
done since the initial head-to-head PGAP comparison.

## Pipeline state (where we are)

1. **9 genomes downloaded and prepared** — `ecoli K-12`, `hpylori`,
   `mgenitalium`, `mjannaschii`, `ecolo157` (E. coli O157:H7 Sakai),
   `bsubtilis` (B. subtilis 168), `mtb` (M. tuberculosis H37Rv),
   `synechocystis sp. PCC 6803`, `paeruginosa` (P. aeruginosa PAO1).
   UniProt proteome FASTA + RefSeq xref + NCBI `_genomic.gff` + NCBI
   `_protein.faa` + NCBI `_genomic.fna` for all 9.

2. **Leave-9-out DIAMOND + MMseqs2 reference DB** built on Swiss-Prot
   with the 9 genomes' accessions excluded (556k proteins).

3. **Predictors** (`DIAMOND blastp` + `HMMER/Pfam`) run on all 9
   genomes against the leave-9-out database. Outputs live at
   `/data/hohndor/gspa/proteomes/<tag>_preds9/`.

4. **Dual ground truth** extracted from `goa_uniprot_all.gaf.gz` in a
   single pass: `*_truth_exp.tsv` (experimental evidence codes only) +
   `*_truth_all.tsv` (all non-NOT evidence including IEA).

5. **Claims parser** (`02b_parse_predictors_to_claims.py`) consumes the
   predictor outputs + a GOA subset filtered to the DIAMOND targets,
   producing `<tag>_claims.jsonl` for every genome.

6. **Baseline and prior-enabled integration** run for all 9 genomes via
   `gspa integrate`, producing `<tag>_integrated.tsv` and
   `<tag>_integrated_priors.tsv`.

7. **Bootstrap F-max + IC-recall** (`benchmark_pgap_v2.py`) computed
   against both `exp` and `all-GOA` truth sets. 200-bootstrap 95% CIs.

8. **GAEF metrics** (Completeness / Coherence / Consistency) computed
   for GSPA and PGAP via `gspa evaluate` against synthetic GFF/GAF
   derived from each annotation TSV.

9. **Ablation study** (`run_ablation.sh` + `compile_ablation.py`) —
   4 configurations per genome (DIAMOND-only, Pfam-only, combined,
   combined+priors) with bootstrap F-max for each. Results in
   `ABLATION_REPORT.txt`.

## Final tables already produced

- `FINAL_REPORT.txt` — F-max (exp + all-GOA) + GAEF + dark-matter
  strip test on the 9 genomes (priors-off configuration).
- `ABLATION_REPORT.txt` — per-predictor and per-prior ablation, all
  9 genomes, both truth sets.

Key numbers (full-GOA truth, 95% CI):
| Genome | GSPA | PGAP |
|---|---|---|
| hpylori | 0.754 [0.730, 0.775] | 0.316 [0.298, 0.336] |
| mgenitalium | 0.913 [0.897, 0.930] | 0.469 [0.446, 0.492] |
| mjannaschii | 0.641 [0.625, 0.668] | 0.285 [0.267, 0.303] |
| ecoli | 0.670 [0.662, 0.679] | — |
| ecolo157 | 0.835 [0.827, 0.845] | — |
| bsubtilis | 0.673 [0.659, 0.689] | — |
| mtb | 0.716 [0.705, 0.726] | — |
| synechocystis | 0.614 [0.599, 0.633] | — |
| paeruginosa | 0.601 [0.588, 0.616] | — |

GSPA 1.9×–2.4× over PGAP on every genome where PGAP had GO
annotations in its RefSeq GFF.

Ablation: DIAMOND alone is the dominant signal; adding Pfam gives
+0.03 to +0.07 F-max on full-GOA truth; priors had zero effect on
this run (see next section).

## Priors: why they did nothing in the ablation

The ablation was run with a toy 17-line `test_pathways.tsv` (5
pathways), no metabolic gaps file, and taxon constraints that were
never wired into `gspa integrate`. Every prior was silently no-op:

| Prior | Required input | Available? |
|---|---|---|
| EssentialityPrior | GO reasoner, essential profile | needs non-`--lite` mode |
| CoherencePrior | Pathway DB + has_part pairs | 5 pathways → never fired |
| ConsistencyPrior | SatConsistencyChecker | never wired |
| GapFillingPrior | MetabolicGap JSONL | no `--gaps` passed |
| GenomicContextPrior | Operons + pathway DB | operons OK, DB too thin |

## Work in flight (fixes the "priors = 0" artifact)

1. **Taxon constraints auto-loaded** from whichever GO OWL the user
   passes. If the file is `go.owl` (core), the SAT checker gets
   whatever constraints that release surfaces; if it's `go-plus.owl`,
   the extended axiom set is used. Landed in `IntegrateCommand`.

2. **KEGG pathway DB** — `build_kegg_pathway_tsv.py` converts KEGG's
   `ec-pathway` link + `pathway` list + GO's `ec2go` into a real
   `kegg_pathways.tsv`: 169 pathways × 5262 reactions × 4060 rows
   carrying GO terms. Replaces the 17-line test file.

3. **gapseq find** running on all 9 genomes in parallel
   (`run_gapseq.sh` → `/data/hohndor/gspa/proteomes/gapseq/<tag>/`).
   Outputs `<tag>-all-Pathways.tbl` + `<tag>-all-Reactions.tbl`. A
   gapseq quirk: the conda share dir is read-only, so we copied the
   install to `/data/hohndor/envs/gapseq-rw/` to make the uniprot
   sequence cache writable. Wall-clock estimate: ~5–7 hours for
   all 9 with 9-way parallelism on the unimatrix01 16-core box.

4. **`parse_gapseq_gaps.py`** parses the gapseq tables into the
   `MetabolicGap` JSONL format that `gspa integrate --gaps` consumes.
   A gap = any reaction in a `Prediction=true` pathway whose status
   is not `good_blast`; `bad_blast` rows flagged as
   `gapseq_guessed=true` so the `GapFillingPrior` can discount them.

5. **Smoke-test confirmed the wiring works end-to-end** — using a
   hand-crafted 3-gap JSONL on hpylori, both `GapFillingPrior`
   (6 function-level boosts) and `GenomicContextPrior` (205 annotation-
   level boosts) fire and produce non-zero movement in the posterior
   (iter 0 ∆p = 0.00087). The smoke test is
   `benchmark/smoke_test_priors.sh`.

6. **Orchestrator** (`orchestrate_priors.sh`) is sitting in the
   background waiting for all 9 gapseq runs to land their
   `*-all-Pathways.tbl` files, then will run:
    - `parse_gapseq_gaps.py` per genome
    - `run_integrate_full_priors.sh` (all 5 priors, go-plus.owl, KEGG
      pathways, operons, gapseq gaps)
    - `run_fmax_full_priors.sh` (bootstrap F-max combined vs full-
      priors vs PGAP)
    - `print_fmax.py` → `FULL_PRIORS_FMAX_SUMMARY.txt`

## Known issues / open loops

- Only ~28 GO `only_in_taxon` axioms come out of `go-plus.owl` via
  our current extractor. The real release has hundreds. The
  extractor uses `SubClassOf ObjectSomeValuesFrom(...)`, which
  doesn't match GO's current encoding of taxon constraints. This is
  a follow-up.
- `EssentialityPrior` still no-ops under `--lite` because it needs
  the ELK reasoner to enumerate descendants of essential GO terms.
  The reasoner cache at `/data/hohndor/gspa/reference/reasoner-cache`
  exists but `--lite` still skips ELK init entirely. For the upcoming
  full-priors run this is acceptable — EssentialityPrior is only
  one of the five.
- Dark-matter strip test on the ablation run returned 0 recovery,
  because the 5-pathway test DB had no pathways whose missing
  functions matched the stripped truth GO terms. Should improve
  substantially once the KEGG pathway DB is in play.

## Files added since the last commit

- `benchmark/STATUS.md` — this document.

## Files added in the last three commits

- `benchmark/FINAL_REPORT.txt`
- `benchmark/ABLATION_REPORT.txt`
- `benchmark/benchmark_pgap_v2.py`
- `benchmark/build_kegg_pathway_tsv.py`
- `benchmark/compile_ablation.py`
- `benchmark/compile_final.py`
- `benchmark/extract_goa_dual.py`
- `benchmark/make_synth_gff_gaf.py`
- `benchmark/orchestrate_priors.sh`
- `benchmark/parse_gapseq_gaps.py`
- `benchmark/print_fmax.py`
- `benchmark/run_fmax_full_priors.sh`
- `benchmark/run_gapseq.sh`
- `benchmark/run_integrate_full_priors.sh`
- `benchmark/smoke_test_priors.sh`
- `benchmark/strip_test.py`
- `gspa-cli/src/main/groovy/gspa/cli/GspaMain.groovy` (`--reasoner-cache` flag on `evaluate`)
- `gspa-cli/src/main/groovy/gspa/cli/IntegrateCommand.groovy` (auto SAT checker wiring)
- `gspa-core/src/main/groovy/gspa/metrics/QualityPipeline.groovy` (Date.format fix + reasoner cache setter)
- `gspa-core/src/main/groovy/gspa/ontology/GoOntology.groovy` (getAxioms fix)
