# GSPA — Genome-Scale Protein Annotation

GSPA is a functional-annotation pipeline for prokaryotic and archaeal
genomes, MAGs, and microbial communities. It is intended as a drop-in
replacement for pipelines like Prokka / PGAP / bakta but extends them
in two directions:

1. **Multi-evidence Bayesian integration.** Rather than emitting the
   union of what each predictor produced, GSPA combines evidence from
   sequence similarity, protein domains, orthology, genomic context,
   metabolic-model gap analysis, and (optionally) structure/DL models
   via a log-odds Noisy-OR combiner with correlation-group collapse.
2. **Ontology-derived priors.** GO-level quality metrics
   (essentiality, coherence, consistency) and metabolic-model gaps
   feed back into annotation as Bayesian priors that boost weak
   evidence when it closes a pathway, fills an essential function, or
   matches a gapseq-identified missing reaction — and downweight
   claims that violate taxon constraints.

It ships with the full upstream predictor stack (DIAMOND, MMseqs2,
HMMER/Pfam, FoldSeek, eggNOG-mapper, InterProScan, AMRFinder, dbCAN,
antiSMASH, gapseq, OperonPredictor, SignalP, DeepTMHMM, …) wired
through a `Predictor` / `GenomePredictor` interface, and three core
quality metrics (Completeness, Coherence, Consistency) drawn from
Hoehndorf et al.'s genome-annotation quality framework.

## Project layout

Multi-module Gradle (Groovy 4 + Java 21):

- **gspa-core** — data model, GO ontology (OWL API + ELK), quality
  metrics (SAT4J for taxon-constraint consistency), I/O, config,
  Phase 7 integration engine (`gspa.integration`), Phase 8 dark-matter
  suggester (`gspa.integration.suggester`).
- **gspa-predictors** — `Predictor` interface, `AbstractToolPredictor`,
  gene callers, all tool wrappers, the `AnnotationPipeline`
  orchestrator, and the community / crossfeeding analyzer.
- **gspa-cli** — picocli CLI: `annotate`, `evaluate`, `compare`,
  `report`, `integrate`.
- **gspa-nf** — Nextflow workflow that runs the full predictor stack in
  containers (Docker + BioContainers).
- **benchmark/** — Standalone Python scripts + shell drivers for the
  9-genome head-to-head benchmark against PGAP (see below). Not part
  of the Gradle build.

## Evidence integration (Phase 7)

Each predictor emits `Annotation` objects that are lifted into
`EvidenceClaim` records (function type, function ID, protein ID,
calibrated raw score, evidence type, source). The integrator:

1. Groups claims by `(protein, function)`.
2. Within each group, collapses claims by evidence-correlation group
   (homology, structure, context, ml-sequence, localization, …) so
   DIAMOND + Pfam + eggNOG don't triple-count the same homology
   signal.
3. Combines via log-odds **Noisy-OR** with per-evidence-type
   reliability weights.
4. Runs up to 6 rounds of iterative refinement, each round:
   - Re-evaluates the priors against the current MAP annotation set.
   - Adds per-claim prior boosts: `L_post = L_lik + Σ λ_k · prior_k`.
   - Damping 0.5 (Jacobi) with divergence rollback.
5. Writes `IntegratedAnnotationSet` with full provenance: which
   predictors fired, which priors contributed how much, convergence
   iteration, final log-odds.

### Priors

| Prior | What it does | Needs |
|---|---|---|
| `EssentialityPrior` | For every uncovered essential function, boost candidate claims whose GO term is a descendant of it | essential profile + GO reasoner |
| `CoherencePrior` | For every process / pathway with partial coverage, boost the missing GO terms by `(1 − fraction_annotated)` | has_part pairs + pathway DB |
| `ConsistencyPrior` | Downweight (default soft, `−3` log-odds) GO terms that appear in any SAT4J UNSAT core against taxon constraints | SatConsistencyChecker + GO taxon axioms |
| `GapFillingPrior` | Boost claims whose EC / GO matches a gapseq-identified missing reaction (gapseq guesses get a 0.7× discount) | gapseq output |
| `GenomicContextPrior` | Per-operon pathway consensus; boost weak claims that match the top pathway's required set, with extra weight if the same GO closes a gapseq gap | operon file + pathway DB |

All priors are gated on input availability — each is silent if its
data is missing, so `gspa integrate` degrades gracefully on incomplete
inputs. Default prior weights come from the benchmark or from
hand-tuned defaults in `IntegrateCommand`.

## Dark-matter suggester (Phase 8)

Bayesian "dark matter" contextual-gap assigner. Given a metabolic gap
`(pathway P, reaction R, target function f_R)`, it:

1. Scores every operon with a Bayes factor `BF(O, P)` that the operon
   participates in pathway P, using the current posteriors of its
   members against `π_0(f)` background.
2. Within each operon passing `BF ≥ BF_min`, computes a per-protein
   log-odds `L_R(p) = L_likelihood + L_operon + L_lm` — plus a
   commitment penalty so proteins already strongly annotated for
   *other* pathway functions aren't over-selected.
3. Softmax over the operon members to get `q(p)`.
4. If `q(top) > 0.5`, emits a `SingletonSuggestion`. Otherwise takes
   the smallest top-k whose cumulative `q ≥ coverage_threshold` (=0.9
   by default) and emits a `DisjunctiveSuggestion`.

Every suggestion carries the motivating gap, the Bayes factor, the
per-protein log-odds decomposition, and the `q` distribution.

## Quality metrics (GAEF)

Ontology-based genome-annotation-quality metrics:

- **Completeness** — essential-function profiles per kingdom. What
  fraction of the expected GO terms for this organism are covered?
- **Coherence** — process / pathway / complex coherence.
  - Process: for every `C SubClassOf has_part some F` where `C` is
    present, is `F` also present? (ELK-extracted)
  - Pathway: fraction of each required pathway's GO terms covered.
  - Complex: fraction of each protein-complex's members covered.
- **Consistency** — taxon constraints as a propositional-SAT problem
  (SAT4J). UNSAT core identifies the conflicting annotations.
- **Information content** — mean IC of the annotated set.
- **Composite score** — weighted combination.

## Benchmark

The benchmark directory contains a self-contained 9-genome head-to-head
evaluation of GSPA vs PGAP, run on a leave-9-out Swiss-Prot reference
(so no genome in the evaluation contributes to its own DIAMOND hits).

Genomes:
- `ecoli` — E. coli K-12 (GCF_000005845.2)
- `hpylori` — H. pylori 26695 (GCF_000008525.1)
- `mgenitalium` — M. genitalium G37 (GCF_000027325.1)
- `mjannaschii` — M. jannaschii DSM 2661 (GCF_000091665.1, Archaea)
- `ecolo157` — E. coli O157:H7 Sakai (GCF_000008865.2)
- `bsubtilis` — B. subtilis 168 (GCF_000009045.1)
- `mtb` — M. tuberculosis H37Rv (GCF_000195955.2)
- `synechocystis` — Synechocystis sp. PCC 6803 (GCF_000009725.1)
- `paeruginosa` — P. aeruginosa PAO1 (GCF_000006765.1)

Two ground-truth sets from a single GOA scan:
- `*_truth_exp.tsv` — experimental evidence codes only (EXP, IDA, IMP,
  IPI, IGI, IEP, HTP, HDA, HMP, HGI, HEP, TAS, IC)
- `*_truth_all.tsv` — all non-NOT evidence codes (IEA included)

### Key results

F-max (bootstrap 200×, 95% CI) against full-GOA truth:

| Genome | GSPA | PGAP | Ratio |
|---|---|---|---|
| hpylori | **0.754** [0.730, 0.775] | 0.316 [0.298, 0.336] | 2.4× |
| mgenitalium | **0.913** [0.897, 0.930] | 0.469 [0.446, 0.492] | 1.9× |
| mjannaschii | **0.641** [0.625, 0.668] | 0.285 [0.267, 0.303] | 2.2× |
| ecoli | 0.670 [0.662, 0.679] | — | — |
| ecolo157 | 0.835 [0.827, 0.845] | — | — |
| bsubtilis | 0.673 [0.659, 0.689] | — | — |
| mtb | 0.716 [0.705, 0.726] | — | — |
| synechocystis | 0.614 [0.599, 0.633] | — | — |
| paeruginosa | 0.601 [0.588, 0.616] | — | — |

PGAP could only be compared where its RefSeq GFF carried GO
annotations (3 of the 9 genomes). On every other genome we report
GSPA-alone numbers.

GAEF comparison on the 3 PGAP-annotated genomes (run via
`gspa evaluate` with full ELK + curated pathways):

| Genome | Method | Completeness | Process coh. | Pathway coh. | Composite |
|---|---|---|---|---|---|
| hpylori | GSPA | 0.844 | 0.951 | 0.667 | 0.902 |
| hpylori | PGAP | 0.719 | 0.918 | 0.396 | 0.824 |
| mgenitalium | GSPA | 0.750 | 0.952 | 0.694 | 0.878 |
| mgenitalium | PGAP | 0.688 | 0.889 | 0.444 | 0.817 |
| mjannaschii | GSPA | 0.840 | 0.887 | 0.639 | 0.889 |
| mjannaschii | PGAP | 0.960 | 0.862 | 0.375 | 0.886 |

Full tables are in:
- [`benchmark/FINAL_REPORT.txt`](benchmark/FINAL_REPORT.txt) — F-max
  (both truth sets), GAEF metrics, dark-matter strip test.
- [`benchmark/ABLATION_REPORT.txt`](benchmark/ABLATION_REPORT.txt) —
  4-way ablation (DIAMOND-only / Pfam-only / combined / combined +
  priors) across all 9 genomes and both truth sets.
- [`benchmark/STATUS.md`](benchmark/STATUS.md) — current state of the
  benchmark pipeline and the in-flight full-priors run.

### Ablation headline

Full-GOA F-max per configuration (95% CI omitted for brevity):

| Genome | DIAMOND | Pfam | Combined | +priors* |
|---|---|---|---|---|
| ecoli | 0.656 | 0.137 | **0.670** | 0.670 |
| hpylori | 0.754 | 0.000 | **0.754** | 0.754 |
| mgenitalium | 0.913 | 0.000 | **0.913** | 0.913 |
| mjannaschii | 0.641 | 0.000 | **0.641** | 0.641 |
| ecolo157 | 0.807 | 0.256 | **0.835** | 0.835 |
| bsubtilis | 0.609 | 0.257 | **0.673** | 0.673 |
| mtb | 0.681 | 0.240 | **0.716** | 0.716 |
| synechocystis | 0.564 | 0.268 | **0.614** | 0.614 |
| paeruginosa | 0.597 | 0.037 | **0.601** | 0.601 |

*The `+priors` column shows zero effect because the ablation ran with
a toy 5-pathway test file and no `--gaps`, so four of the five priors
had no data to act on. A full-priors run with the real KEGG pathway DB,
go-plus.owl taxon constraints, and live gapseq output is in flight —
see [`benchmark/STATUS.md`](benchmark/STATUS.md).

## Build & test

```bash
./gradlew build            # compile + unit tests
./gradlew clean test       # fresh build
./gradlew :gspa-cli:shadowJar   # fat jar at gspa-cli/build/libs/gspa-0.1.0-SNAPSHOT.jar
```

Requires Java 21+ (the Gradle wrapper is included).

## Running

```bash
# 1. Full end-to-end annotation (runs predictors + quality metrics)
java -jar gspa-cli/build/libs/gspa-0.1.0-SNAPSHOT.jar annotate \
  -i genome.fna -o out/ --kingdom bacteria --go-owl go.owl

# 2. Quality evaluation only (from a GFF3 + GAF you already have)
java -jar ... evaluate \
  -i genome.gff -a annotations.gaf --go-owl go-plus.owl \
  --ec2go ec2go.txt --pathways kegg_pathways.tsv \
  --reasoner-cache ./reasoner-cache \
  -o quality_report.json

# 3. Phase 7 integration from pre-parsed claims (benchmark / BO driver)
java -jar ... integrate \
  --claims claims.jsonl \
  --out integrated.tsv \
  --go-owl go-plus.owl \
  --essential-profile bacteria \
  --pathways kegg_pathways.tsv --ec2go ec2go.txt \
  --operons operons.tsv --gaps gapseq_gaps.jsonl \
  --reasoner-cache ./reasoner-cache \
  --enable-priors essentiality,coherence,consistency,gap_filling,genomic_context \
  --dark-matter --suggestions-out suggestions.tsv
```

## Conventions

- Groovy source, Spock tests, Gradle Kotlin DSL build files
- Package root: `gspa.*`
- YAML config with hierarchical merging:
  `defaults → kingdom preset → user file → CLI overrides`
- All predictors implement `Predictor`; genome-scale ones implement
  `GenomePredictor`; external tools wrap `AbstractToolPredictor`
  with a command + parse-pattern.

## License

See `LICENSE` once added. Uses OWL API (LGPL), ELK (Apache), SAT4J
(LGPL/EPL), picocli (Apache), Jackson (Apache), Spock (Apache).

## Citation

Paper in preparation. If you use GSPA, please cite this repository
and the quality-metrics paper the GAEF framework is drawn from.
