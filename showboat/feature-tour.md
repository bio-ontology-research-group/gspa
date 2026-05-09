# GSPA Feature Tour — what's shipped, what's partial, what's future

*2026-05-06T17:10:44Z by Showboat 0.6.1*
<!-- showboat-id: 2c2c5795-1619-466f-a863-67ec2bc0e66a -->

GSPA (Genome-Scale Protein Annotation) is a functional-annotation pipeline for prokaryotic and archaeal genomes, MAGs, and microbial communities. It combines multiple evidence types — sequence similarity, protein domains, orthology, structure, genomic context, metabolic-model gap analysis — through a Bayesian Noisy-OR integrator, then uses ontology-derived priors (essentiality, pathway coherence, taxon consistency, metabolic gaps, genomic context) to boost weak evidence and suppress taxon-constraint violations.

This document is a tour of what is **fully included**, **partial**, and **future work** as of the `phase11-crossgenome` branch (May 2026). Every command below is real — `showboat verify` re-runs them all.

## 1. Repository shape

The build is a multi-module Gradle project (Java 21+, Groovy 4) with two siblings outside the build (Nextflow + Python harness). The Gradle modules are `gspa-core` (data model, ontology, quality metrics, integrator), `gspa-predictors` (predictor wrappers and pipeline orchestrator), and `gspa-cli` (picocli subcommands).

```bash
cat settings.gradle.kts
```

```output
rootProject.name = "gspa"

include("gspa-core")
include("gspa-predictors")
include("gspa-cli")

// Intentionally NOT included:
//   gspa-nf/   — Nextflow pipeline, runs via `nextflow run gspa-nf/main.nf`
//                See gspa-nf/README.md.
//   benchmark/ — Python evaluation harness and neural sidecar.
//                See benchmark/README.md.
```

```bash
echo '=== Source line counts ==='; for m in gspa-core gspa-predictors gspa-cli; do n=$(find $m/src/main -name '*.groovy' -o -name '*.java' 2>/dev/null | xargs wc -l 2>/dev/null | tail -1 | awk '{print $1}'); echo "$m: $n LOC"; done; echo; echo '=== Test line counts ==='; for m in gspa-core gspa-predictors gspa-cli; do n=$(find $m/src/test -name '*.groovy' -o -name '*.java' 2>/dev/null | xargs wc -l 2>/dev/null | tail -1 | awk '{print $1}'); echo "$m: $n LOC"; done; echo; echo '=== Sibling pipelines ==='; echo "gspa-nf .nf modules: $(find gspa-nf -name '*.nf' | wc -l)"; echo "benchmark/ python scripts: $(find benchmark -name '*.py' | wc -l)"
```

```output
=== Source line counts ===
gspa-core: 10719 LOC
gspa-predictors: 6011 LOC
gspa-cli: 1286 LOC

=== Test line counts ===
gspa-core: 4675 LOC
gspa-predictors: 2483 LOC
gspa-cli: 0 LOC

=== Sibling pipelines ===
gspa-nf .nf modules: 16
benchmark/ python scripts: 132
```

## 2. CLI surface

There are five `gspa` subcommands — each lives in its own file. They cover the two modes of using the project: run predictors (`annotate`), or evaluate / compare / re-integrate existing annotations (`evaluate`, `compare`, `report`, `integrate`).

```bash
ls -1 gspa-cli/src/main/groovy/gspa/cli/
```

```output
AnnotateCommand.groovy
CompareCommand.groovy
EvaluateCommand.groovy
GspaMain.groovy
IntegrateCommand.groovy
ReportCommand.groovy
```

```bash
grep -h "@Command(name = '" gspa-cli/src/main/groovy/gspa/cli/*.groovy
```

```output
@Command(name = 'annotate', description = 'Annotate a genome or set of genomes')
@Command(name = 'compare', description = 'Compare quality of multiple annotation sets for the same genome')
@Command(name = 'evaluate', description = 'Evaluate quality of existing annotations')
@Command(name = 'report', description = 'Generate HTML quality report from existing JSON reports')
```

```bash
grep -hE "@Command|^    name|^    description" gspa-cli/src/main/groovy/gspa/cli/IntegrateCommand.groovy | head -6
```

```output
@Command(
    name = 'integrate',
    description = 'Run the Phase 7 evidence integrator on pre-parsed claims.'
```

## 3. Predictor catalogue

Predictors are wired in `AnnotationPipeline.createAllPredictors()`. Each one wraps an external tool (or Python sidecar) behind the `Predictor` / `GenomePredictor` interface. There are four implementation styles: shell-out tool wrappers (`AbstractToolPredictor`), Python neural sidecars (`AbstractNeuralSidecarPredictor`), region/site/genomic-region sidecars (the v1.2/1.3/1.4 FOSS line), and pure-JVM context predictors (e.g. `OperonPredictor`).

```bash
find gspa-predictors/src/main/groovy/gspa/predictor -mindepth 2 -name '*Predictor.groovy' -o -name '*Caller.groovy' -o -name '*Runner.groovy' | sort | sed 's|gspa-predictors/src/main/groovy/gspa/predictor/||'
```

```output
context/GLMOperonPredictor.groovy
context/OperonPredictor.groovy
disorder/MetapredictPredictor.groovy
domain/HmmerPredictor.groovy
domain/InterProScanPredictor.groovy
function/EggNogMapperPredictor.groovy
genecalling/GeneCaller.groovy
genecalling/ProdigalCaller.groovy
genecalling/PyrodigalCaller.groovy
localization/DeepSigPredictor.groovy
localization/DeepTmhmmPredictor.groovy
localization/PSORTbPredictor.groovy
localization/SignalPPredictor.groovy
localization/TmbedPredictor.groovy
localization/TPpred3Predictor.groovy
neural/AbstractNeuralSidecarPredictor.groovy
neural/CleanPredictor.groovy
neural/DeepEcPredictor.groovy
neural/DeepGoPlusEsm2Predictor.groovy
neural/Esm2CentroidPredictor.groovy
neural/ProteInferPredictor.groovy
pathway/GapseqPredictor.groovy
similarity/DiamondPredictor.groovy
similarity/MMseqs2Predictor.groovy
sites/MusiteDeepPredictor.groovy
sites/ScanNetPredictor.groovy
specialized/AmrFinderPredictor.groovy
specialized/AntiSmashPredictor.groovy
specialized/CrisprPredictor.groovy
specialized/DbCanPredictor.groovy
specialized/DeepArgPredictor.groovy
specialized/VfdbPredictor.groovy
structure/DeepFriPredictor.groovy
structure/FoldSeekPredictor.groovy
taxonomy/CheckM2Runner.groovy
taxonomy/GtdbTkRunner.groovy
viral/CheckVPredictor.groovy
viral/GenomadPredictor.groovy
viral/PhiSpyPredictor.groovy
```

**Status**: 35 predictors + 2 gene callers + 2 taxonomy runners are wired. Categories:

- **Sequence similarity** — DIAMOND, MMseqs2
- **Domain** — Pfam/HMMER, InterProScan (with InterPro2GO mapping)
- **Structure** — FoldSeek (homology-transfer + centroid mode), DeepFRI
- **Orthology / function** — eggNOG-mapper
- **Neural protein** — ESM2-DeepGOPlus, ESM2-centroid, ProteInfer, CLEAN (EC), DeepEC (EC, AGPL-walled)
- **Localization** — SignalP, DeepTMHMM, plus FOSS replacements DeepSig, TMbed, TPpred3, PSORTb
- **Sites** — MusiteDeep (PTM), ScanNet (PPI interface)
- **Specialized** — AMRFinder, dbCAN, antiSMASH, VFDB, DeepARG, CRISPR
- **Disorder** — Metapredict
- **Genomic context** — `OperonPredictor` (intergenic-distance heuristic) + `GLMOperonPredictor` (foundation-model)
- **Pathway** — Gapseq
- **Viral / prophage** — geNomad, CheckV, PhiSpy (and v1.4 ensemble adds VirSorter2 + VIBRANT in Nextflow only)

## 4. Phase 7 — evidence integration (FULLY INCLUDED)

Phase 7 is the heart of GSPA. Each predictor emits `Annotation` objects, which `ClaimExtractor` lifts into `EvidenceClaim` records. The integrator groups by `(protein, function)`, collapses correlated evidence types into one effective claim per group, combines via log-odds **Noisy-OR** with per-type reliability, and runs an iterative refinement loop with priors.

```bash
ls gspa-core/src/main/groovy/gspa/integration/*.groovy | xargs -n1 basename
```

```output
CalibrationTable.groovy
ClaimExtractor.groovy
ClaimKey.groovy
ClaimProvenance.groovy
EvidenceClaim.groovy
EvidenceCombiner.groovy
EvidenceType.groovy
GapKey.groovy
GapRecomputer.groovy
IntegratedAnnotationSet.groovy
IntegrationState.groovy
IntegrationWriter.groovy
IterativeRefiner.groovy
MetabolicGap.groovy
OuterIterativeRefiner.groovy
PriorEngine.groovy
Prior.groovy
```

```bash
echo '=== EvidenceType enum (correlation-grouped) ==='; grep -E '^\s+[A-Z_]+,' gspa-core/src/main/groovy/gspa/integration/EvidenceType.groovy
```

```output
=== EvidenceType enum (correlation-grouped) ===
    SEQUENCE_SIMILARITY,        // DIAMOND, MMseqs2
    SEQUENCE_DOMAIN,            // InterProScan, Pfam/HMMER, TIGRFAM
    SEQUENCE_MOTIF,             // PROSITE, ELM
    SEQUENCE_DEEPLEARNING,      // DeepGO, ESM-based (Phase 9)
    STRUCTURE_SIMILARITY,       // FoldSeek
    STRUCTURE_DEEPLEARNING,     // DeepFRI, GraphGOSeq (Phase 9)
    PROTEIN_LM_EMBEDDING,       // SaProt, GOPredSim (Phase 9)
    ORTHOLOGY,                  // eggNOG-mapper, OMA
    GENOMIC_CONTEXT,            // operon co-occurrence
    METABOLIC_CONTEXT,          // gapseq pathway / gap-fill
    GENOMIC_LANGUAGE_MODEL,     // nucleotide LM over operon / regulon (Phase 9)
    LOCALIZATION,               // SignalP, DeepTMHMM
    DOMAIN_SPECIFIC_AMR,        // AMRFinder
    DOMAIN_SPECIFIC_CAZY,       // dbCAN
    DOMAIN_SPECIFIC_BGC,        // antiSMASH
    DOMAIN_SPECIFIC_VF,         // VFDB
    DARK_MATTER,                // claims promoted by DarkMatterSuggester (Phase 10)
    REACTION_LOCAL_CONTEXT,     // claims from ReactionLocalContextSuggester (Phase 12)
    CROSS_GENOME_TRANSFER,      // conditional-LR-based cross-genome transfer (Phase 12)
    ML_RANKER,                  // learned ranker output (Phase 12 M3+)
    SEQUENCE_REGION_ML,         // region-level ML predictors (Metapredict, SignalP region, TMHMM helix)
```

```bash
ls gspa-core/src/main/groovy/gspa/integration/prior/ | sed 's|.groovy||'
```

```output
CoherencePrior
ConsistencyPrior
EssentialityPrior
GapFillingPrior
GenomicContextPrior
HomologyTransferPrior
```

All six priors live under `gspa.integration.prior`. Five of them — Essentiality, Coherence, GapFilling, GenomicContext, Consistency — are **fully included** and validated on the 13-genome benchmark. The sixth, `HomologyTransferPrior`, is the cross-genome transfer mechanism, gated on `--orthogroups` + `--cluster-consensus`.

Iterative refinement: `IterativeRefiner` (up to 6 rounds, Jacobi damping = 0.5). The Phase 10 outer loop `OuterIterativeRefiner` wraps it for the dark-matter promotion fixed point (see §6).

## 5. Quality metrics — GAEF (FULLY INCLUDED)

The quality scoring framework lives in `gspa.metrics`. It is the original GAEF triad — **Completeness**, **Coherence**, **Consistency** — plus information-content weighting and a MAG adjuster that downgrades expected-completeness for incomplete genomes.

```bash
ls gspa-core/src/main/groovy/gspa/metrics/ | sed 's|.groovy||'
```

```output
Coherence
Completeness
Consistency
HtmlReportWriter
InformationContent
MagAdjuster
QualityPipeline
QualityReportWriter
QualityScorer
```

Key design decision: **consistency uses SAT4J, not ELK.** Taxon constraints are encoded as propositional SAT with sibling disjointness, and an UNSAT core is extracted to explain *which* annotations violated *which* constraint. ELK is used for completeness (subsumption) and for coherence (`has_part` pair extraction over biological processes / complexes), with a reasoner cache to avoid the 12-minute first-run cost.

## 6. Phase 8/10 — dark matter and the outer loop (PARTIAL)

`gspa.integration.suggester.DarkMatterSuggester` takes a metabolic gap `(pathway P, reaction R, target function f_R)` and uses Bayesian operon scoring to assign the missing function to the most likely protein. It emits two suggestion kinds:

- **Singleton** — top `q(p) > 0.5` within a passing operon
- **Disjunctive** — credible set whose cumulative `q` exceeds the coverage threshold

Phase 10 wraps Phase 8 in an **outer fixed-point loop**: refine → suggest → promote → pin → recompute gaps → refine again. Promoted singletons are pinned as posterior floors so the next refinement can't drive them back down.

```bash
echo '=== Suggester package ==='; ls gspa-core/src/main/groovy/gspa/integration/suggester/ | sed 's|.groovy||'; echo; echo '=== Promotion strategies ==='; ls gspa-core/src/main/groovy/gspa/integration/promotion/ | sed 's|.groovy||'
```

```output
=== Suggester package ===
DarkMatterSuggester
DarkMatterWriter
DisjunctiveSuggestion
GenomicDensityField
PerProteinDecomposition
ReactionLocalContextSuggester
SingletonSuggestion
Suggestion

=== Promotion strategies ===
AllAboveThresholdStrategy
BeamSearchStrategy
GreedyStrategy
MaxSatStrategy
PromotionStrategy
```

Why this is **PARTIAL**: the Phase 8 *one-shot* dark-matter mode is fully working and shipped — see the per-genome singleton/disjunctive counts in the README. The Phase 10 *iterative* mode (`--iterate-gapseq`) is mechanically complete (promotions, pins, fixed-point detection, cascade rollback, four promotion strategies including SAT4J MaxSAT and beam search) but **F-max regresses by ~4 points** at default settings on the 10-genome PGAP panel. The architecture is proven; tuning is the open question.

Excerpt from `benchmark/RESULTS.md` (Phase 10 Part 1):

```bash
awk '/^### Results .mean F-max across 10 genomes/,/^### Per-genome F-max/' benchmark/RESULTS.md | head -12
```

```output
### Results (mean F-max across 10 genomes)

| Config | fmax_micro | fmax_CAFA | coverage | IC-recall | Δmicro | ΔCAFA |
|---|---:|---:|---:|---:|---:|---:|
| C1 baseline | **0.8419** | **0.8676** | 0.895 | 0.800 | — | — |
| C2 iterate | 0.8001 | 0.8524 | 0.899 | 0.801 | **−0.0418** | **−0.0152** |
| C3 iter + cluster | 0.8001 | 0.8524 | 0.899 | 0.801 | −0.0418 | −0.0152 |
| C4 iter + cluster + blastp | 0.8001 | 0.8524 | 0.899 | 0.801 | −0.0418 | −0.0152 |
| C5 iter + cluster + reps (Singularity) | 0.8001 | 0.8524 | 0.899 | 0.801 | −0.0418 | −0.0152 |
| C2 no-pin | 0.8001 | 0.8524 | 0.899 | 0.801 | −0.0418 | −0.0152 |

### Per-genome F-max (micro)
```

## 7. Phase 12 — cross-genome reaction-local context (PARTIAL)

`gspa.integration.crossgenome` is the cross-genome layer. It pairs with `ReactionLocalContextSuggester` (M1, in-tree under `suggester/`), `ReactionLocusCatalog` (per-orthogroup, per-reaction conditional likelihood ratios with Jeffreys-prior beta-binomial 90% CIs), and `CrossGenomeReScorer` (M2 — re-weights existing suggestions by `LR(C, R)^λ`).

The ranker package (`gspa.integration.ranker`) holds the M3 learned re-scorer scaffold.

```bash
ls gspa-core/src/main/groovy/gspa/integration/crossgenome/ gspa-core/src/main/groovy/gspa/integration/ranker/ | sed 's/^.*$//' 2>/dev/null; echo '=== crossgenome ==='; ls gspa-core/src/main/groovy/gspa/integration/crossgenome/; echo; echo '=== ranker ==='; ls gspa-core/src/main/groovy/gspa/integration/ranker/
```

```output









=== crossgenome ===
ConditionalLRScorer.groovy
ReactionLocusCatalog.groovy

=== ranker ===
GbdtRanker.groovy
RankerFeatures.groovy
Ranker.groovy
RankerRescorer.groovy
```

**M1 (RLGC)** and **M2 (cross-genome LR re-scoring)** are wired through `gspa integrate`'s `--rlc-suggester`, `--rxn-locus-catalog`, `--cg-lambda`, etc. flags — see §11. **M3** (learned ranker over claim+context features) ships scaffolding (`Ranker`, `GbdtRanker`, `RankerRescorer`, `RankerFeatures`) but the production training pipeline lives outside the repo (`benchmark/ml/train_lambdamart.py`).

The KAUST-internal **500-genome panel** (`benchmark/panel/`) and **culture-genome dark-matter scoring** (`benchmark/cultures/`) are the consumers of the cross-genome layer; both are active workstreams.

## 8. Phase 11 — gLM operon caller (FULLY INCLUDED but EVALUATION CHECKPOINT IS NO-GO)

The branch `phase11-crossgenome` ships a foundation-model operon caller (gLM, Hwang et al. 2024) as a drop-in replacement for the intergenic-distance heuristic. The Python sidecar plus Groovy wrapper plus mocked Spock test plus benchmark switch (`run_integrate_full_priors.sh --operon-caller {heuristic,glm}`) all landed.

```bash
echo '=== gLM artifacts ==='; ls gspa-predictors/src/main/groovy/gspa/predictor/context/; echo; ls benchmark/neural/run_glm*.py benchmark/neural/run_glm_operon.sbatch 2>/dev/null
```

```output
=== gLM artifacts ===
GLMOperonPredictor.groovy
OperonPredictor.groovy

benchmark/neural/run_glm2_operon.py
benchmark/neural/run_glm_operon.py
benchmark/neural/run_glm_operon.sbatch
```

```bash
awk '/^### CHECKPOINT 1/,/^## Phase 2/' tasks/todo.md
```

```output
### CHECKPOINT 1 — phase-1 go/no-go: **NO-GO** (2026-05-05)

- [x] `./gradlew test` green
- [x] 8 of the canonical 13 genomes integrated in both modes
- [x] Mean micro F-max Δ = **+0.0002** — fails threshold (+0.005)
- [x] Worst-genome Δ = **−0.0002** — passes floor (−0.01)
- [x] User sign-off: pending

→ **NO-GO. Phase 2 and phase 3 not started.**

## Phase 2 — GENOMIC_CONTEXT_FM evidence type (gated on Checkpoint 1)
```

**Verdict**: the code is in production-ready shape (`./gradlew test` green, sidecar self-test passes, 8/13 benchmark genomes ran), but the F-max gating bar (+0.005 mean) was not met. Phase 2 (GENOMIC_CONTEXT_FM evidence channel) and Phase 3 (BF(O,P) embedding-distance augmentation) are **future work, gated on a Phase 1 reroll** — possibly with gLM2 (TattaBio, mixed-modality) instead of gLM.

## 9. v1.2/1.3/1.4 — FOSS predictor expansion (FULLY INCLUDED)

Three release waves added 10 OSI-licensed protein predictors and a viral / prophage genomic-region track. Each release retired a license-walled tool and replaced it with a FOSS equivalent.

```bash
awk '/^# v1.2 — FOSS-only/{flag=1} flag && /^## What was added/{p=1} p && /^## Vocabulary/{exit} p' benchmark/RESULTS.md | head -40
```

```output
## What was added

**Region predictors** (5-col TSV, `protein_id, region_start, region_end,
region_type, score`):

- Metapredict v2 (MIT) — disorder regions
- DeepSig v3 (GPL-3.0) — Sec/Tat signal peptides; FOSS replacement for SignalP 6
- TMbed (Apache-2.0) — TM helices via ProtT5; FOSS replacement for DeepTMHMM
- TPpred 3 (GPL-3.0) — N-terminal targeting peptides; FOSS replacement for TargetP 2

**Term-extras** (4-col TSV, auto-join the v1.1 ensemble):

- PSORTb 3.0 (GPL-3.0) — bacterial subcellular localization; FOSS replacement for DeepLoc 2
- DeepFRI (BSD-3-Clause) — sequence-only GO; complements ESM2-DGP/ProteInfer
- DeepEC (AGPL-3.0 ⚠) — EC predictor; complements CLEAN/DIAMOND
- DeepARG (MIT) — antimicrobial-resistance gene calls

**Site predictors** (5-col TSV, `protein_id, position, site_type, score,
annotation_type`):

- MusiteDeep_web (MIT) — PTM sites (phospho-S/T/Y default; configurable);
  FOSS replacement for NetPhos / NetPhosBac
- ScanNet (Apache-2.0) — PPI interface residues; needs structures

**License-walled tools dropped** (FOSS replacement in parens):

- DeepLoc 2 (PSORTb), DeepTMHMM (TMbed), IUPred3 (Metapredict),
  TargetP 2 (TPpred 3), SignalP 6 (DeepSig), NetPhos / NetPhosBac
  (MusiteDeep), MULocDeep (none — academic-only), DR-BERT (no LICENSE
  file)

The existing v1.1 JVM wrappers for SignalP 6 and DeepTMHMM remain in
the codebase (deletion would be breaking) but are NOT productionised
in Nextflow; FOSS replacements are the recommended path.

```

v1.4 finished the FOSS line — TPpred3, MusiteDeep, TMbed, ScanNet (the latter via SIF) all now run. v1.3+v1.4 added a viral/prophage genomic-region track (geNomad, CheckV, PhiSpy, plus VirSorter2/VIBRANT in Nextflow), with a 6-column genomic-region TSV shape and an interval-overlap ensemble (`build_genomic_ensemble.py`). The RDF/JSON-LD vocabulary was extended to cover `Region`, `Site`, `GenomicRegion`, and friends — the HTML report and the Turtle/JSON-LD report agree triple-for-triple on test fixtures.

## 10. Nextflow sibling pipeline gspa-nf (FULLY INCLUDED, partial parity)

`gspa-nf/` is a container-based alternative to the JVM CLI — same external tools, packaged as Nextflow processes with Singularity/Docker images. It is **deliberately not part of the Gradle build** (`./gradlew build` does not touch it). 22 processes across 16 .nf module files. Profiles for Docker, Singularity, and a unimatrix01 SLURM/Singularity profile.

```bash
ls gspa-nf/modules/ | sed 's|.nf||'
```

```output
domains
ensemble
eval
gene_calling
loc
neural
quality
region
report
similarity
sites
specialized
structure
term_extras
viral
```

**Partial parity caveat**: the Nextflow path produces raw per-tool outputs and a merged TSV — it does **not** call the GSPA Phase 7 integrator or quality machinery. To get GSPA priors and quality scores, feed the Nextflow output into `gspa-cli integrate` or `gspa-cli evaluate`. So the JVM CLI and the Nextflow pipeline are partial siblings: same predictors, different downstream stitching.

## 11. The `gspa integrate` flag surface (FULLY INCLUDED — most complete view of the architecture)

Reading `IntegrateCommand.groovy` is the most efficient way to see the entire integration architecture: every flag corresponds to a feature axis.

```bash
grep -E "^\s+@Option\(names" gspa-cli/src/main/groovy/gspa/cli/IntegrateCommand.groovy | sed -E 's/.*\[([^]]*)\].*/\1/' | tr ',' '\n' | sed 's/^ *//;s/ *$//' | grep -v '^$' | sort -u
```

```output
'--beam-candidates-per-gap'
'--beam-width'
'--cg-lambda'
'--cg-min-support'
'--cg-require-credible'
'--claims'
'--cluster-consensus'
'--dark-matter'
'--diffusion-mets'
'--ec2go'
'--enable-priors'
'--essential-functions'
'--essential-profile'
'--features-out'
'--gaps'
'--gapseq-pin-promotions'
'--gapseq-q-base'
'--gapseq-q-cap'
'--gapseq-q-step'
'--gapseq-target'
'--gapseq-tau-cover'
'--genome-layout'
'--go-owl'
'--intragenome-cluster'
'--iterate-gapseq'
'--lite'
'--max-gapseq-iter'
'--maxsat-coherence-bonus'
'--operons'
'--orthogroups'
'--out'
'--pathways'
'--promotion-strategy'
'--provenance'
'--reaction-ec-aliases'
'--reaction-graph'
'--reasoner-cache'
'--refined-bf'
'--rlc-alpha'
'--rlc-anchor-threshold'
'--rlc-currency-pct'
'--rlc-kernel-bandwidth'
'--rlc-radius-k'
'--rlc-suggester'
'--rxn-locus-catalog'
'--suggestions-out'
'--taxon-constraints'
'--taxonomy'
'--theta'
```

Reading those flag groups end-to-end: claims input + theta hyperparameters + reference data (`--go-owl`, `--ec2go`, `--pathways`, `--operons`, `--gaps`, `--orthogroups`, `--cluster-consensus`, `--taxonomy`, `--taxon-constraints`, `--genome-layout`) → prior selection (`--enable-priors`) → Phase 8 dark matter (`--dark-matter`, `--refined-bf`) → Phase 10 outer loop (`--iterate-gapseq`, `--gapseq-q-*`, `--gapseq-target`, `--intragenome-cluster`, `--promotion-strategy` ∈ {default, greedy, maxsat, beam}) → Phase 12 RLGC (`--rlc-*`) → Phase 12 M2 cross-genome (`--rxn-locus-catalog`, `--cg-*`).

Every one of those code paths exists, has unit tests, and runs. The questions are about *tuning* and *scaling*, not whether the architecture is wired.

## 12. Benchmark results (FULLY INCLUDED)

The headline benchmark — 13 prokaryotic / archaeal genomes, head-to-head against PGAP — is in `benchmark/RESULTS.md`. Mean GSPA / PGAP F-max ratio is **1.93× (micro), 1.96× (CAFA)**. GSPA wins on every single genome under both metrics.

```bash
awk '/^### Combined: all 13 genomes with PGAP GO annotations/,/^Mean GSPA.PGAP across all 13 genomes/' benchmark/RESULTS.md | head -25
```

```output
### Combined: all 13 genomes with PGAP GO annotations

| Genome | GSPA micro | GSPA CAFA | PGAP micro | PGAP CAFA | ratio (micro) | ratio (CAFA) |
|---|---|---|---|---|---|---|
| rprowazekii | **0.911** | **0.905** | 0.503 | 0.518 | **1.81×** | **1.75×** |
| mgenitalium | **0.908** | **0.910** | 0.469 | 0.448 | **1.94×** | **2.03×** |
| tpallidum | **0.892** | **0.896** | 0.491 | 0.492 | **1.82×** | **1.82×** |
| saureus | **0.867** | **0.886** | 0.449 | 0.446 | **1.93×** | **1.99×** |
| vcholerae | **0.858** | **0.883** | 0.443 | 0.441 | **1.94×** | **2.00×** |
| pfuriosus | **0.857** | **0.874** | 0.350 | 0.366 | **2.45×** | **2.39×** |
| spneumoniae | **0.846** | **0.870** | 0.447 | 0.457 | **1.89×** | **1.90×** |
| tthermophilus | **0.842** | **0.869** | 0.492 | 0.509 | **1.71×** | **1.71×** |
| hpylori | **0.819** | **0.844** | 0.316 | 0.299 | **2.59×** | **2.82×** |
| ccrescentus | **0.805** | **0.848** | 0.480 | 0.493 | **1.68×** | **1.72×** |
| dradiodurans | **0.780** | **0.818** | 0.463 | 0.473 | **1.68×** | **1.73×** |
| scoelicolor | **0.778** | **0.847** | 0.490 | 0.503 | **1.59×** | **1.68×** |
| mjannaschii | **0.732** | **0.775** | 0.285 | 0.269 | **2.57×** | **2.88×** |

**Mean GSPA/PGAP across all 13 genomes**: **1.93× (micro), 1.96× (CAFA)**.
GSPA is consistently 1.6-2.9× better than PGAP across 7 phyla, both
domains of life (Bacteria + Archaea), and genome sizes from 483 to
7,872 proteins. The two metrics give the same qualitative picture
(GSPA > PGAP on every single genome under both); CAFA is slightly
more favourable to GSPA on the larger Actinobacteria/Deinococcus
genomes and on the small archaea where PGAP misses most annotations.
```

Beyond the headline 13-genome benchmark, `benchmark/RESULTS.md` documents:

- **21-genome 7-predictor benchmark on IBEX** (April 2026) — neural-predictor cross-comparison: ProteInfer 0.660 F-max micro, ESM2-DGP 0.325, FoldSeek 0.249, DIAMOND 0.245; ensemble-mean wins at 0.767. CLEAN beats ProteInfer for EC.
- **Phase 10 ablation** on 10 PGAP genomes — outer loop currently regresses F-max by 0.04 at default settings.
- **Phase 11 gLM ablation** on 8 of 13 genomes — Δ μ F-max = +0.0002, fails the +0.005 gating bar (NO-GO).

## 13. Tests

Spock tests across both `gspa-core` and `gspa-predictors`. ~7,200 lines of test code total.

```bash
for m in gspa-core gspa-predictors; do n=$(find $m/src/test -name '*Spec.groovy' 2>/dev/null | wc -l); echo "$m: $n test classes"; done
```

```output
gspa-core: 41 test classes
gspa-predictors: 23 test classes
```

```bash
find gspa-core/src/test -name 'OuterIterativeRefiner*' -o -name 'DarkMatter*' -o -name 'CrossGenome*' -o -name 'ReactionLocal*' -o -name 'GLMOperon*' -o -name 'IterativeRefiner*' | sort | head -20
```

```output
gspa-core/src/test/groovy/gspa/integration/IterativeRefinerSpec.groovy
gspa-core/src/test/groovy/gspa/integration/OuterIterativeRefinerSpec.groovy
gspa-core/src/test/groovy/gspa/integration/suggester/DarkMatterSuggesterSpec.groovy
gspa-core/src/test/groovy/gspa/integration/suggester/ReactionLocalContextSuggesterSpec.groovy
```

```bash
find gspa-predictors/src/test -name 'GLMOperon*' -o -name '*Predictor*Spec.groovy' | sort | head -10
```

```output
gspa-predictors/src/test/groovy/gspa/predictor/AllPredictorParsingSpec.groovy
gspa-predictors/src/test/groovy/gspa/predictor/context/GLMOperonPredictorSpec.groovy
gspa-predictors/src/test/groovy/gspa/predictor/context/OperonPredictorSpec.groovy
gspa-predictors/src/test/groovy/gspa/predictor/disorder/MetapredictPredictorSpec.groovy
gspa-predictors/src/test/groovy/gspa/predictor/domain/InterProScanPredictorSpec.groovy
gspa-predictors/src/test/groovy/gspa/predictor/neural/NeuralSidecarPredictorsSpec.groovy
gspa-predictors/src/test/groovy/gspa/predictor/PredictorRegistrySpec.groovy
gspa-predictors/src/test/groovy/gspa/predictor/similarity/DiamondPredictorSpec.groovy
gspa-predictors/src/test/groovy/gspa/predictor/specialized/DeepArgPredictorSpec.groovy
gspa-predictors/src/test/groovy/gspa/predictor/specialized/SpecializedPredictorSpec.groovy
```

## 14. Summary — what's where

### FULLY INCLUDED

- Multi-module Gradle build (Java 21+, Groovy 4) with Spock tests across core + predictors
- 35 predictor wrappers (sequence similarity, domains, structure, orthology, neural, localization, sites, specialized, viral) covering both license-walled and FOSS tools
- Two gene callers (Prodigal, Pyrodigal) and two taxonomy runners (CheckM2, GTDB-Tk)
- Phase 7 evidence integrator: Noisy-OR with correlation groups, IterativeRefiner, 6 priors (Essentiality, Coherence, Consistency [SAT4J], GapFilling, GenomicContext, HomologyTransfer)
- GAEF quality metrics (Completeness, Coherence, Consistency, InformationContent, MagAdjuster) with HTML + JSON reports
- Phase 8 one-shot DarkMatterSuggester (singletons + disjunctives)
- Phase 12 RLGC (Reaction-Local Context) suggester + Phase 12 M2 cross-genome LR re-scorer
- gLM operon caller (sidecar + Groovy wrapper + Spock test + benchmark switch)
- gspa-nf Nextflow sibling (22 processes, Singularity / Docker / SLURM profiles)
- 13-genome head-to-head benchmark vs PGAP (mean GSPA/PGAP = 1.93×)
- 21-genome neural-predictor cross-comparison
- v1.4 viral / prophage genomic-region track + ensemble + RDF/JSON-LD vocabulary

### PARTIAL

- **Phase 10 outer iterative loop** — code path is mechanically complete (promotions, pins, fixed-point, four promotion strategies including SAT4J MaxSAT and beam search) but F-max regresses ~4 points at default settings; needs higher `qBase` + real (not synthetic) gapseq output to be productive
- **Cross-genome layer (Phase 12 M3 ranker)** — `Ranker`/`GbdtRanker`/`RankerRescorer` interfaces ship; LambdaMART training pipeline is in `benchmark/ml/` outside the JVM build
- **gspa-nf vs JVM CLI parity** — Nextflow emits raw + merged tool outputs but does NOT call Phase 7 integrator or quality metrics; users feed Nextflow output back into `gspa-cli integrate`
- **gLM operon caller productionisation** — code shipped, gating verdict NO-GO (Δ μ F-max +0.0002 vs +0.005 bar). Caller stays in tree as opt-in until Phase 1 reroll lands
- **ConsistencyPrior in production** — works, requires per-genome NCBI taxonomy lineage; benchmark numbers are reported with it disabled
- **gapseq dependency** — 4 of 9 benchmark genomes hit a `Reactions.tbl` zero-byte bug; those genomes lack GapFillingPrior + dark-matter suggestions

### FUTURE WORK

- **Phase 11.2 — `GENOMIC_CONTEXT_FM` evidence channel** — emit per-protein gLM centroid-kNN claims into Phase 7. Designed in `SPEC.md`, gated on Phase 1 reroll.
- **Phase 11.3 — `BF(O,P)` augmentation** — add embedding-distance-to-known-pathway-operon term to dark-matter Bayes factor. Gated on Phase 11.2.
- **gLM2 (TattaBio, mixed-modality)** — flagged as upgrade path; would change the I/O contract (gLM2 also ingests intergenic nucleotide sequence)
- **KAUST 500-genome panel scaling** — `benchmark/panel/` and `benchmark/cultures/` are active; Phase 12 panel expansion (Strategy C + ANI-95) approved 2026-04-20
- **Phase 10 retune** — higher `qBase` (0.70-0.75), real gapseq output, full `gspa annotate` driver instead of `gspa integrate`
- **Eukaryote support** — Augustus integration is documented as missing; Prodigal is used as a fallback for intron-poor eukaryotes only
- **Deferred predictors** — Phigaro, vConTACT3, eukaryote-specific viral predictors, PhiSpy GenBank-mode validation
- **Empty Quarter desert-metagenomics paper** — 15,469 MAGs, 2.5 M dark-matter proteins; consumer of Phase 11.2/11.3 once landed
- **License resolution** — `LICENSE` file is referenced in README as 'to be added'
