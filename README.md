# GSPA v3.0.0 — Genome-Scale Protein Annotation

[![Release](https://img.shields.io/badge/release-v3.0.0-blue.svg)](https://github.com/bio-ontology-research-group/gspa/releases/tag/v3.0.0)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![Java 21+](https://img.shields.io/badge/Java-21%2B-orange.svg)](https://adoptium.net/)

**GSPA (Genome-Scale Protein Annotation)** is a stand-alone neuro-symbolic and flux-certified stoichiometric framework for whole-proteome protein function prediction, multi-genome cohort context inference, and genome-scale quality assessment across **Bacteria**, **Archaea**, and **Eukaryota**.

Whereas standard annotation pipelines and protein language model (PLM) predictors annotate each protein in isolation, **GSPA v3.0.0** formulates function prediction as a **joint genome- and cohort-wide inference problem** executed through one shared `IntegrationWorkflow` across CLI, batch, web portal, and Nextflow routes:
1. **Stage 1 — Iterative Bayesian Evidence Integration**: Combines calibrated per-protein evidence across sequence homology (`DIAMOND`, `MMseqs2`), domain profiles (`Pfam`, `InterProScan`), structural similarity (`FoldSeek`, `DeepFRI`), and deep protein language models (`ESM-2`, `ProtT5`, `DeepGOWeb`, `ProtInfer`, `CLEAN`) via group-wise Noisy-OR pooling, then iteratively refines log-odds posteriors across the entire proteome using **six genome-scale priors**:
   - `EssentialityPrior`: Core domain-universal essential function completeness (`bacteria`, `archaea`, `eukaryote`).
   - `CoherencePrior`: Metabolic pathway/module stoichiometry (KEGG, MetaCyc) and obligate heteromeric protein complex completion (`complex_terms.tsv`).
   - `ConsistencyPrior`: SAT4J UNSAT-core taxonomic consistency (`only_in_taxon`, `never_in_taxon`).
   - `AnatomicalPrior`: Cross-ontology GO–Uberon–Plant Ontology (`PO`)–Cell Ontology (`CL`) anatomical and developmental structure constraints (`2,714` verified logical definitions linked to `8,330` taxon constraints).
   - `GapFillingPrior` & `RLGC S2`: Metabolic network gap-filling (`gapseq`) coupled to operon pathway hypergeometric Bayes factors, Pfam coarse EC compatibility, and `ESM-2` reaction-centroid similarity.
   - `GenomicContextPrior`: Operon and synteny neighborhood functional coupling.
2. **Stage 2A & 2B — Deterministic Symbolic Enforcement & Flux-Certified Stoichiometric Locus Refinement (`StoichiometricLocusRefiner`)**:
   - **Stage 2A (`--enforce-stage2`)**: Guarantees $100\%$ taxonomic and anatomical consistency via **semantic ancestor relaxation** (`--enforce-consistency-mode relax`) and prunes unsupported singleton subunits of obligate heteromeric complexes.
   - **Stage 2B (`--enable-fba-refiner`)**: Compiles Gene–Protein–Reaction (GPR) mechanisms with `SparseRoleAggregator` (Top-$M$ isozyme OR and geometric-mean subunit AND), solves positive-cost flux-certified network repair with two-step lexicographic parsimonious FBA (pFBA) in HiGHS (`0.0%` zero-flux padding vs. `99.3%` in negative-cost MILPs), and redistributes flux support back to specific genomic loci via `LocusResponsibilityUpdater` against an explicit unassigned-gene hypothesis $h_\varnothing$.
3. **Stage 3 — Optional Neuro-Symbolic Agentic Annotation Completion (`gspa-agent`, `--enable-agent`)**: Selectively routes remaining post-Stage-2 genome-scale anomalies (`PATHWAY_GAP`, `INCOMPLETE_COMPLEX`, `BORDERLINE_ESSENTIAL`, `BORDERLINE_CONFLICT` — typically $<2.5\%$ of the proteome) to a tool-using local/cluster language model (`Qwen 3.8 27B` or fast 1-step SLM mode) gated by a deterministic **Symbolic Guardrail** that blocks any hallucinated or clade-violating GO terms.

📖 **Full CLI & Stand-Alone User Guide:** [`docs/CLI_REFERENCE.md`](docs/CLI_REFERENCE.md)  
📄 **Manuscript & Figures:** [`paper/main.pdf`](paper/main.pdf) | [`results/`](results/) | [`benchmark/manifests/`](benchmark/manifests/)

---

## Quickstart (Stand-Alone Tool)

### 1. Build or Download Stand-Alone Release (`v3.0.0`)

```bash
git clone https://github.com/bio-ontology-research-group/gspa.git
cd gspa
./gradlew clean test :gspa-cli:shadowJar
./gspa --version
```

### Choose a workflow preset

```bash
./gspa presets
./gspa integrate --preset bacteria --claims claims.jsonl --out integrated.tsv
./gspa annotate --preset eukaryote --proteins proteome.faa --db reference.dmnd -o out/
./gspa annotate --preset metagenomics --print-config
```

Presets cover prokaryotes, bacteria, archaea, eukaryotes, individual MAGs and
mixed metagenomes. Explicit command-line options override their defaults.
Sequence annotation needs installed predictors and their data; integration can
use the bundled references. See [workflow presets](docs/scenarios/presets.md)
for setup, biological assumptions and overrides.

### Download references and models

```bash
./gspa resources
./gspa download --group annotation --dir /data/gspa
./gspa download deepgo-plusplus-models --dir /data/gspa
./gspa download --all --dry-run
# To attempt every optional resource as well:
./gspa download --all --dir /data/gspa
```

Individual IDs and groups use one catalog. Downloads write reusable paths and
integrity receipts; unavailable or licensed assets are reported explicitly.
See [resource downloads](docs/reference/resources.md) for prerequisites,
large-download sizes, model coverage and custom bundles.

### 2. Run Genome-Scale Integration (`gspa integrate`)

```bash
# Minimal stand-alone run (using bundled taxonomic, anatomical, complex, and essentiality resources)
./gspa integrate \
  --claims claims.jsonl \
  --out integrated.tsv \
  --taxon NCBITaxon_562 \
  --essential-profile bacteria \
  --enforce-stage2 \
  --enforce-consistency-mode relax \
  --lite

# Full run with Stage 1 Priors, Stage 2A Symbolic Enforcement, Stage 2B StoichiometricLocusRefiner, and Stage 3 Agent
./gspa integrate \
  --claims ecoli_claims.jsonl \
  --out ecoli_integrated.tsv \
  --go-owl go.owl \
  --taxon NCBITaxon_562 \
  --essential-profile bacteria \
  --pathways kegg_pathways.tsv \
  --modules kegg_modules.tsv \
  --ec2go ec2go.txt \
  --operons ecoli_operons.tsv \
  --gaps ecoli_gaps.jsonl \
  --enable-priors essentiality,coherence,consistency,anatomical,gap_filling,genomic_context \
  --enforce-stage2 \
  --enforce-consistency-mode relax \
  --enable-fba-refiner \
  --fba-data-dir benchmark/manifests/reference_data/psb_revision/code_snapshot/data \
  --fba-gram neg --fba-medium minimal \
  --dark-matter --suggestions-out ecoli_dark_matter.tsv \
  --enable-agent --agent-model qwen3.8-27b \
  --lite
```

---

## v3 development evidence and outstanding validation

The follow-up audit reopened all twelve acceptance requirements after the
`v3.0.0` tag. That tag does not establish completion of the requested biological
validation. The [current acceptance record](tasks/v3-acceptance.md) distinguishes
implemented mechanisms, verified execution, empirical results and remaining work.
AllTheBacteria is excluded from this evaluation.

- **Learned context:** production logistic and tree models have been fitted and
  compared through Python and JVM inference. Their present target is reference
  annotation support under an explicit closed-world policy. Unconfirmed
  alternatives are not experimentally proven negatives; current models have no
  trained locus/null layer. Predictor identity and input-range checks make
  unsupported inputs abstain. The [family-excluded results](docs/benchmarks/v3-family-excluded-training.md)
  and [observed-gap comparison](docs/benchmarks/v3-observed-gap-cohort.md) retain
  low precision, harms, abstentions and negative comparisons.
- **EC prediction:** the matched Price development run found no top-1 gain from
  current pathway priors: CLEAN and the four-predictor pool each retain 60/148
  correct targets. The historical transfer of homology-trained context models
  to EC-specialist scores degraded ranking substantially. Those results do not
  support a large improvement over CLEAN. See the [complete Price reports](docs/benchmarks/v3-price-development.md)
  and the [input-coverage audit](docs/benchmarks/v3-price-context-data.md).
  The [zero-baseline audit](docs/benchmarks/v3-price-calibration-audit.md) also shows
  why lower Brier/ECE alone does not establish better enzyme recovery.
- **Reaction context and recurrence:** matched production comparisons now include
  observed search gaps, whole-cohort context and explicit pathway, universal-network
  and organism-specific modes. Verified execution and duplicate/sequence-loss
  controls do not by themselves establish accurate missing-enzyme recovery.
  The [completed context comparison](docs/benchmarks/v3-context-comparison-results.md)
  finds no consistent context advantage in its controlled homology-source arms.
  Only 15/39 experimental-reference EC functions enter that candidate set, and
  5/39 after family exclusion. The [complete candidate trace](docs/benchmarks/v3-ec-candidate-coverage.md)
  locates the losses before integration; routine reports now retain missing
  candidates and both cutoffs.
- **Metabolism and communities:** published auxotroph, nutrient and coculture
  observations are evaluated separately from independent solver agreement.
  Several phenotype comparisons underperform constant baselines. Growth rescue
  is conditional on model, medium and transport assumptions and is not a gene
  assignment. See the [auxotroph](docs/benchmarks/v3-metabolic-development.md),
  [nutrient phenotype](docs/benchmarks/v3-nutrient-phenotypes.md) and
  [community](docs/benchmarks/v3-nutrient-community-results.md) reports.
- **Agent comparison:** the initial Stage 3 benchmark was invalidated because
  public inputs leaked answers and included invented learned probabilities.
  Its reported gains are not accepted evidence. The [replacement production comparison](docs/benchmarks/v3-stage3-production-comparison.md)
  executes real unresolved cases, but lacks matched experimental outcome labels
  and therefore cannot yet establish agent accuracy.
- **Condition planning:** the [equal-budget planner experiment](docs/benchmarks/v3-planner-evaluation.md)
  tests condition selection for known assignments. Independent function recovery
  with every target assignment masked, calibrated null/role decisions, and
  measured secondary-product validation remain outstanding. The [masked-assignment implementation](docs/benchmarks/v3-masked-assignment-experiments.md)
  records current mechanisms and biological limitations.

---

## Measured v3 results (release-candidate benchmarks, 2026-10-31)

The v2.1.0 results tables that stood here are withdrawn. The Price-149 "+ GSPA Full" gain (about +0.10 F-max)
came from a stand-alone script that probed each protein's own assay substrate, which leaks the answer; the BacDive,
Protocol H and Stage 3 agent figures could not be traced to runs or results files
(`paper/validation_report_fable_2026-10-19.md`). The table below lists what the current code does on established
benchmarks; every number links to a report with its inputs, receipts and negative results.

| Benchmark | Result | Report |
|---|---|---|
| Six enzyme predictors on the 148 Price et al. (2018) experimental enzymes (22 bacteria) | Top-1 correct: CLEAN 60, EnzBERT 53, MAPred 50, DeepProZyme 46, GraphEC 18, TopEC 5; GSPA integration without the metabolic step changes none of them | [Price evaluation](docs/benchmarks/v3-price-rc-evaluation.md) |
| Stage 2B flux-certified refinement on the same proteins | No predictor improves (0 corrections; CLEAN top-3 87 → 80); it only repairs the few reactions needed for growth on one medium | [Stage 2B](docs/benchmarks/v3-price-stage2b.md) |
| Six pre-registered follow-ups to Stage 2B | No gain on the Price enzymes; a top-10 candidate limit (`--fba-candidate-top-k 10`, off by default) removes Stage 2B's one harm; gene assignments are at chance against RB-TnSeq fitness | [Follow-ups](docs/benchmarks/v3-stage2b-followups.md) |
| Taxon consistency for enzyme claims (`--taxon-lineage`) | Calls forbidden in bacteria above threshold: 465 → 0 across six predictors; no Price call lost | [Taxon fix](docs/benchmarks/v3-price-taxon-fix.md) |
| Biology checks of real predictions (tRNA charging, core enzymes, taxon consistency, truth agreement) | Core enzymes found in 21/22 genomes for four predictors; most charging failures are the GatCAB asparagine label | [Biology checks](docs/benchmarks/v3-biology-checks.md) |
| Off-growth enzyme discovery (pre-registered; targeted gap-filling) | 23 of 178 hidden experimental enzymes recovered vs 0 for biomass gap-filling; most calls on labelled proteins wrong | [Off-growth](docs/benchmarks/v3-offgrowth-discovery.md) |

The per-requirement verdicts are in the [acceptance record](tasks/v3-acceptance.md).

---

## License

GSPA is licensed under the **GNU General Public License v3.0 or later (`GPL-3.0-or-later`)**. See [`LICENSE`](LICENSE) for full details.


### Biological tutorials

Start with the [biological tutorial collection](docs/tutorials/index.md): identify
candidate threonine-biosynthesis genes, investigate transferred cellular locations,
and learn how missing sequence or mixed organisms change the interpretation.

```bash
./tutorial setup
./tutorial run bacteria
./tutorial browser bacteria
```

Short lessons cover bacteria, archaea, yeast, incomplete bins, mixed communities,
phage, and unknown prokaryotic domains. Each connects observed results with a
biological question and a worked interpretation. An optional
[neural lesson](docs/tutorials/neural.md) explores ESM-2 candidates on CPU and GPU.
