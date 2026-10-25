# GSPA v2.1.0 — Genome-Scale Protein Annotation

[![Release](https://img.shields.io/badge/release-v2.1.0-blue.svg)](https://github.com/bio-ontology-research-group/gspa/releases/tag/v2.1.0)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![Java 21+](https://img.shields.io/badge/Java-21%2B-orange.svg)](https://adoptium.net/)

**GSPA (Genome-Scale Protein Annotation)** is a stand-alone, three-stage neuro-symbolic and flux-certified stoichiometric framework for whole-proteome protein function prediction and genome-scale quality assessment across **Bacteria**, **Archaea**, and **Eukaryota**.

Whereas standard annotation pipelines and protein language model (PLM) predictors annotate each protein in isolation, **GSPA v2.1.0** formulates function prediction as a **joint genome-wide inference problem**:
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

### 1. Build or Download Stand-Alone Release (`v2.1.0`)

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

## Empirical Benchmark Results Across 6 Pillars & 8 Suites (`GSPA v2.1.0`)

Full tables, immutable manifests (`MANIFEST_CHECKSUMS.sha256`), and reports are in [`results/`](results/) and [`benchmark/manifests/`](benchmark/manifests/).

### 1. 108-Genome Stoichiometric Audit & Re-Audited `Price-149` RB-TnSeq Benchmark (`Pillar 2`, `Suite 6`)

- **108-Genome Mathematical Invariant Audit:** Negative-cost MILP gap-filling (`METEOR`) adds a mean of `313.35` reactions/genome of which `311.19` (**99.31%**) carry zero flux ($|v_j| \le 10^{-7}$). GSPA's positive-cost + lexicographic pFBA `StoichiometricLocusRefiner` adds `8.05` reactions/genome, **100.0%** (`8.05/8.05`) of which are flux-certified (**0.0% padding**).

| Baseline Predictor (`Price-149`, $n=149$) | Baseline Top-1 (Acc.) | Baseline $F_{\max}$ | `+ METEOR` $F_{\max}$ | `+ SLR (Cognate)` $F_{\max}$ | `+ GSPA Full` Top-1 (Acc.) | `+ GSPA Full` $F_{\max}$ ($\Delta$) |
|---|---:|---:|---:|---:|---:|---:|
| **CLEAN** | 60 (40.3%) | 0.5592 | 0.5635 | 0.6132 | **81 (54.4%)** | **0.6680 (`+0.1088`)** |
| **DeepProZyme** | 46 (30.9%) | 0.5506 | 0.5689 | 0.6186 | **71 (47.7%)** | **0.6734 (`+0.1228`)** |
| **EnzBERT** | 53 (35.6%) | 0.4773 | 0.4805 | 0.5302 | **74 (49.7%)** | **0.5850 (`+0.1077`)** |
| **MAPred** | 50 (33.6%) | 0.5497 | 0.5480 | 0.5994 | **70 (47.0%)** | **0.6542 (`+0.1045`)** |
| **GraphEC** | 18 (12.1%) | 0.2640 | 0.3402 | 0.3899 | **42 (28.2%)** | **0.4447 (`+0.1807`)** |
| **TopEC** | 5 (3.4%) | 0.0336 | 0.0336 | 0.0403 | **10 (6.7%)** | **0.0588 (`+0.0252`)** |

### 2. Large-Scale `BacDive` Experimental Phenotype & Enzyme Validation (`Pillar 3`, `Suites 7 & 8`, 4,788 Strains)

| Evaluation Track / Panel (`BacDive`) | Taxonomic Split / Scope | Assays ($n$) | DPZ Baseline | METEOR (MILP) | GSPA Stage 2B (SLR) | GSPA Full (`v2.1.0`) |
|---|---|---:|---:|---:|---:|---:|
| **10 Diagnostic Enzyme Assays ($F_1$)** | All 10 EC classes pooled | 99,672 | 0.6443 | 0.6133 | 0.7186 | **0.7342** |
| **Heteromeric Urease `EC 3.5.1.5` ($F_1$)** | UreABG ($\alpha\beta\gamma$) AND-complex | 9,965 | 0.5748 | 0.5610 | 0.6716 | **0.7051** |
| **`Track 1`: Blind Substrate Growth (AUC)** | 4,788 Strains (Genus-Blocked) | 72,312 | 0.5540 | 0.5310 | 0.7180 | **0.7460** |
| **`Track 1`: Blind Substrate Growth (AUC)** | 1,000 RefSeq Panel (Genus-Blocked) | 15,840 | 0.5650 | 0.5440 | 0.7290 | **0.7580** |
| **`Track 2`: 30%-Assisted Growth (AUC)** | Held-out 70% Assays (Genus-Blocked) | 11,088 | 0.5650 | 0.5720 | 0.7710 | **0.7940** |
| **`Track 2`: Dark-Matter Gene Attribution** | Catabolic Locus Hit@1 (Genus-Blocked) | 1,000 | 0.2840 | 0.2840 | 0.6150 | **0.7480** |

### 3. Role-Stratified LRO (`Suites 3–5`), Heteromeric Complexes, Eukaryotic Anatomy (`Pillar 4`) & Stage 3 Agent (`Pillar 5`)

- **15-Genome Role-Stratified LRO (`Suites 3–5`):** In the $\le 30\%$ identity dark proteome (`Protocol H`), where sequence/PLM homology and METEOR achieve `0.0%` Hit@1, GSPA achieves **100.0% Hit@1** on `Stratum A` (Core Biomass) and **69.2% Hit@1** (`100.0%` conditioned on $\mathcal{C}_{\text{unique}}=1$) on `Stratum B2` (Cognate Substrate Catabolic Flux).
- **Obligate Heteromeric Protein Complexes (175 LSO Trials across 6 Prokaryotes):** Raises held-out subunit recovery $F_1$@0.50 from **0.000** (`Noisy-OR`) to **0.742** (`Stage 1 + Stage 2 + Stage 3 Agent`, **100.0% Recall**, **98.3% Hit@1**).
- **Multicellular Eukaryotic Anatomical Constraints (6 Species):** Eliminates **100.0%** ($2{,}377/2{,}377$) of cross-clade anatomical violations vs. **4.4%** for direct GO taxon constraints, while semantic ancestor relaxation (`relax`) improves CAFA $F_{\max}$ over `remove` across all 6 species (`+0.009` to `+0.029`).
- **Stage 3 Agent (`gspa-agent`, 120 Anomalies):** Outperforms 5-fold CV non-agentic Logistic Regression and GBDT rerankers (`F1 = 0.750, Precision = 0.600`) with **F1 = 0.939, Precision = 0.885, Recall = 1.000, Accuracy = 94.2%**, and **0.0% ontology violation rate** (`0/120` vs. `6.7%` unguarded).

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
