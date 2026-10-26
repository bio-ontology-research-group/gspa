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

## GSPA v3.0.0 Empirical Benchmark Suites & Biological Validation

All v3.0.0 benchmarks run through the shared production inference engine (`gspa integrate` / `IntegrationWorkflow`) and sidecar engines with SHA-256-recorded inputs, manifests, and stage ledgers. Full reports, case-level accounting, and honest positive and negative findings are documented in [`docs/benchmarks/`](docs/benchmarks/) and [`benchmark/v3/results/`](benchmark/v3/results/):

1. **Family-Masked Learned Context Integration, 10-Arm Ablation & Structural Perturbation Controls (`swissprot_label_rich_v1`, 18 Prokaryotic Genomes, 4,262 candidates):**
   - Out-of-fold / held-out family-disjoint evaluation (`train` = 12 genomes, `calibration` = 3, `test` = 3; `family_disjoint_verified = true`).
   - Production logistic model (`model_logistic.json`, held-out AUROC = `0.7934`, Brier = `0.0997`, ECE = `0.0577`) improves held-out GO $F_{\max}$ under homology claims from `0.2881` (`stage0`) to **`0.3221`** (`+0.0340`, precision rising from `0.2733` to `0.3529`) and under family-excluded claims from `0.1841` to **`0.1921`**, while GBDT (`model_gbdt.json`) lifts family-excluded EC top-1 from `0.1081` to **`0.1351`**.
   - Perturbation controls verify biological mechanism: shuffling gene neighborhoods drops GO $F_{\max}$ by `-0.0111` (precision `0.3529` $\to$ `0.2878`), $50\%$ assembly fragmentation drops GO $F_{\max}$ by `-0.0177` and EC top-1 by `-0.0541`, and $5\times$ strain duplication leaves GO $F_{\max}$ (`0.3221`) and EC top-1 (`0.3784`) invariant (`max_abs_probability_diff = 0.0017`).
   - Report: [`docs/benchmarks/v3-training-cohort.md`](docs/benchmarks/v3-training-cohort.md) | [`docs/benchmarks/v3-training-data.md`](docs/benchmarks/v3-training-data.md).
2. **Full-Proteome Price-148 RB-TnSeq Context, Cohort, Learned-Model & Locus-Decoy Benchmark (22 NCBI Assemblies, 115,542 Proteins, 42 Production Runs):**
   - Evaluates 7 predictor variants (`clean`, `dpz`, `enzbert`, `graphec`, `mapred`, `topec`, `ensemble_pool4`) across 6 inference arms (`stage0`, `stage1`, `stage1_rlc`, `stage1_rlc_cohort`, `learned_logistic`, `learned_gbdt`).
   - On isolated 4-digit EC numbers lacking curated `ec2go` mappings, heuristic `stage1` priors leave raw EC rankings unchanged, whereas `learned_logistic` / `learned_gbdt` revise candidate probabilities directly (`ensemble_pool4` top-1 = `61/148`, `0.4122`; `graphec` top-1 improves from `18/148` to `23/148`, `+3.38%`) and `stage1_rlc_cohort` generates `5,693–6,942` stoichiometric gap-rescue EC suggestions per variant.
   - Report: [`docs/benchmarks/v3-price-development.md`](docs/benchmarks/v3-price-development.md) | [`docs/benchmarks/v3-price-context-data.md`](docs/benchmarks/v3-price-context-data.md).
3. **Metabolic Knockout & Auxotroph Phenotype Validation Across 3 Registered Protocols (61 Measured Cases, `iML1515`):**
   - Evaluates `legacy_code_defaults`, `manuscript_unregulated`, and `manuscript_regulatory_closures` (with independent COBRApy/GLPK cross-verification, max growth difference $< 10^{-9}$). Regulatory reaction closures (`R_TPI` on D-alanine) rescue `tpiA5` (`25/44` $\to$ `26/44` binary agreement, `59.1%`).
   - Report: [`docs/benchmarks/v3-metabolic-development.md`](docs/benchmarks/v3-metabolic-development.md).
4. **Coupled SBML Community Cross-Feeding Benchmark (Wintermute et al. 2010, 1,035 Synthetic Lethal Pairs × 6 Regimes = 6,210 Coupled LPs):**
   - Evaluates joint two-strain stoichiometric FBA against matched no-transfer controls across 46 *E. coli* Keio auxotrophs (`1,035` pairs). Under trace-rescue aerobic conditions, best-F1 reaches `0.340` (`precision = 0.219`, `recall = 0.765`, `AUROC = 0.562`), quantifying the structural limits of mass-balance stoichiometry without metabolite-specific exporter kinetics.
   - Report: [`docs/benchmarks/v3-community-development.md`](docs/benchmarks/v3-community-development.md).
5. **Stage 3 Anomaly Resolution: `qwen3.8-27b` vs. Equally Informed Simple Models (40 Held-Out Test Anomalies):**
   - Compares `qwen3.8-27b` (`http://unimatrix01.kaust.edu.sa:8000/v1`, with and without symbolic guardrails) against 6 deterministic and learned comparators on identical held-out anomaly packets. Learned classifiers (`learned_gbdt` net gain `+23`, precision `0.8966`; `learned_logistic` net gain `+22`, precision `0.8667`) outperform `qwen3.8-27b` (`qwen_agent_no_guardrail` net gain `+15`, precision `1.0000`; `qwen_agent_guarded` net gain `+7`, precision `1.0000`) at $10^6\times$ lower latency (`0.0001 s` vs. `131.6 s`).
   - Report: [`docs/benchmarks/v3-stage3-agent-evaluation.md`](docs/benchmarks/v3-stage3-agent-evaluation.md).
6. **Equal-Budget Condition-Planner Benchmark (`qwen3.8-27b` vs. 5 Deterministic Planners, 44 Genes × 24 Carbon Sources):**
   - Compares `fixed_paired`, `random_paired`, `evidence_paired`, `random_item`, `evidence_item`, and `qwen` (`qwen3.8-27b`) under identical solver budgets ($B \in \{8, 16\}$) alongside a secondary metabolism (`MIBiG` / Westhoff 2021) condition-dependence audit.
   - At budget $B=8$ (LP cap 9), strict Qwen-attributable discovery reaches **`85.4%`** (`35/41` eligible genes, `100.0%` on clean runs) vs. `39.0%` (`16/41`) for `evidence_paired` (exact McNemar $p = 0.0000$), `19.5%` for `random_paired`, `17.1%` for `evidence_item`, and `12.2%` for `fixed_paired`. At budget $B=16$ (LP cap 17), `qwen` achieves **`80.5%`** strict Qwen-attributable (`87.8%` hybrid) vs. `75.6%` for `evidence_paired`, and **`75.9%`** (`22/29`) binary agreement on visited measured pairs.
   - Report: [`docs/benchmarks/v3-planner-evaluation.md`](docs/benchmarks/v3-planner-evaluation.md).

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
