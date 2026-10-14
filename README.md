# GSPA v2.0.0 — Genome-Scale Protein Annotation

[![Release](https://img.shields.io/badge/release-v2.0.0-blue.svg)](https://github.com/bio-ontology-research-group/gspa/releases/tag/v2.0.0)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![Java 21+](https://img.shields.io/badge/Java-21%2B-orange.svg)](https://adoptium.net/)

**GSPA (Genome-Scale Protein Annotation)** is a stand-alone, three-stage neuro-symbolic framework for whole-proteome protein function prediction and genome-scale quality assessment across **Bacteria**, **Archaea**, and **Eukaryota**.

Whereas standard annotation pipelines and protein language model (PLM) predictors annotate each protein in isolation, **GSPA v2.0.0** formulates function prediction as a **joint genome-wide inference problem**:
1. **Stage 1 — Iterative Bayesian Evidence Integration**: Combines calibrated per-protein evidence across sequence homology (`DIAMOND`, `MMseqs2`), domain profiles (`Pfam`, `InterProScan`), structural similarity (`FoldSeek`, `DeepFRI`), and deep protein language models (`ESM-2`, `ProtT5`, `DeepGOWeb`, `ProtInfer`, `CLEAN`) via group-wise Noisy-OR pooling, then iteratively refines log-odds posteriors across the entire proteome using **six genome-scale priors**:
   - `EssentialityPrior`: Core domain-universal essential function completeness (`bacteria`, `archaea`, `eukaryote`).
   - `CoherencePrior`: Metabolic pathway/module stoichiometry (KEGG, MetaCyc) and obligate heteromeric protein complex completion (`complex_terms.tsv`).
   - `ConsistencyPrior`: SAT4J UNSAT-core taxonomic consistency (`only_in_taxon`, `never_in_taxon`).
   - `AnatomicalPrior`: Cross-ontology GO–Uberon–Plant Ontology (`PO`)–Cell Ontology (`CL`) anatomical and developmental structure constraints (`2,714` verified logical definitions linked to `8,330` taxon constraints).
   - `GapFillingPrior` & `RLGC S2`: Metabolic network gap-filling (`gapseq`) coupled to operon pathway hypergeometric Bayes factors, Pfam coarse EC compatibility, and `ESM-2` reaction-centroid similarity.
   - `GenomicContextPrior`: Operon and synteny neighborhood functional coupling.
2. **Stage 2 — Deterministic Symbolic Constraint & Coherence Enforcement**: Guarantees $100\%$ taxonomic and anatomical consistency via **semantic ancestor relaxation** (`--enforce-consistency-mode relax`, replacing clade-invalid specific GO terms with their most specific valid GO ancestor while preserving unpenalized evidence log-odds) and prunes unsupported singleton subunits of obligate heteromeric complexes.
3. **Stage 3 — Optional Neuro-Symbolic Agentic Annotation Completion (`gspa-agent`, `--enable-agent`)**: Selectively routes remaining post-Stage-2 genome-scale anomalies (`PATHWAY_GAP`, `INCOMPLETE_COMPLEX`, `BORDERLINE_ESSENTIAL`, `BORDERLINE_CONFLICT` — typically $<2\%$ of the proteome) to a tool-using local/cluster language model (`Qwen 3.8 27B` or fast 1-step SLM mode) gated by a deterministic **Symbolic Guardrail** that blocks any hallucinated or clade-violating GO terms.

📖 **Full CLI & Stand-Alone User Guide:** [`docs/CLI_REFERENCE.md`](docs/CLI_REFERENCE.md)  
📄 **Manuscript & Figures:** [`paper/main.pdf`](paper/main.pdf) | [`results/`](results/)

---

## Quickstart (Stand-Alone Tool)

### 1. Download Stand-Alone Release (`v2.0.0`)

Download the pre-compiled stand-alone distribution (`gspa-2.0.0-standalone.tar.gz`) or the self-contained executable JAR (`gspa-2.0.0.jar`) directly from [GitHub Releases](https://github.com/bio-ontology-research-group/gspa/releases/tag/v2.0.0) (requires **Java 21+**):

```bash
curl -LO https://github.com/bio-ontology-research-group/gspa/releases/download/v2.0.0/gspa-2.0.0-standalone.tar.gz
tar -xzf gspa-2.0.0-standalone.tar.gz
cd gspa-2.0.0

./gspa --version
# 2.0.0
```

Or build from source with the included Gradle wrapper:

```bash
git clone https://github.com/bio-ontology-research-group/gspa.git
cd gspa
./gradlew clean test :gspa-cli:shadowJar
./gspa --version
```

### 2. Run Genome-Scale Integration (`gspa integrate`)

Because `gspa-2.0.0.jar` bundles all cross-ontology constraint tables (`go-taxon-constraints.tsv`, `ncbi-taxon-hierarchy.tsv`, `go-anatomical-constraints.tsv`, `uberon-po-taxon-constraints.tsv`, `uberon-po-hierarchy.tsv`, `complex_terms.tsv`) and essentiality presets (`bacteria`, `archaea`, `eukaryote`), you can run Stage 1 + Stage 2 immediately out of the box:

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

# Full run with GO OWL, KEGG pathways/modules, operons, metabolic gaps, dark-matter RLGC, and Stage 3 Agent
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
  --enforcement-actions-out ecoli_stage2_actions.tsv \
  --dark-matter --suggestions-out ecoli_dark_matter.tsv \
  --enable-agent --agent-model qwen3.8-27b \
  --lite
```

### 3. Run GAEF Genome-Scale Quality Evaluation (`gspa evaluate`)

```bash
# Convert integrated.tsv + proteome.faa to GFF3 + GAF2.2
uv run python benchmark/make_synth_gff_gaf.py \
  --fasta bsubtilis.faa --tsv integrated.tsv \
  --gff-out synth.gff --gaf-out synth.gaf

# Evaluate Completeness, Pathway Coherence, Complex Coherence, Taxon/Anatomy Consistency & IC
./gspa evaluate \
  -i synth.gff -a synth.gaf \
  --go-owl go.owl \
  --taxon NCBITaxon_1423 -k bacteria \
  --pathways kegg_pathways.tsv --modules kegg_modules.tsv --ec2go ec2go.txt \
  --lite -o quality_report.json
```

### 4. Interactive HTML5 Genome & Pathway Browser (`gspa visualize`)

```bash
./gspa visualize \
  --workdir . \
  --genome-id bsubtilis \
  --out bsubtilis_browser.html
```

---

## Empirical Benchmark Results (GSPA v2.0.0)

GSPA v2.0.0 is evaluated across **four benchmark pillars** spanning 21 reference genomes, 13 PGAP/DeepFRI benchmark genomes, 97 uncultured MAGs (16 phyla), 1,248 leave-reaction-out (LRO) metabolic gap-filling trials, 175 leave-subunit-out (LSO) heteromeric complex trials, 6 eukaryotic anatomical constraint benchmarks, and a 120-anomaly Stage 3 agent benchmark. Full tables and reports are in [`results/`](results/).

### Pillar 1 & 2: Multi-Predictor Integration, 9-Genome Ablation & Cross-Clade Perturbation Defense

- **21-Genome 7-Predictor Benchmark:** Integrating DIAMOND, InterProScan, Pfam, ESM-2, ProtT5, DeepGOWeb, and ProtInfer via GSPA achieves **0.869 micro $F_{\max}$** (**0.822 CAFA $F_{\max}$**), outperforming every individual predictor (`ESM-2`: `0.858`, `ProtT5`: `0.853`, `DIAMOND`: `0.840`, `ProtInfer`: `0.751`, `DeepGOWeb`: `0.692`).
- **13-Genome PGAP / DeepFRI Comparison:** 3-source GSPA (`DIAMOND` + `InterProScan` + `Pfam`) achieves **0.818 micro $F_{\max}$** (**0.806 CAFA $F_{\max}$**), surpassing **NCBI PGAP** (`0.782` micro $F_{\max}$, `0.768` CAFA $F_{\max}`), **InterProScan** (`0.722`), **Pfam** (`0.506`), and **metagenomic-deepFRI** (`0.585`).

| Ablation Configuration (9 Reference Genomes, 81 Runs) | Micro $F_{\max}$ | CAFA $F_{\max}$ | Essential Completeness | Pathway Coherence | Complex Coherence | GAEF Composite Score | Perturbation Suppression (900 Injected Violations) |
|---|---:|---:|---:|---:|---:|---:|---:|
| `A1_noisy_or` (Unconstrained Baseline) | 0.8705 | 0.8407 | 0.9400 | 0.4975 | 0.8636 | 0.8459 | 0.0% |
| `A2_essentiality` (`+ EssentialityPrior`) | 0.8705 | 0.8407 | 0.9447 | 0.4975 | 0.8636 | 0.8473 | 0.0% |
| `A3_coherence` (`+ CoherencePrior`) | 0.8707 | 0.8406 | 0.9400 | 0.5146 | 0.8636 | 0.8475 | 0.0% |
| `A6_consistency_anatomical` (`+ Taxon & Anatomy`) | **0.8708** | 0.8405 | 0.9400 | 0.4975 | 0.8636 | 0.8542 | **99.7%** |
| `A7_stage1_all` (All 6 Stage 1 Priors) | 0.8705 | 0.8398 | **0.9479** | **0.5146** | **0.8852** | **0.8643** | 99.7% |
| `A9_stage2_relax` (Stage 1 + Stage 2 `relax`) | 0.8705 | 0.8398 | **0.9479** | **0.5146** | **0.8852** | **0.8643** | **99.7%** (`100%` direct taxon) |

*(Full report: [`results/pillar1_2_ablation_report.md`](results/pillar1_2_ablation_report.md))*

### Pillar 3, Scenario A: Leave-Reaction-Out (LRO) Metabolic Gap-Filling & 97-MAG Dark Matter

| Benchmark Dataset | Method / Scoring Engine | Cases ($n$) | Hit@1 | Hit@3 | Hit@5 | MRR |
|---|---|---:|---:|---:|---:|---:|
| ***E. coli* K-12 All LRO** | Baseline Homology Suggester | 850 | 1.8% | 5.4% | 8.9% | 0.051 |
| ***E. coli* K-12 All LRO** | GSPA RLGC S2 (Catalytic + Operon BF + Pfam) | 850 | 63.5% | 85.2% | 91.6% | 0.758 |
| ***E. coli* K-12 All LRO** | **GSPA RLGC S2 + ESM-2 Reaction Centroid** | 850 | **63.7%** | **85.5%** | **91.8%** | **0.760** |
| ***E. coli* K-12 Twilight Zone ($<25\%$ id.)** | Baseline Homology Suggester | 80 | 5.0% | 11.3% | 17.5% | 0.104 |
| ***E. coli* K-12 Twilight Zone ($<25\%$ id.)** | GSPA RLGC S2 (Catalytic + Operon BF + Pfam) | 80 | 25.0% | 53.8% | 72.5% | 0.442 |
| ***E. coli* K-12 Twilight Zone ($<25\%$ id.)** | **GSPA RLGC S2 + ESM-2 Reaction Centroid** | 80 | **30.0%** | **56.3%** | **72.5%** | **0.478** |
| ***B. subtilis* 168 Cross-Species LRO** | Baseline Homology Suggester | 398 | 2.5% | 7.3% | 11.0% | 0.064 |
| ***B. subtilis* 168 Cross-Species LRO** | **GSPA RLGC S2 (Default Parameters)** | 398 | **25.1%** | **51.8%** | **68.8%** | **0.431** |
| **Ambiguous Multi-Candidate LRO** | Stage 1 RLGC S2 Alone | 50 | 52.0% | 82.0% | 100.0% | 0.701 |
| **Ambiguous Multi-Candidate LRO** | **Stage 1 RLGC S2 + Stage 3 `gspa-agent`** | 50 | **86.0%** | **96.0%** | **100.0%** | **0.915** |

- **97 Uncultured Prokaryotic MAGs (16 Phyla):** Resolves **2,175** missing metabolic pathway steps at $q \ge 0.30$ (`22.4` gaps/MAG) and **568** high-confidence pathway steps at $q \ge 0.50$ (`5.9` gaps/MAG).

### Pillar 3, Scenario B & C: Heteromeric Complex Completion & Eukaryotic Anatomical Constraints

- **Obligate Heteromeric Protein Complexes (175 LSO Trials across 6 Prokaryotes):** Raises held-out subunit recovery $F_1$@0.50 from **0.000** (`Noisy-OR`) to **0.668** (`Stage 1 Coherence + Operon Context`, `96.0%` Hit@1) and **0.742** (`+ Stage 2 Pruning + Stage 3 Agent`, **100.0% Recall**, **98.3% Hit@1**).
- **Multicellular Eukaryotic Anatomical Constraints (6 Species: *S. cerevisiae*, *C. elegans*, *D. melanogaster*, *P. patens*, *A. thaliana*, *M. musculus*):** Direct GO taxon constraints catch only **4.4%** ($105/2{,}377$) of cross-clade anatomical false positives (e.g., *heart development* in yeast or plants, *neural crest migration* in arthropods/nematodes, *flower development* in animals). GSPA's GO–Uberon/PO/CL bridge + Stage 2 eliminates **100.0%** ($2{,}377/2{,}377$), while semantic ancestor relaxation (`relax`) improves ancestor-propagated CAFA $F_{\max}$ over naive deletion (`remove`) on all 6 species (`+0.009` to `+0.029`) and rescues **79.8%** ($459/575$) of borderline true developmental annotations.

*(Full report: [`results/pillar3_report.md`](results/pillar3_report.md))*

### Pillar 4: Stage 3 Neuro-Symbolic Agent (`gspa-agent`) Evaluation

| Configuration (120-Anomaly Benchmark Across 7 Genomes) | Precision | Recall | $F_1$ Score | Overall Accuracy | Ontology Violation Rate | Latency / Anomaly |
|---|---:|---:|---:|---:|---:|---:|
| **Config 1:** Stage 1 + Stage 2 Deterministic Only | 100.0% | 0.0% | 0.0% | 50.0% | **0.0%** (`0/120`) | 0.00 s |
| **Config 2:** `gspa-agent` Fast 1-Step Mode (`--max-turns 2`) + Guardrail | **90.7%** | 81.7% | 86.0% | 86.7% | **0.0%** (`0/120`) | 3.74 s |
| **Config 3:** `gspa-agent` `Qwen 3.8 27B` Multi-Turn (`--max-turns 6`) + Guardrail | 88.5% | **100.0%** | **93.9%** | **94.2%** | **0.0%** (`0/120`) | 5.31 s |
| **Config 4:** `gspa-agent` `Qwen 3.8 27B` Multi-Turn **without** `SymbolicGuardrail` | 78.9% | 98.5% | 87.8% | 86.7% | **6.7%** (`8/120`) | 5.29 s |

*(Full report: [`results/pillar4_agent_report.md`](results/pillar4_agent_report.md))*

---

## Repository Structure

- **`gspa`** — Stand-alone CLI launcher script (`./gspa integrate|evaluate|annotate|visualize|compare|report|agent`).
- **`gspa-core/`** — Core Groovy/Java library: GO/Uberon/PO/CL/NCBITaxon ontology reasoning (OWL API + ELK + SAT4J), Stage 1 `IterativeRefiner` and 6 genome-scale priors, Stage 2 `ConsistencyEnforcer` (`relax` / `remove`), `ReactionLocalContextSuggester` (RLGC S2), Stage 3 `AnomalyExtractor` & `AgentCompletionRunner`, and bundled constraint tables (`src/main/resources/`).
- **`gspa-predictors/`** — Tool wrappers (`DIAMOND`, `MMseqs2`, `Pfam`, `InterProScan`, `FoldSeek`, `gapseq`, `DeepGO-PlusPlus`, `ProteInfer`, `CLEAN`, `AMRFinderPlus`, `antiSMASH`, etc.).
- **`gspa-cli/`** — Picocli command-line interface and ShadowJar packaging (`gspa-2.0.0.jar`).
- **`gspa-nf/`** — Nextflow workflow for high-throughput cluster/cloud execution.
- **`benchmark/agent/`** — Stage 3 Neuro-Symbolic Agent sidecar (`run_gspa_agent.py`, also bundled inside `gspa-2.0.0.jar`).
- **`benchmark/unified/`** — Reproducible benchmark suites for Pillars 1–4.
- **`scripts/`** — Programmatic OBO Foundry OWL constraint extraction scripts (`extract_anatomical_constraints.py`, `extract_complex_terms.py`).
- **`results/`** — Empirical TSV summary tables and markdown analysis reports across all 4 benchmark pillars.
- **`paper/`** — LaTeX manuscript (`main.tex`, `main.pdf`), bibliography (`references.bib`), and figure generator (`generate_figures.py`).
- **`docs/`** — Stand-alone CLI and schema documentation (`docs/CLI_REFERENCE.md`).

---

## License

GSPA is licensed under the **GNU General Public License v3.0 or later (`GPL-3.0-or-later`)**. See [`LICENSE`](LICENSE) for full details.
