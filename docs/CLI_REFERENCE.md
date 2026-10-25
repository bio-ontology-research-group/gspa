# GSPA v2.0.0 — Stand-Alone User Guide & CLI Reference

**GSPA (Genome-Scale Protein Annotation)** is a unified, stand-alone three-stage neuro-symbolic software tool for whole-proteome protein function prediction and genome-scale quality evaluation across Bacteria, Archaea, and Eukaryota.

---

**Workflow defaults:** `gspa presets` lists the presets for `annotate` and
`integrate`. Use `--preset NAME --print-config` to inspect settings; see the
[workflow preset guide](scenarios/presets.md).

**Resources and models:** `gspa resources` lists downloadable IDs;
`gspa download ID`, `--group NAME`, or `--all` retrieves them.
See [resource downloads](reference/resources.md).

## Table of Contents

1. [System Requirements & Installation](#1-system-requirements--installation)
2. [Architecture Overview](#2-architecture-overview)
3. [Subcommands Overview](#3-subcommands-overview)
4. [`gspa integrate` — Three-Stage Genome-Scale Evidence Integration](#4-gspa-integrate--three-stage-genome-scale-evidence-integration)
5. [`gspa evaluate` — GAEF Genome-Scale Quality Assessment](#5-gspa-evaluate--gaef-genome-scale-quality-assessment)
6. [`gspa agent` — Stand-Alone Stage 3 Neuro-Symbolic Agent](#6-gspa-agent--stand-alone-stage-3-neuro-symbolic-agent)
7. [`gspa annotate` — End-to-End Proteome Predictor Execution](#7-gspa-annotate--end-to-end-proteome-predictor-execution)
8. [`gspa visualize` — Interactive Self-Contained HTML5 Genome Browser](#8-gspa-visualize--interactive-self-contained-html5-genome-browser)
9. [Input & Output File Schemas](#9-input--output-file-schemas)
10. [Bundled Ontology & Constraint Resources](#10-bundled-ontology--constraint-resources)

---

## 1. System Requirements & Installation

### Requirements
- **Java 21+** (OpenJDK / Temurin 21 or later) for the stand-alone `gspa-2.0.0.jar` engine.
- **Python 3.9+** (standard library only; `uv` used automatically if installed) for `gspa visualize` and the optional Stage 3 agent (`--enable-agent` / `gspa agent`).

### Option A: Download Pre-Built Stand-Alone Release (Recommended)

Download `gspa-2.0.0-standalone.tar.gz` or `gspa-2.0.0.jar` directly from the [GitHub Release v2.0.0](https://github.com/bio-ontology-research-group/gspa/releases/tag/v2.0.0):

```bash
# Download and extract the stand-alone distribution
curl -LO https://github.com/bio-ontology-research-group/gspa/releases/download/v2.0.0/gspa-2.0.0-standalone.tar.gz
tar -xzf gspa-2.0.0-standalone.tar.gz
cd gspa-2.0.0

# Verify installation
./gspa --version
# Output: 2.0.0
```

Or run the self-contained fat JAR directly from any directory:

```bash
curl -LO https://github.com/bio-ontology-research-group/gspa/releases/download/v2.0.0/gspa-2.0.0.jar
java -Xmx8g -jar gspa-2.0.0.jar --version
```

### Option B: Build from Source

```bash
git clone https://github.com/bio-ontology-research-group/gspa.git
cd gspa
./gradlew clean test :gspa-cli:shadowJar
./gspa --version
```

The self-contained executable JAR is written to `gspa-cli/build/libs/gspa-2.0.0.jar`.

---

## 2. Architecture Overview

GSPA v2.0.0 operates in three composable stages:

1. **Stage 1 — Iterative Bayesian Evidence Integration (`IterativeRefiner`)**:
   - Ingests per-protein functional claims (`claims.jsonl`) across sequence similarity (`DIAMOND`, `MMseqs2`), domain profiles (`Pfam`, `InterProScan`), structural similarity (`FoldSeek`, `DeepFRI`), and protein language models (`ESM-2`, `ProtT5`, `DeepGOWeb`, `ProtInfer`, `CLEAN`).
   - Collapses correlated evidence channels within each evidence group and combines independent groups via log-odds Noisy-OR pooling.
   - Iteratively refines log-odds posteriors (up to 6 rounds, Jacobi damping $\lambda = 0.5$) across the entire proteome using **six genome-scale priors**:
     - `essentiality` (`EssentialityPrior`): Boosts borderline candidates covering missing domain-essential core functions (`bacteria`, `archaea`, `eukaryote`).
     - `coherence` (`CoherencePrior`): Coupled pathway/module completion (KEGG pathways, KEGG Modules, MetaCyc) and obligate heteromeric protein complex completion (`complex_terms.tsv`).
     - `consistency` (`ConsistencyPrior`): SAT4J UNSAT-core taxonomic constraint penalty (`only_in_taxon`, `never_in_taxon`).
     - `anatomical` (`AnatomicalPrior`): Cross-ontology GO–Uberon/PO/CL anatomical constraint enforcement and active-structure developmental boost for multicellular eukaryotes.
     - `gap_filling` (`GapFillingPrior` + `ReactionLocalContextSuggester` RLGC S2): Metabolic gap-filling driven by draft metabolic network gaps (`gapseq`), operon pathway hypergeometric Bayes factors, Pfam coarse EC compatibility, and optional ESM-2 reaction-centroid similarity.
     - `genomic_context` (`GenomicContextPrior`): Operon/synteny neighborhood functional coupling.
2. **Stage 2 — Deterministic Symbolic Constraint & Coherence Enforcement (`ConsistencyEnforcer`)**:
   - Guarantees $100\%$ taxonomic and anatomical consistency on the post-Stage-1 annotation set.
   - Supports semantic **ancestor relaxation** (`--enforce-consistency-mode relax`), which replaces a clade-invalid specific GO term (e.g., *neural crest cell migration* in *D. melanogaster* or *camera-type eye development* in *P. patens*) with its most specific taxonomically and anatomically valid GO ancestor while recovering the unpenalized evidence log-odds.
   - Prunes isolated singleton subunit false positives for obligate heteromeric complexes (`--prune-singleton-complexes`).
3. **Stage 3 — Optional Neuro-Symbolic Agentic Annotation Completion (`gspa-agent`, `--enable-agent`)**:
   - Extracts remaining genome-scale functional anomalies (`PATHWAY_GAP`, `INCOMPLETE_COMPLEX`, `BORDERLINE_ESSENTIAL`, `BORDERLINE_CONFLICT` — typically $<2\%$ of the proteome).
   - Routes each anomaly to an OpenAI-compatible local/cluster SLM/LLM endpoint (default: `qwen3.8-27b`) equipped with 5 read-only genomic inspection tools.
   - Verifies every proposed `PROMOTE`, `RELAX`, or `DEMOTE` action through a deterministic `SymbolicGuardrail` (`SatConsistencyChecker` + `AnatomicalConstraints` + GO vocabulary check), guaranteeing zero hallucinated or clade-violating annotations.

---

## 3. Subcommands Overview

```text
Usage: gspa [-hV] [COMMAND]
Genome-Scale Protein Annotation pipeline
  -h, --help      Show this help message and exit.
  -V, --version   Print version information and exit.
Commands:
  annotate   Run external predictors and annotate a genome or proteome
  evaluate   Evaluate GAEF genome-scale quality of annotations
  compare    Compare two or more annotation sets on a genome
  report     Generate summary quality and annotation reports
  integrate  Run Stage 1/2/3 genome-scale evidence integration
  visualize  Emit a single self-contained HTML browser for a GSPA workspace
  agent      Run the Stage 3 Neuro-Symbolic Agentic Completion sidecar directly
```

---

## 4. `gspa integrate` — Three-Stage Genome-Scale Evidence Integration

### Basic Usage (Bundled Constraints Only — Zero External Databases Required)

Because `gspa-2.0.0.jar` bundles `go-taxon-constraints.tsv`, `ncbi-taxon-hierarchy.tsv`, `go-anatomical-constraints.tsv`, `uberon-po-taxon-constraints.tsv`, `uberon-po-hierarchy.tsv`, `complex_terms.tsv`, and domain essentiality profiles (`bacteria`, `archaea`, `eukaryote`), you can run Stage 1 + Stage 2 immediately on any `claims.jsonl`:

```bash
./gspa integrate \
  --claims claims.jsonl \
  --out integrated.tsv \
  --taxon NCBITaxon_562 \
  --essential-profile bacteria \
  --enforce-stage2 \
  --enforce-consistency-mode relax \
  --lite
```

### Full Prokaryotic Usage (with GO OWL, KEGG Pathways/Modules, Operons, Metabolic Gaps & Stage 3 Agent)

```bash
./gspa integrate \
  --claims ecoli_claims.jsonl \
  --out ecoli_integrated.tsv \
  --provenance ecoli_provenance.json \
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
  --dark-matter \
  --suggestions-out ecoli_dark_matter.tsv \
  --enable-agent \
  --agent-model qwen3.8-27b \
  --agent-api-base http://localhost:8000/v1 \
  --lite
```

### Multicellular Eukaryotic Usage (with Anatomical & Developmental Constraints + Ancestor Relaxation)

```bash
./gspa integrate \
  --claims dmel_claims.jsonl \
  --out dmel_integrated.tsv \
  --go-owl go.owl \
  --taxon NCBITaxon_7227 \
  --essential-profile eukaryote \
  --enable-priors essentiality,coherence,consistency,anatomical \
  --enforce-stage2 \
  --enforce-consistency-mode relax \
  --enforcement-actions-out dmel_relaxed_actions.tsv \
  --lite
```

### Key Options for `gspa integrate`

| Option | Description | Default |
|---|---|---|
| `--claims <FILE>` | **(Required)** Input JSONL file of calibrated per-protein evidence claims | — |
| `--out <FILE>` | **(Required)** Output TSV of integrated posterior probabilities and provenance | — |
| `--provenance <FILE>` | Optional detailed JSON provenance file | — |
| `--go-owl <FILE>` | Path to `go.owl` for DAG ancestor propagation, IC calculation, and relaxation | — |
| `--taxon <TAXON>` | Organism NCBI Taxon ID (`NCBITaxon_562`, `562`, `bacteria`, `archaea`, `eukaryote`) | — |
| `--infer-taxon` | Automatically infer organism taxon from high-confidence claims | `false` |
| `--essential-profile <NAME>` | Built-in core essentiality preset (`bacteria`, `archaea`, `eukaryote`) | — |
| `--pathways <FILE>` | Pathway definitions TSV (`pathway_id`, `name`, `ec_numbers`, `go_terms`) | — |
| `--modules <FILES>` | Additional pathway/module TSVs (e.g., `kegg_modules.tsv`), comma-separated | `[]` |
| `--ec2go <FILE>` | EC-to-GO mapping file (`ec2go.txt`) | — |
| `--operons <FILE>` | Operon assignments TSV (`operon_id <TAB> p1,p2,...` or tab-separated protein IDs) | — |
| `--gaps <FILE>` | Draft metabolic gaps JSONL (from `parse_gapseq_gaps.py`) | — |
| `--enable-priors <LIST>` | Comma-separated active Stage 1 priors | `essentiality,coherence,consistency,anatomical,gap_filling,genomic_context` |
| `--enforce-stage2` | Enable Stage 2 Symbolic Constraint & Coherence Enforcement | `false` |
| `--enforce-consistency-mode` | Stage 2 resolution mode: `relax` (ancestor relaxation), `remove`, `downrank`, `flag` | `relax` |
| `--prune-singleton-complexes` | Prune isolated singleton subunit annotations in obligate heteromeric complexes | `true` |
| `--enforcement-actions-out <FILE>` | Write Stage 2 audit log of `REMOVE`, `RELAX`, and `PRUNE_COMPLEX` actions | — |
| `--dark-matter` | Run the metabolic dark-matter / RLGC S2 enzyme suggester | `false` |
| `--suggestions-out <FILE>` | Output TSV for dark-matter reaction candidate rankings | — |
| `--enable-agent` | Enable Stage 3 Neuro-Symbolic Agentic Annotation Completion (`gspa-agent`) | `false` |
| `--anomalies-out <FILE>` | Export post-Stage-2 anomalies to JSONL (with or without `--enable-agent`) | — |
| `--agent-resolutions-out <FILE>` | Output JSONL of Stage 3 agent resolutions and guardrail decisions | — |
| `--agent-model <NAME>` | Stage 3 model identifier (`qwen3.8-27b` or `qwen-small`) | `qwen3.8-27b` |
| `--agent-api-base <URL>` | OpenAI-compatible VLLM/Ollama API base URL | `$GSPA_AGENT_API_BASE` |
| `--lite` | Skip full ELK `has_part` precomputation (fast mode; recommended) | `false` |

---

## 5. `gspa evaluate` — GAEF Genome-Scale Quality Assessment

Computes the Genome Annotation Evaluation Framework (GAEF) metrics on any genome GFF3 + GAF annotation set:
- **Completeness**: Fraction of domain-essential core functions (`bacteria`: 50 GO terms; `archaea`: 44; `eukaryote`: 65) covered by the annotated proteome.
- **Pathway Coherence**: Mean completeness across triggered metabolic pathways and KEGG modules.
- **Complex Coherence**: Fraction of triggered obligate heteromeric protein complexes (`complex_terms.tsv`) with coherent multi-subunit support.
- **Taxonomic & Anatomical Consistency**: Formal SAT4J verification against direct GO taxon constraints and GO–Uberon/PO/CL anatomical structure constraints.
- **Mean Information Content (IC)** & **GAEF Composite Score**.

```bash
# Convert integrated.tsv + proteome.faa into synthetic GFF3 + GAF2.2
uv run python benchmark/make_synth_gff_gaf.py \
  --fasta bsubtilis.faa \
  --tsv integrated.tsv \
  --gff-out synth.gff \
  --gaf-out synth.gaf

# Run GAEF quality evaluation
./gspa evaluate \
  -i synth.gff \
  -a synth.gaf \
  --go-owl go.owl \
  --taxon NCBITaxon_1423 \
  -k bacteria \
  --pathways kegg_pathways.tsv \
  --modules kegg_modules.tsv \
  --ec2go ec2go.txt \
  --lite \
  -o quality_report.json
```

---

## 6. `gspa agent` — Stand-Alone Stage 3 Neuro-Symbolic Agent

You can run the Stage 3 agent either inline during `gspa integrate --enable-agent` or as a decoupled two-step workflow:

```bash
# Step 1: Extract post-Stage-2 anomalies during integration
./gspa integrate \
  --claims claims.jsonl \
  --out stage2.tsv \
  --go-owl go.owl \
  --taxon NCBITaxon_562 \
  --essential-profile bacteria \
  --pathways kegg_pathways.tsv \
  --ec2go ec2go.txt \
  --operons operons.tsv \
  --enforce-stage2 \
  --anomalies-out anomalies.jsonl \
  --lite

# Step 2: Run the Stage 3 Neuro-Symbolic Agent with SymbolicGuardrail
./gspa agent \
  --anomalies anomalies.jsonl \
  --out agent_resolutions.jsonl \
  --model qwen3.8-27b \
  --api-base http://localhost:8000/v1 \
  --max-concurrency 8 \
  --max-turns 6
```

Environment variables supported by `gspa agent`:
- `GSPA_AGENT_API_BASE`: OpenAI-compatible API endpoint (e.g., `http://localhost:8000/v1`)
- `GSPA_AGENT_API_KEY` or `~/.config/borg-llm/api-key`: Bearer token for authenticated endpoints
- `GSPA_AGENT_MODEL`: Override model name sent to the endpoint

---

## 7. `gspa annotate` — End-to-End Proteome Predictor Execution

`gspa annotate` wraps external homology, domain, structure, and neural predictors (`DIAMOND`, `MMseqs2`, `HMMER/Pfam`, `InterProScan`, `FoldSeek`, `eggNOG-mapper`, `gapseq`, `AMRFinderPlus`, `antiSMASH`, `DeepGO-PlusPlus`, `ProteInfer`, `CLEAN`) from a single configuration file or CLI flags:

```bash
./gspa annotate \
  -i genome.fna \
  -o gspa_workspace/ \
  --kingdom bacteria \
  --threads 16
```

---

## 8. `gspa visualize` — Interactive Self-Contained HTML5 Genome Browser

Generates a single zero-dependency HTML5 file combining protein tables, GO DAG views, `igv.js` genome tracks, operon maps, KEGG pathway/module completion diagrams, and GAEF quality diagnostics:

```bash
./gspa visualize \
  --workdir gspa_workspace/ \
  --genome-id bsubtilis \
  --out bsubtilis_browser.html
```

---

## 9. Input & Output File Schemas

### Input: `claims.jsonl` (One JSON object per line)
```json
{
  "protein_id": "P00561",
  "function_id": "GO:0004072",
  "function_type": "GO",
  "go_aspect": "F",
  "source": "diamond",
  "evidence_type": "SEQUENCE_SIMILARITY",
  "raw_score": 0.912,
  "calibrated_prob": 0.895
}
```

### Output: `integrated.tsv` (Tab-separated)
| Column | Type | Description |
|---|---|---|
| `protein_id` | `string` | Protein accession / locus tag |
| `type` | `string` | `GO` or `EC` |
| `function_id` | `string` | GO term (`GO:XXXXXXX`) or EC number (`1.1.1.3`) |
| `go_aspect` | `string` | `F` (Molecular Function), `P` (Biological Process), `C` (Cellular Component) |
| `posterior_prob` | `float` | Final posterior probability after Stage 1/2/3 ($[0, 1]$) |
| `likelihood_logodds` | `float` | Combined Noisy-OR evidence likelihood log-odds |
| `final_logodds` | `float` | Final posterior log-odds after priors, Stage 2 relaxation, and Stage 3 agent |
| `n_supporting` | `int` | Number of supporting per-protein evidence claims |
| `priors_fired` | `string` | Comma-separated audit trail of prior/stage contributions (e.g., `coherence:0.850,stage2_relax_from_GO:0007507:0.000`) |
| `convergence_iter` | `int` | Number of Stage 1 coordinate-ascent iterations |

---

## 10. Bundled Ontology & Constraint Resources

All ontology constraint files are pre-compiled into `gspa-2.0.0.jar` (`gspa-core/src/main/resources/`) and can be regenerated from fresh OBO Foundry OWL files using the scripts in `scripts/`:

| Bundled Resource Path | Records | Source Ontologies & Extraction Script |
|---|---|---|
| `/taxon-constraints/go-taxon-constraints.tsv` | 1,076 constraints | `go.owl` + `go-computed-taxon-constraints.obo` (`RO:0002160` / `RO:0002161`) |
| `/taxon-constraints/ncbi-taxon-hierarchy.tsv` | 11,466 edges | `ncbitaxon.owl` + `go.owl` taxon hierarchy (`scripts/extract_anatomical_constraints.py`) |
| `/anatomical-constraints/go-anatomical-constraints.tsv` | 2,714 axioms (2,673 GO terms $\to$ 1,294 structures) | `go.owl` logical definitions across 10 RO/BFO relations (`scripts/extract_anatomical_constraints.py`) |
| `/anatomical-constraints/uberon-po-taxon-constraints.tsv` | 8,330 constraints | `uberon.owl`, `po.owl`, `cl.owl` (`scripts/extract_anatomical_constraints.py`) |
| `/anatomical-constraints/uberon-po-hierarchy.tsv` | 29,211 edges | `uberon.owl`, `po.owl`, `cl.owl` `is_a` + `part_of` closure (`scripts/extract_anatomical_constraints.py`) |
| `/complex-terms/complex_terms.tsv` | 1,923 complexes (27 with explicit subunit mappings) | `go.owl` + EBI Complex Portal (`scripts/extract_complex_terms.py`) |
| `/essential-functions/{bacteria,archaea,eukaryote}.tsv` | 50 / 44 / 65 core GO terms | Curated domain-universal core biological processes and molecular functions |
