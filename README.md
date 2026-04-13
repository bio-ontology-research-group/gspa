# GSPA — Genome-Scale Protein Annotation

GSPA is a functional-annotation pipeline for prokaryotic and archaeal
genomes, MAGs, and microbial communities. It combines multiple
**evidence types** (sequence similarity, protein domains, orthology,
genomic context, metabolic-model gap analysis) through Bayesian
integration, then uses ontology-derived **priors** (essentiality,
pathway coherence, taxon consistency, metabolic gaps) to boost weak
evidence and suppress taxon-constraint violations.

It ships with wrappers for ~16 external tools (DIAMOND, MMseqs2,
HMMER/Pfam, FoldSeek, eggNOG-mapper, InterProScan, gapseq, AMRFinder,
dbCAN, antiSMASH, SignalP, DeepTMHMM, …) and three built-in quality
metrics (Completeness, Coherence, Consistency) from the GAEF framework.

## Complete example: annotating a bacterial genome

This walkthrough annotates *B. subtilis* 168 from scratch using the
benchmark workflow. Estimated wall-clock times are for a 16-core
server; most steps parallelize trivially.

### Prerequisites

Download once (total ~17 GB):

```bash
# GO ontology + taxon constraints
curl -L -o go.owl http://purl.obolibrary.org/obo/go.owl                                    # 130 MB
curl -L -o go-taxon-constraints.obo \
  http://current.geneontology.org/ontology/imports/go-computed-taxon-constraints.obo         # 2 MB

# EC → GO mapping
curl -L -o ec2go http://current.geneontology.org/ontology/external2go/ec2go                 # 350 KB

# Swiss-Prot reference (for leave-one-out DIAMOND DB)
curl -L -o uniprot_sprot.fasta.gz \
  https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz  # 93 MB

# RefSeq ↔ UniProt ID mapping (for operon protein-ID bridging)
curl -L -o gene_refseq_uniprotkb_collab.gz \
  https://ftp.ncbi.nlm.nih.gov/refseq/uniprotkb/gene_refseq_uniprotkb_collab.gz           # 1.2 GB

# GOA annotations (ground-truth evaluation only; not needed for annotation itself)
curl -L -o goa_uniprot_all.gaf.gz \
  https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz                # 15 GB

# KEGG pathway DB (built from KEGG REST + ec2go)
python3 benchmark/build_kegg_pathway_tsv.py \
  --ec-pathway <(curl -sL https://rest.kegg.jp/link/pathway/ec) \
  --pathway-names <(curl -sL https://rest.kegg.jp/list/pathway) \
  --ec2go ec2go --out kegg_pathways.tsv                                                     # 1 min
```

### Step 1 — Download the genome (~10 s)

```bash
# UniProt proteome FASTA (normalized to bare accessions)
curl -sL "https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000001570&format=fasta" \
  | awk '/^>/{split($0,a,"|"); print ">"a[2]" "substr($0,2); next}{print}' > bsubtilis.faa

# NCBI genomic FASTA + GFF (for gapseq and operon extraction)
curl -sL ".../GCF_000009045.1_ASM904v1_genomic.fna.gz" | gunzip > bsubtilis.fna
curl -sL ".../GCF_000009045.1_ASM904v1_genomic.gff.gz" | gunzip > bsubtilis.gff
```

### Step 2 — Build leave-one-out DIAMOND DB (~2 min)

```bash
# Exclude this genome's proteins from Swiss-Prot
grep '^>' bsubtilis.faa | awk '{print substr($1,2)}' > exclude.txt
python3 benchmark/filter_fasta_by_exclude.py \
  --fasta <(gunzip -c uniprot_sprot.fasta.gz) --exclude exclude.txt --output reference.fasta
diamond makedb --in reference.fasta --db reference --threads 16
```

### Step 3 — Run predictors (~40 min)

```bash
# DIAMOND blastp — sequence similarity (SEQUENCE_SIMILARITY evidence type)
diamond blastp --db reference --query bsubtilis.faa \
  --out diamond.tsv --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore stitle \
  --evalue 1e-5 --max-target-seqs 50 --query-cover 50 --subject-cover 50 --id 30 \
  --threads 16                                                                     # ~10 s

# HMMER/Pfam — protein domain families (SEQUENCE_DOMAIN evidence type)
hmmsearch --domtblout pfam.domtbl --noali -E 1e-5 --cpu 16 Pfam-A.hmm bsubtilis.faa  # ~15 min

# InterProScan — multi-database domain/family/motif scan with InterPro2GO
#   (SEQUENCE_DOMAIN evidence type; runs Pfam, TIGRFAM, CDD, SUPERFAMILY,
#    Panther, ProSite, HAMAP, FunFam, etc. and maps hits to GO via InterPro)
interproscan.sh -i bsubtilis.faa -o interproscan.tsv -f TSV --goterms --cpu 8 -dp  # ~25 min
```

### Step 4 — Parse predictor outputs → claims.jsonl (~3 min)

```bash
python3 benchmark/02b_parse_predictors_to_claims.py \
  --results-dir ./preds \
  --goa goa_uniprot_all.gaf.gz \
  --pfam2go pfam2go.txt \
  --interproscan interproscan.tsv \
  --test-accs exclude.txt \
  --output claims.jsonl
```

### Step 5 — Extract operons from GFF (~1 s)

```bash
python3 benchmark/make_operons.py bsubtilis.gff operons_refseq.tsv 300

# Remap RefSeq → UniProt using the collab file
python3 benchmark/build_refseq_uniprot_map.py \
  --collab gene_refseq_uniprotkb_collab.gz \
  --genomes bsubtilis:refseq_ids.txt:uniprot_accs.txt \
  --out-dir maps/
# Then remap the operon file using the resulting TSV
```

### Step 6 — Run gapseq (~8-10 h)

```bash
gapseq find -p all bsubtilis.fna
# → bsubtilis-all-Pathways.tbl, bsubtilis-all-Reactions.tbl

python3 benchmark/parse_gapseq_gaps.py \
  --pathways-tbl bsubtilis-all-Pathways.tbl \
  --reactions-tbl bsubtilis-all-Reactions.tbl \
  --ec2go ec2go --tag bsubtilis --out gaps.jsonl
```

### Step 7 — Integrate with priors + dark matter (~30 s)

```bash
java -jar gspa.jar integrate \
  --claims claims.jsonl \
  --out integrated.tsv \
  --go-owl go.owl --lite \
  --essential-profile bacteria \
  --pathways kegg_pathways.tsv --ec2go ec2go \
  --operons operons.tsv \
  --gaps gaps.jsonl \
  --enable-priors essentiality,coherence,gap_filling,genomic_context \
  --dark-matter --suggestions-out suggestions.tsv
```

Output:
- `integrated.tsv` — per-protein × per-GO-term posterior probabilities
  with provenance (which predictors, which priors, convergence iteration)
- `suggestions.tsv` — dark-matter singleton/disjunctive suggestions for
  metabolic gaps

### Step 8 — Evaluate quality (~20 s with reasoner cache, ~12 min first run)

```bash
# Convert integrated TSV to synthetic GFF + GAF for gspa evaluate
python3 benchmark/make_synth_gff_gaf.py \
  --fasta bsubtilis.faa --tsv integrated.tsv \
  --gff-out synth.gff --gaf-out synth.gaf

java -jar gspa.jar evaluate \
  -i synth.gff -a synth.gaf --go-owl go.owl \
  --ec2go ec2go --pathways kegg_pathways.tsv \
  --reasoner-cache ./reasoner-cache \
  -k bacteria -o quality.json
```

### Estimated total time

| Step | Wall clock | Notes |
|---|---|---|
| Download genome | 10 s | UniProt + NCBI FTP |
| Build DIAMOND DB | 2 min | One-time per reference set |
| DIAMOND blastp | 10 s | 16 threads |
| HMMER/Pfam | 15 min | 16 threads, ~4k proteins |
| InterProScan | 25 min | 8 threads, ~4k proteins; runs Pfam + TIGRFAM + CDD + SUPERFAMILY + Panther + ... |
| Parse claims | 3 min | GOA scan dominates; ~10s with pre-built subset |
| Extract operons | 1 s | Intergenic distance from GFF |
| gapseq find | **8-10 h** | MetaCyc pathway detection via tblastn |
| Integrate + dark matter | 30 s | All 4 priors + Phase 8 suggester |
| Quality evaluation | 20 s | With reasoner cache |
| **Total** | **~9-11 h** | **Dominated by gapseq** |

Without gapseq (skip step 6, lose GapFillingPrior and dark-matter
suggestions): **~45 min total**.
Without InterProScan (lose TIGRFAM/CDD/SUPERFAMILY/Panther GO hits):
**~20 min total**.

## Project layout

Multi-module Gradle (Groovy 4 + Java 21):

- **gspa-core** — data model, GO ontology (OWL API + ELK), quality
  metrics (SAT4J for taxon-constraint consistency), I/O, config,
  Phase 7 integration engine (`gspa.integration`), Phase 8 dark-matter
  suggester (`gspa.integration.suggester`).
- **gspa-predictors** — `Predictor` interface, `AbstractToolPredictor`,
  gene callers, all tool wrappers, `AnnotationPipeline` orchestrator,
  community / crossfeeding analyzer.
- **gspa-cli** — picocli CLI: `annotate`, `evaluate`, `compare`,
  `report`, `integrate`.
- **gspa-nf** — Nextflow workflow (Docker + BioContainers).
- **benchmark/** — 9-genome head-to-head vs PGAP. See
  [`benchmark/RESULTS.md`](benchmark/RESULTS.md).

## Evidence integration (Phase 7)

Each predictor emits `Annotation` objects lifted into `EvidenceClaim`
records. The integrator:

1. Groups claims by `(protein, function)`.
2. Collapses by evidence-correlation group (homology, structure,
   context, ml-sequence, …) so DIAMOND + Pfam don't double-count.
3. Combines via log-odds **Noisy-OR** with per-type reliability.
4. Iterative refinement (up to 6 rounds, Jacobi damping 0.5):
   `L_post = L_lik + Σ λ_k · prior_k(protein, function, state)`
5. Full provenance per annotation.

### Priors

| Prior | What it does | Needs |
|---|---|---|
| `EssentialityPrior` | Boost claims for uncovered essential functions (exact-match or descendant-expanded with ELK) | Essential profile |
| `CoherencePrior` | Boost claims that close triggered-but-incomplete pathways, scaled by `(1 − coverage)` | KEGG pathway DB |
| `ConsistencyPrior` | Penalize (−3 log-odds) GO terms in SAT4J UNSAT core against taxon constraints | `--taxonomy` lineage + taxon constraints OBO |
| `GapFillingPrior` | Boost claims matching gapseq-identified missing reactions (0.7× for gapseq guesses) | gapseq gaps JSONL |
| `GenomicContextPrior` | Boost weak claims in operons with pathway consensus; extra weight if claim closes a gap | Operons + pathway DB |

All priors degrade gracefully — each is silent when its inputs are
missing. Validated to fire correctly on 9 benchmark genomes with
50+ claims crossing the 0.5 threshold per genome.

## Dark-matter suggester (Phase 8)

For each gapseq-identified metabolic gap `(pathway P, reaction R,
target function f_R)`:

1. **Bayes factor** `BF(O, P)` for each operon — does this operon
   participate in pathway P? Based on current posteriors of members.
2. **Per-protein log-odds** `L_R(p) = L_lik + L_op + commitment_penalty`
   within passing operons (BF ≥ 10).
3. **Softmax** over operon → `q(p)` per protein.
4. **Singleton** if `q(top) > 0.5`; **disjunctive** credible set
   otherwise (cumulative q ≥ 0.9).

Benchmark results on 5 genomes with gapseq data:

| Genome | Gaps | Singleton | Disjunctive | Total |
|---|---|---|---|---|
| ecoli | 368 | 913 | 2,903 | 3,816 |
| ecolo157 | 386 | 1,095 | 3,516 | 4,611 |
| paeruginosa | 427 | 1,083 | 3,835 | 4,918 |
| mjannaschii | 124 | 193 | 377 | 570 |
| synechocystis | 277 | 314 | 1,015 | 1,329 |

## Quality metrics (GAEF)

- **Completeness** — essential-function coverage per kingdom.
- **Coherence** — process (ELK has_part), pathway, complex coverage.
- **Consistency** — taxon constraints via SAT4J; UNSAT core for
  violation explanation.
- **Information content** — mean IC of annotated terms.
- **Composite score** — weighted combination.

## Benchmark results

9-genome head-to-head against PGAP. GSPA outperforms PGAP 1.9–2.4×
on every genome where PGAP has GO annotations:

| Genome | GSPA F-max | PGAP F-max | Ratio |
|---|---|---|---|
| hpylori | **0.754** | 0.316 | 2.4× |
| mgenitalium | **0.913** | 0.469 | 1.9× |
| mjannaschii | **0.641** | 0.285 | 2.2× |

Priors hold F-max stable (±0.003) while improving coverage by +0.001
to +0.068 and IC-recall by +0.001 to +0.027. Largest gains on
mjannaschii (Archaea, where homology is weakest): coverage
0.208 → 0.276 (+33% relative).

Full tables: [`benchmark/RESULTS.md`](benchmark/RESULTS.md) |
[`benchmark/ABLATION_REPORT.txt`](benchmark/ABLATION_REPORT.txt) |
[`benchmark/STATUS.md`](benchmark/STATUS.md)

## Build & test

```bash
./gradlew build                  # compile + unit tests
./gradlew clean test             # fresh build
./gradlew :gspa-cli:shadowJar    # fat jar
```

Requires Java 21+ (Gradle wrapper included).

## License

See `LICENSE` once added. Uses OWL API (LGPL), ELK (Apache), SAT4J
(LGPL/EPL), picocli (Apache), Jackson (Apache), Spock (Apache).

## Citation

Paper in preparation. If you use GSPA, please cite this repository
and the quality-metrics paper the GAEF framework is drawn from.
