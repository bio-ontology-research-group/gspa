# GSPA Benchmark Status

Last updated: 2026-04-12 14:30 UTC+3

## What works end-to-end

The full GSPA integration pipeline is operational on 9 genomes across
4 phyla (Proteobacteria, Firmicutes, Actinobacteria, Cyanobacteria,
Euryarchaeota). Every component has been validated:

### Predictors
- **DIAMOND blastp** against a leave-9-out Swiss-Prot reference (556k
  proteins). ~10s per genome.
- **HMMER/Pfam** (hmmsearch against Pfam-A). ~10-15 min per genome.
- Both parsed into unified `claims.jsonl` via `02b_parse_predictors_to_claims.py`.

### Evidence integration (Phase 7)
- **Noisy-OR combiner** with correlation-group collapse (DIAMOND + Pfam
  in the homology group).
- **4 active priors** (all validated to fire and produce real posterior
  changes):

| Prior | What fires | Typical activity |
|---|---|---|
| EssentialityPrior | 8-17 uncovered essentials per genome (exact-match, no ELK) | Boosts claims for essential functions not yet covered |
| CoherencePrior | 1,305-2,042 pathway-missing terms per genome | Boosts claims that would close triggered KEGG pathways |
| GapFillingPrior | 135-271 functions per genome (where gapseq data available) | Boosts claims matching gapseq-identified missing reactions |
| GenomicContextPrior | 55-802 claims per genome | Boosts weak claims in operons with pathway consensus |

- **ConsistencyPrior**: architecture validated (SAT4J, UNSAT core
  extraction, 13k taxon constraints from OBO file), but gated on
  `--taxonomy` flag because without per-genome taxon lineage it
  over-penalizes. Not included in the current F-max numbers.
- **Iterative refinement**: converges in 1 iteration without
  ConsistencyPrior (priors are additive, no inter-prior conflict);
  6 iterations when ConsistencyPrior is active.
- **50+ claims cross the 0.5 posterior threshold** per genome due to
  prior boosts (verified on E. coli: likelihood 0.3-0.5 → posterior
  0.5-0.6 after coherence/gap/context boosts).

### Dark matter suggester (Phase 8)
- Fully operational on all genomes with gapseq gap data + operons.
- Namespace bridging fix: MetaCyc gap pathway IDs → KEGG pathway DB
  via GO-term reverse index.
- Synechocystis fix: NCBI `gene_refseq_uniprotkb_collab.gz` mapping
  file recovers 2,913 RefSeq→UniProt mappings where UniProt API
  returned zero.

### Quality metrics (GAEF)
- `gspa evaluate` runs with full ELK (via reasoner cache) or lite mode.
- Completeness, process/pathway/complex coherence, consistency, IC,
  composite score all computed for GSPA and PGAP per genome.

## Current numbers

### F-max (full-GOA truth, 200-bootstrap 95% CI)

| Genome | GSPA | GSPA +priors | PGAP | Ratio |
|---|---|---|---|---|
| ecoli | 0.670 | 0.670 | — | — |
| ecolo157 | 0.835 | 0.835 | — | — |
| bsubtilis | 0.673 | 0.674 | — | — |
| mtb | 0.716 | 0.715 | — | — |
| synechocystis | 0.614 | 0.614 | — | — |
| paeruginosa | 0.601 | 0.598 | — | — |
| hpylori | 0.754 | 0.753 | 0.316 | **2.4×** |
| mgenitalium | 0.913 | 0.912 | 0.469 | **1.9×** |
| mjannaschii | 0.641 | 0.639 | 0.285 | **2.2×** |

Priors hold F-max stable (±0.003) while improving coverage (+0.001 to
+0.068) and IC-recall (+0.001 to +0.027).

### Dark matter suggestions

| Genome | Gaps | Singleton | Disjunctive | Total |
|---|---|---|---|---|
| ecoli | 368 | 913 | 2,903 | 3,816 |
| ecolo157 | 386 | 1,095 | 3,516 | 4,611 |
| paeruginosa | 427 | 1,083 | 3,835 | 4,918 |
| mjannaschii | 124 | 193 | 377 | 570 |
| synechocystis | 277 | 314 | 1,015 | 1,329 |
| hpylori | — | — | — | pending |
| bsubtilis | — | — | — | pending |
| mtb | — | — | — | pending |
| mgenitalium | — | — | — | pending |

## In flight

4 gapseq re-runs (hpylori, bsubtilis, mtb, mgenitalium) are running
on unimatrix01 in `/tmp` (local disk, avoiding GlusterFS `sed -i`
corruption). Started ~07:38, expected completion ~15:00-16:00.
Orchestrator is polling and will auto-chain: parse gaps → integrate
with dark matter → update results.

## Bugs fixed during this benchmark session

1. **`Date.format()` → `SimpleDateFormat`** (Groovy 4 compat in
   QualityPipeline)
2. **`ontology.axioms()` → `ontology.getAxioms()`** (OWL API method
   resolution in GoOntology)
3. **`cls.IRI` → `cls.getIRI()`** (Groovy property access in
   GoOntology, 6 occurrences)
4. **`nodeSet.flattened()` → `nodeSet.getFlattened()`** (ELK return
   type in GoReasoner, 4 occurrences)
5. **EssentialityPrior disabled under `--lite`** — added fallback to
   exact-match boosting without GO descendant expansion
6. **CoherencePrior pathway branch killed by goReasoner guard** —
   split guard so pathway coherence works without ELK
7. **ConsistencyPrior over-penalizing without taxonomy** — gated on
   `--taxonomy` flag
8. **SAT4J UNSAT core `int[]` vs `IVecInt`** — `Xplain.minimalExplanation()`
   returns `int[]`, not `IVecInt`
9. **DarkMatterSuggester pathway-ID namespace mismatch** — MetaCyc
   gap IDs vs KEGG pathway IDs; bridged via GO-term reverse index
10. **Synechocystis zero RefSeq→UniProt mappings** — UniProt API
    `xref_refseq` empty; fixed via NCBI collab file
11. **GlusterFS `sed -i` corruption** — gapseq Reactions.tbl
    zero-byte files; re-run from `/tmp`
12. **Conda Java 21 SIGSEGV** — CDS archive corruption; replaced
    with Eclipse Temurin JDK 21

## Commit history (this session)

```
727f66a Fix synechocystis via NCBI RefSeq-UniProt collab file
26be810 Fix DarkMatterSuggester pathway-ID namespace bridging
8ddc080 Gate ConsistencyPrior on --taxonomy
1a59bb8 Fix SAT4J UNSAT core extraction: int[] not IVecInt
4a321d7 Fix all five priors to actually fire; add OBO taxon constraint loader
8020d3d Benchmark orchestration scripts for the full-priors run
1b3eb0e Wire taxon constraints into gspa integrate; KEGG pathway + gapseq parsers
910b16b Ablation study on 9 genomes
7004cb3 9-genome benchmark: GSPA vs PGAP with bootstrap CIs, dual truth, GAEF
e551e8e PGAP head-to-head benchmark scripts
```
