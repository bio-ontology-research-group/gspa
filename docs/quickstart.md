# Quickstart

This page runs GSPA Stages 1 and 2A on a tiny hand-written set of per-protein predictions. It needs only Java and the GSPA jar; all taxonomic, anatomical, complex, and essentiality resources are bundled.

## 1. Per-protein evidence

GSPA does not predict functions from sequence by itself. It integrates *claims*: calibrated per-protein predictions from any predictor, one JSON object per line. Create `claims.jsonl`:

```json
{"protein_id": "P1", "function_id": "GO:0004072", "function_type": "GO", "go_aspect": "F", "source": "diamond", "evidence_type": "SEQUENCE_SIMILARITY", "raw_score": 0.91, "calibrated_prob": 0.89}
{"protein_id": "P1", "function_id": "GO:0004072", "function_type": "GO", "go_aspect": "F", "source": "interproscan", "evidence_type": "SEQUENCE_DOMAIN", "raw_score": 0.80, "calibrated_prob": 0.75}
{"protein_id": "P2", "function_id": "GO:0007507", "function_type": "GO", "go_aspect": "P", "source": "esm2", "evidence_type": "SEQUENCE_DEEPLEARNING", "raw_score": 0.70, "calibrated_prob": 0.66}
{"protein_id": "P3", "function_id": "2.7.2.4", "function_type": "EC", "source": "clean", "evidence_type": "SEQUENCE_DEEPLEARNING", "raw_score": 0.60, "calibrated_prob": 0.55}
```

Protein P1 has aspartate kinase activity supported by two sources. Protein P2 carries a language-model prediction of *heart development*, which cannot occur in a bacterium. In practice you produce claims with [`gspa annotate`](reference/cli-annotate.md) or the parsers described in the [usage scenarios](scenarios/index.md); the format, including the accepted `evidence_type` values, is specified in [input and output formats](reference/file-formats.md).

## 2. Integrate

```bash
./gspa integrate \
  --claims claims.jsonl \
  --out integrated.tsv \
  --taxon NCBITaxon_562 \
  --essential-profile bacteria \
  --enforce-stage2 \
  --enforce-consistency-mode relax \
  --enforcement-actions-out actions.tsv \
  --lite
```

`--taxon` asserts that the genome is *Escherichia coli*, which activates taxonomic and anatomical consistency checking. `--essential-profile bacteria` loads the bacterial essential core used by the essentiality prior. `--enforce-stage2` runs symbolic enforcement after the iterative refinement, and `--lite` skips the full ontology reasoner, which is the recommended mode for most runs.

## 3. Read the output

`integrated.tsv` lists every surviving annotation with its posterior:

```text
protein_id  type  function_id  go_aspect  posterior_prob  likelihood_logodds  final_logodds  n_supporting  priors_fired  convergence_iter
P3          EC    2.7.2.4                 0.357500        -0.5862             -0.5862        1                            3
P1          GO    GO:0004072   F          0.623000        0.5023              0.5023         2                            3
```

`actions.tsv` records what Stage 2A did and why:

```text
dimension    action  protein_id  term        reason
consistency  remove  P2          GO:0007507  anatomical-inconsistent for NCBITaxon_562: GO:0007507 (heart development)
                                             [results in development of UBERON:0007100 (primary circulatory organ)]
                                             violates only_in_taxon(NCBITaxon_33213) ... (no valid non-root ancestor)
```

The heart-development claim was already down-weighted during Stage 1 by the consistency prior and was then removed in Stage 2A. With `--go-owl go.owl`, `relax` mode replaces such a term by its most specific ancestor that is valid in the organism instead of removing it; here no GO ontology was loaded, so no ancestor was available. See [symbolic enforcement](concepts/stage2a-enforcement.md).

## Next steps

- Annotate a real genome from FASTA: [bacterial isolate genome](scenarios/bacterial-isolate.md).
- Add operons, pathways and metabolic gaps: [metabolic gaps and dark matter](scenarios/metabolic-gaps-dark-matter.md).
- Every option of `gspa integrate`: [reference](reference/cli-integrate.md).
