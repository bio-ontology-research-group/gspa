# 6th place — stefanistefanov — 4-model NN ensemble

Code: https://github.com/stefanistefanov/CAFA6-Solution . A clean **all-neural**
multi-input fusion ensemble (no boosting, no kNN) — useful as the "single NN that
fuses every modality" reference, and a third independent vote for **NEA + PubMed text**.

## Inputs (per protein)
- **ESM2-3B** (`esm2_t36_3B_UR50D`) sequence embeddings.
- **ProtT5 / ProstT5** sequence embeddings.
- **ProstT5 of 3Di** sequences (3Di generated from AlphaFoldDB structures via foldseek)
  — a structure channel without needing a folding GPU per query.
- **PubMedBERT** (`NeuML/pubmedbert-base-embeddings`) over title+abstract+chemicals+MeSH
  of the **5 most recent** PubMed articles (PMIDs from UniProt ID-mapping); **mean ⊕ max**
  of their embeddings.
- **Non-experimental GO terms + evidence codes** (IBA/IBD/IEA/IGC/IKR/IRD/ISA/ISM/ISO/
  ISS/NAS/ND/RCA) as inputs (`test_terms_evidence.tsv`).
- **Taxon lineage** (`ete3.NCBITaxa`), and **normalized sequence length** (len/10000, clipped 1).

## Architecture
Heterogeneous fusion:
- **TermEvidenceEncoder** — 2-layer Conv1D over the term-evidence one-hot matrix → term
  logits (idea from **its7171's CAFA5 3rd-place** NonExperimental module — the common
  ancestor of all the NEA approaches).
- **TermEvidenceEmbeddingEncoder** — 1-layer transformer over (term ⊕ evidence)
  embedding sequences → mean-pooled 512-dim; term/evidence embeddings learned.
- **TaxonEmbeddingEncoder** — 1-layer transformer over the taxon lineage → 512-dim.
- Each PLM/structure/text embedding → FF layer → 1024-dim; all concatenated with the two
  encoder outputs + seq-length.
- **Per-aspect heads** (BP/CC/MF), FF+GELU+Dropout(0.5)+linear; hidden 3072/1024/1536.
- A **logits gate** per head controls how much the TermEvidenceEncoder logits influence
  the final prediction (gated by the fused feature vector).

## Training
4 models differing in input mix and GOA release (226 May-2025 / 229 Dec-2025) and
protein set (SwissProt vs +TrEMBL). 15 epochs, AdamW (wd 0.01), OneCycleLR, batch 128,
**BCEWithLogitsLoss with label weights = sqrt(IA) capped at 0.5**. Last-epoch checkpoint.

## Lessons for DG++
- **Third independent confirmation that NEA (term-evidence) + PubMed text are core
  signals** — and that the term-evidence encoder (its7171) is the canonical way to
  consume electronic annotations in a NN.
- **ProstT5-of-3Di** is a neat way to inject structure as an embedding (no per-query fold)
  — complements our existing ProstT5 head and the foldseek component.
- Loss weight `sqrt(IA)` capped — yet another point on the "how to weight by IA" spectrum
  (4th: frequency, not IA; 5th: 1+IA; 6th: sqrt(IA) capped 0.5). Treat IA-weighting of the
  loss strictly as a tunable ablation, never a default assumption.
