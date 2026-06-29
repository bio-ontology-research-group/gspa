# 1st place — GOAlpha (ZhuLab, Fudan) — private wFmax **0.5243**

Team: Hancheng Liu, Tianyang Huang, Huiying Yan (ZhuLab, Prof. Shanfeng Zhu).
GOAlpha = the lab's **GOLabeler / NetGO / NetGO-3.0 Learning-to-Rank** lineage
re-run on CAFA6. All references are their own prior work.

## Method
Seven heterogeneous component scores → **XGBoost Learning-to-Rank (LTR)** → final.

| component | category | source |
|---|---|---|
| BLAST-KNN | sequence | weighted kNN vote over BLAST homologs |
| SVM-ESM2 | sequence | per-GO SVM on ESM2 embeddings |
| LR-InterPro | domain | LogReg on 118,677-dim binary InterPro vector |
| FoldSeek-KNN | structure | GO transfer from FoldSeek structural neighbours |
| Net-KNN | network | guilt-by-association over STRING PPI |
| GOXML | literature | extreme multi-label over protein–literature units |
| GORetriever | literature | retrieve GO candidates from similar descriptions/lit, **cross-encoder rerank** |

Plus **GO-term frequency** and a **21-dim species encoding**, all fed to the LTR.

## The decisive move — a test-matched, temporally-clean validation set
They built a validation benchmark that **mimics the species distribution of the
CAFA6 test set**: sampled 1,000 experimentally-annotated proteins to the same
species mix, excluded from training, used *only* for model selection. And: "all
evidence must be restricted to information available **before** the test protein
annotations were obtained" — domain DBs, PPI, GOA, and literature included, not
just sequence.

## Figures (read directly)

**Fig 2 — performance on the species-matched validation benchmark (weighted F-max):**
Blast-KNN ≈ **0.59**, LR-InterPro ≈ 0.58, SVM-ESM2 ≈ 0.575, FoldSeek-KNN ≈ 0.56,
GOXML ≈ 0.55, GORetriever ≈ 0.54, **Net-KNN ≈ 0.45 (clearly the weakest single
component)**, **GOAlpha ensemble ≈ 0.68** (≈ +0.09 over best single).

**Fig 3 — public vs private leaderboard, per component:**
- Blast-KNN 0.29→**0.38**, FoldSeek 0.30→**0.385**, SVM-ESM2 0.305→**0.405**,
  LR-InterPro 0.335→**0.385** — sequence/structure/PLM all **rise** on private.
- **Net-KNN 0.29→0.26 — the ONLY component that DROPS on private.**
- GORetriever 0.425→0.42 (flat); **GOXML 0.48→0.42 (best on *public*, drops on
  private)** — literature tops public but does not generalize.
- GOAlpha 0.49→**0.524**.

## Lessons that transfer to DG++
1. **Integration, not new models, is the lever.** Best single private component ≈
   0.40 (SVM-ESM2); LTR ensemble 0.524 = **+0.12 from integration**. Our whole 2025
   submission (0.377) ≈ one mediocre single component — the +0.12 ensemble lift is
   what we left on the table.
2. **Public LB lies.** GOXML wins public, drops private. Chasing a single LB-like
   number selects the wrong model. Hence LTR over a generalization-faithful val set.
3. **Net-KNN is the weakest and the only anti-generalizing component** — independent
   confirmation of our own clean-eval finding that `diam` dominates `net`. Even a
   correctly-built Net-KNN tops out ~0.45 val / 0.26 private. Temper net expectations.
4. **A combination of homology search + PLM "recovers a substantial fraction of the
   final performance"; structure and literature add complementary signal for hard
   proteins** (their A7, verbatim intent).
