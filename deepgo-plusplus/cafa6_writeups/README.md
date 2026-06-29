# CAFA6 winning-solution writeups — distilled

Source: the public Kaggle solution writeups for **CAFA6 Protein Function
Prediction** (saved HTML under `~/Downloads/*Place*Kaggle*.html`, fetched
2026-06-19/25). One markdown per ranked solution plus a cross-cutting
[`SYNTHESIS.md`](SYNTHESIS.md) that turns the recurring patterns into an
actionable backlog for **DeepGO-PlusPlus** (and DG++-Light).

| file | team | private wFmax | one-line |
|---|---|---|---|
| [`01_goalpha_1st.md`](01_goalpha_1st.md)   | GOAlpha (ZhuLab, Fudan) | **0.5243** | 7 components → XGBoost **Learning-to-Rank**; species-matched validation |
| [`02_pyboost_gcn_2nd.md`](02_pyboost_gcn_2nd.md) | btbpanda et al. | ~0.45–0.48 | py-boost + LogReg + NN → **GCN stack**; **PubMed TF-IDF** = the new win; GOA electronic-evidence features |
| [`03_u900_repro_3rd.md`](03_u900_repro_3rd.md)  | yuanllo153 | **0.4464** | reproduction of U900 CAFA5-2nd (py-boost+LR+NN → GCN); propagate=True + NaN targets |
| [`04_frozen_4th.md`](04_frozen_4th.md)      | FrOzen777 | **0.4388** | multi-method max/avg ensemble; **awarded run was standalone py-boost 0.335** (ensemble overfit public) |
| [`05_aypy_5th.md`](05_aypy_5th.md)          | AyPy | gold | 3 supervised (MCM-NN+LR+pyboost) + 3 kNN (Foldseek/DIAMOND/**PPI BP-only**); **NEA** + MedCPT text + taxon constraints |
| [`06_stefanov_6th.md`](06_stefanov_6th.md)  | stefanistefanov | — | 4-model NN ensemble fusing ESM2-3B + ProstT5(3Di) + **PubMedBERT** + **term-evidence encoder** + taxon lineage |

**The three signals every top team used that DG++ does not yet:**
1. **Non-experimental GOA annotations + evidence codes as input features**
   ("NEA"; 2nd, 5th, 6th, and its7171's CAFA5 3rd). ~30 % of electronic
   labels later become experimental → a strong, cheap, CPU-only feature.
2. **Literature / PubMed text** as its own evidence stream (1st GOXML/
   GORetriever, 2nd TF-IDF, 5th MedCPT, 6th PubMedBERT). DG++ has a BM25
   `lit` component — the top teams use learned text embeddings instead.
3. **Taxon encoding** (one-hot ~30 taxa / lineage) into the supervised
   models and the integrator (1st 21-dim, 2nd ~30, 5th hierarchy, 6th lineage).

See [`SYNTHESIS.md`](SYNTHESIS.md) for the full prioritised backlog.
