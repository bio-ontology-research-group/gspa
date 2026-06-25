# 2nd place — py-boost + GCN + articles (btbpanda, Fironov, Chervov) — ~0.45–0.48

Direct evolution of their **CAFA5 2nd-place** solution. Code:
https://github.com/btbpanda/CAFA6-updates . The author's own verdict (comment):
the new updates over CAFA5 are worth **≤ 0.04**, and **"70–80 % or more of that is
from the articles."** Articles (literature) is the single biggest new win.

## Pipeline (Fig "overview")
Features → base models (GBDT/LogReg/NN) → **GCN stacker (per BP/MF/CC)** →
postprocessing, with parsed GO annotations fed into the stacker.

## Features
1. **PLM embeddings** — T5 + esm2-small, **pretrained only** ("to our surprise,
   fine-tuning was not helpful").
2. **Articles TF-IDF** (the new win) — for each protein, PMIDs from UniProt
   `lit_pubmed_id` → titles+abstracts via NCBI Entrez → **5,000-dim TF-IDF**. Used
   *cautiously* (leakage risk): only **one** PyBoost model in the final stack used it
   (TF-IDF ⊕ T5).
3. **Taxon one-hot** — ~30 well-represented taxa + one merged bucket. ≈ **+0.02–0.03**
   on base models (CAFA5 estimate).
4. **Cross-ontology SVD** — for limited/partial-knowledge proteins: one-hot the known
   MF/CC labels, **SVD-512**, feed as features to predict BP (and symmetrically). Small
   but "consistently not worse."

## Base models
- **Py-Boost** — their own GPU multi-output GBDT (NeurIPS); 4.5k outputs
  (BP 3000 / MF 1000 / CC 500). Best on CV and LB.
- **Logistic Regression** — 13.5k outputs (10000/2000/1500); weaker on popular terms,
  better on rare ones.
- **Neural Net** — public-notebook architecture, 13.5k targets, T5+ESM.

## Conditional ("alternative") modeling
Targets ∈ {0, 1, NaN}: a term is NaN if no parent is positive; NaN masked in loss.
Model learns **P(term | ≥1 parent exists)**; reconstruct raw probability in DAG order:

    P(GO:N) = P_cond(GO:N) · (1 − ∏_{K ∈ parents(N)} (1 − P(GO:K)))

Lets you predict terms unseen in training (prior mean). They average **raw + conditional**.

## GCN stacker
Node-classification over GO graph (each protein = a graph, shared adjacency). Node
features = 11 base-model predictions × 4 channels (prior-flag, logit, two propagation
logits) + **38 GOA-evidence binary features** + learnable term embedding (dim 8) =
**90 channels**. Trained per ontology, BCE.

## GOA electronic-evidence features (38 binary) — recurring high-value signal
From the GOA `(protein, GO, evidence_code)` triples (electronic: IEA/IBA/ISO/ISS…):
≈ **30 % of electronic labels later become experimental** (≥3 evidences → >50 %). They
encode: "protein has ≥1/2/3/4/5 evidences for term" (propagated and raw), relation type
(part_of / acts_upstream / is_active_in), evidence-code type. Used 2nd, 5th, 6th places.

## Stacker-truncation trick
The metric only scores (protein, ontology) pairs that are experimentally found →
predicting a term in an absent aspect costs nothing. So estimate **P(term | aspect
exists)**: fit the per-ontology stacker only on proteins that have ≥1 experimental term
in that ontology. Small boost + speedup.

## Consistency postprocessing
Final term score = **average of (term prob, max propagated children prob, min
propagated parents prob)** — cheap hierarchy-consistency smoothing.

## Validation
Base models: plain 5-fold CV (better-correlating schemes hurt LB). Stacker: trained on
100 % data, relying on **strong CAFA5 prior** + **stochastic weight averaging** for
robustness (they admit no local val for the stacker is "bad practice").
