# 3rd place — yuanllo153 — private wFmax **0.4464**

An undergraduate's **reproduction of the U900 team's CAFA5 2nd-place** codebase
(`create_helpers.py` + py-boost/LR/NN → GCN), with only small modifications. Local
val: BP 0.398 / MF 0.712 / CC 0.610. Honest and well-explained; the value is in the
crisp explanation of the U900 recipe and the **target-construction** details.

## Feature extraction
- PLMs: **esm2_t33_650M** + **prot_t5_xl_uniref50**, final layer, mean-pooled. (ESM-IF
  did not help.)
- Training set merged with CAFA5 → **145,382 proteins** (CAFA6 ⊂ CAFA5, roughly).

## Target construction (`create_helpers.py`) — "significant impact on the score"
Convert `train_terms.tsv` → protein × GO matrix with values:
- **1** if annotated, **0** if not, **NaN** if the term *and all its ancestors* are 0.
- **`propagate=True`** — up-propagate positives through the DAG (a labelled child sets
  all ancestors to 1). Toggling this "significantly" changes the local score.
- **NaN for unknown negatives** — a missing term whose whole branch is unannotated is
  *unknown*, not negative; NaN is masked in the loss (don't train it as a negative).
- **raw vs conditional** targets: raw learns P(term) (NaN→0); conditional keeps NaN and
  learns P(term | parent), reconciled to full-DAG probabilities by parent propagation.

## Model (U900 architecture)
Base scores from **Py-Boost + Logistic Regression + Neural Network**, used as input
features to a **GCN stacker** (per BP/MF/CC). Py-Boost has 4 variants:
`pb-t5-4500-{cond,raw}`, `pb-t5esm-4500-{cond,raw}`; inputs T5(1024)+taxon(32) or
T5+ESM(1280)+taxon, outputs 4500 GO (BP 3000 / MF 1000 / CC 500).

## Lessons for DG++
- Confirms the **propagate=True + NaN-masked targets** recipe (our plan §0.3; matches
  the C-HMCNN / MCM idea we already implemented in the ESM2-35M head).
- 3rd ≈ 4th ≈ "CAFA5 redux": **execution and target/validation hygiene**, not novel
  architecture, separated the medal band. Reinforces that our gap is protocol +
  integration, not models.
