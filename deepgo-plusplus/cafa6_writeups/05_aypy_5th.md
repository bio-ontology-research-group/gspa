# 5th place — AyPy — private wFmax **0.43468** (gold)

The solution **most structurally identical to DG++-Light**: 3 supervised models + 3
kNN label-transfer modalities (sequence / structure / PPI), merged by a NaN-aware
weighted average with shared postprocessing. Read the pipeline diagram (fig1) — it is
essentially a DG++ blueprint with three signals we lack: **NEA**, **MedCPT text**, and
**taxon constraints in postprocessing**.

## Final ensemble = NaN-aware weighted nanmean of 6 sources
| source | category | weight | aspects |
|---|---|---|---|
| Set 1 (PLM + taxon-Phylum + NEA) | supervised | 0.20 | BP/MF/CC |
| Set 2 (PLM + taxon-Order + NEA) | supervised | 0.16 | BP/MF/CC |
| Set 3 (PLM + taxon-Phylum + **MedCPT PubMed** + NEA) | supervised | 0.16 | BP/MF/CC |
| Foldseek-kNN | label-transfer | 0.16 | BP/MF/CC |
| DIAMOND-kNN | label-transfer | 0.16 | BP/MF/CC |
| **PPI-kNN (STRING)** | label-transfer | 0.16 | **BP only** |

Each "Set" is itself `MCM-NN ×0.4 + LogReg ×0.3 + PyBoost ×0.3`. NaN-aware merge:
where a source has no prediction (e.g. PPI on MF/CC), its weight is dropped and the rest
renormalize.

## Features (the 4 families)
- **PLM**: ESM2-3B (2560) ⊕ ProtT5-XL (1024), mean-pooled = 3584-dim.
- **Taxon hierarchy one-hot**: raw leaf taxon is too sparse in CAFA6 → climb to
  Phylum/Order level; rare buckets → `[UNK]`. (Sets differ by taxon level.)
- **PubMed MedCPT** (Set 3 only): per-protein PubMed title+abstract → MedCPT-Article
  encoder (768-dim, HF `ncbi/MedCPT-Article-Encoder`).
- **NEA (Non-Experimental Annotation)** — GOA `(protein, GO, evidence_code)` with
  non-experimental codes (IBA/IEA/ISO/ISS/ND…); after propagation, **top-1000 most
  frequent non-experimental GO terms** as features. For NN: kept as a **structured
  tensor preserving the evidence-code axis** → dedicated NonExperimental Module. For
  LR/PyBoost: mean over evidence axis → flat vector.

## MCM-NN (= SPROF-GO Max-Constraint Module — same as our esm2_head)
Two-branch: Main MLP (flat features) + NonExperimental Module (NEA tensor); logits summed
→ sigmoid → **MCM**. MCM = parent score ≥ max over descendant scores, via a precomputed
ancestor matrix `M ∈ {0,1}^{K×K}`. **Asymmetric**: at inference max over *all*
descendants; at training, if parent y=1 max over *positive* descendants only (prevents
the "high negative child satisfies parent loss" shortcut). Loss = BCE on post-MCM scores
(MCLoss) — folds hierarchy into the objective, removes max-propagation post-step.
Training: AdamW (wd 0.02), cosine 1e-3→1e-6, EarlyStop on EMA val loss, **per-label loss
weight `1+IA`**, mixup (p0.8 α0.2), EMA 0.99, label smoothing 1e-3. **IA=0 terms kept**
for MCM-NN (auxiliary task → richer hierarchy backbone), dropped for LR/PyBoost. BP needs
tiny batch (8 vs 256) because the MCM `[B,K,K]` tensor scales as K² (BP K=4747).

## kNN scoring formula (DIAMOND-kNN & PPI-kNN) — note vs our `net`
    weight       = bitscore (DIAMOND) | combined_score/1000 (PPI)
    score(term)  = Σ_neighbour (weight · IC(term)) / Σ_neighbour weight
    IC(term)     = −log2(p(term)) / max_IC      (training-set frequency)
**Two differences from our `net`/`build_net_component.py`:** (a) they multiply each
vote by the term **IC**, and (b) they normalize by **Σ weight** (a weighted mean), not by
the per-protein **max**. Foldseek-kNN uses ProFun's similar-in-spirit scheme.

## PPI-kNN is BP-only — external validation of our net finding
> "Although BP terms correlate strongly with the functions of neighbouring proteins,
> little gain was observed on MF and CC in my experiments. So PPI-kNN is used only for BP."

## CV — leak-free, group-aware
`MultilabelStratifiedGroupKFold` (5 folds): **identical sequences = same group** (even
across different accessions) so duplicates never straddle train/val; rare classes get
positives in every fold. (They note CD-HIT/MMseqs2 clustering would be even safer.)

## Postprocessing (shared by every component)
GO max-propagation → scorable-GO filter → **drop pairs already in `train_terms`** →
**taxon constraint (forbidden-GO drop)** → top-N (200) per (protein, aspect) → score≥0.01.

## Lessons for DG++ (high value)
1. **NEA** is the biggest missing cheap feature, and it slots straight into our integrator.
2. **PPI/net → restrict to BP** (matches GOAlpha Fig 3 and our clean-eval).
3. **kNN IC-weighting + Σweight normalization** may improve our `diam`/`net`/`foldseek`.
4. **taxon-constraint + drop-known-labels postprocessing** = our SAT differentiator,
   here in cheap form; we can do it *better* with gspa-core SAT4J.
5. Group-by-identical-sequence CV is stricter than our GroupKFold-by-protein.
