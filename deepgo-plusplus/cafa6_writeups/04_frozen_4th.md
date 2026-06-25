# 4th place — FrOzen777 (Wenbo Dai) — private wFmax **0.4388**

Individual entry. A multi-method ensemble whose **awarded submission was a single
standalone Py-Boost (0.335 public)**, not the tuned ensemble (0.404 public) — the
ensemble **overfit the public LB**. The clearest cautionary tale in the field.

## Component methods (with public-LB scores)
- **FoldSeek-KNN** (CAFA5-1st style) — 0.233; long seqs with no AFDB entry were dropped.
- **Sprof-GO** (inference only) — 0.248.
- **Py-Boost** (CAFA5-2nd codebase, adapted) — 0.295–0.335. ESM2-650M → **ESM3-3B**
  swap; training set expanded to >150k per CAFA6 rules. The `go-basic.obo` update
  introduced duplicate edges / self-loops needing code fixes.
- **TRGO** — self-developed, DeepGO-SE-inspired **hypersphere matching** multi-label
  model with stacked attention — 0.308–0.334.
- **GOA** — predictions taken from the competition forum.

## Ensembles tried
max / min / median / mean / learned-weight blend → best **0.404 public (25th)**.

## Key findings (verbatim intent)
1. **Ensembling to public overfits.** Awarded = standalone Py-Boost 0.335; the 0.404
   ensemble lost on private. *Have a generalization-faithful internal val before you
   ensemble* (the 1st-place lesson from the other side).
2. **Data quality > quantity.** Core ~150k experimental-only beat the 550k SwissProt-only
   and 1.88M SwissProt+TrEMBL extended sets.
3. **ESM2-3B was the stablest sequence backbone** across public-LB validation runs;
   ESM2-650M also "decent." Structure features (esm3-sm 3Di, Pfam-domain-fused for long
   chains) performed **slightly worse than pure sequence** here.
4. (Earlier analysis) **IA-weighting the *loss* hurt**; frequency weighting helped — so
   keep IA for the *metric*, not necessarily the training objective. (5th place tempered
   this: per-label weight `1+IA` worked for their MCM-NN — treat as an ablation.)

## Lessons for DG++
Our 2025 max-combine submission is exactly this failure mode (0.359 < mlp-alone 0.447 on
no-knowledge): naive max destroyed MF. The fix is a **learned, validation-selected
integrator** — which we shipped (`cafa-baseline` 0.483).
