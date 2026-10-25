# Essentiality rescue: legacy vs targeted

Evaluation of the rewritten essentiality prior (`gspa.integration.prior.EssentialityPrior`, rescue mode
`targeted`) against the v2.1.0 behaviour (`legacy`) on the nine reference genomes of the ablation (P02).

## Why

The AllTheBacteria macro-evolution analysis (`benchmark/macroevolution/results/t10_rescue_floor.tsv`)
showed that the v2.1.0 rescue leaves the rescued essential claims straddling the decision threshold: the
refiner caps a prior lift on a sub-threshold claim at log-odds 0.12 (posterior 0.53), damps each iteration
by 0.5, and Stage D ran three iterations, so a rescued claim ends between 0.47 and 0.53 and only 44 % of
the rescues cross 0.50. Whether an essential function counts as present was therefore close to a coin flip.

## What changed

* **Targeted rescue.** For every uncovered essential category the candidates (Stage 0 likelihood in
  [0.16, 0.50] on a term under the category) are scored: Stage 0 log-odds, corroboration across claims
  (`wSupport * ln n_sup`), a bonus when the protein is otherwise unannotated (`wFree`), a penalty when it
  already carries a confident annotation of the same aspect outside the category (`wCompeting`), and the
  fraction of its operon partners annotated to the essential domain (`wContext`). The top candidate is
  pinned to the posterior `rescueTarget` (default 0.60) through the refiner's floor mechanism, so it is
  present at the threshold by construction.
* **Rescue log.** Every rescue is written next to the output as `<out>.rescues.tsv` (category, claim,
  likelihood, score, target), so rescues can be audited against curated annotations.
* `legacy` mode reproduces v2.1.0 exactly (tier ranking, boost `max(1.30 alpha, 0.14 - L0)`, no floor).
  Theta keys: `ess_rescue_mode`, `ess_rescue_target`, `ess_min_candidate`, `ess_max_rescue`,
  `ess_w_support`, `ess_w_free`, `ess_w_competing`, `ess_w_context`.

## Run

`run_bench9.py` runs on unimatrix01 (`/data/hohndor/gspa/ess_rescue/`): nine genomes x {A2 essentiality
only, A7 all Stage 1 priors, A9 Stage 1 + Stage 2A relax} x {legacy, targeted}, then evaluates rescue
precision against `truth_all` (term level and category level), essential completeness at 0.50, the number
of essential-domain claims in [0.45, 0.55), and micro and CAFA F-max (`benchmark_pgap_v2.py`).
Results in `results/`.
