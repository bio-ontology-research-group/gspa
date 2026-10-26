# v3 community benchmark: diagnosis of the all-negative first run

Run under diagnosis: `community-development-0501c14` (1,035 Wintermute & Silver 2010 pairs x 6 arms, iML1515).
Every primary prediction was negative (TP 0, FN 178, TN 430, FP 0 per arm); in the trace-present arms six
strains (guaB, leuB, ppc, pyrD, serA, serC) grew alone and the only growing communities were the 15 pairs
of those six; in the trace-absent arms nothing grew. The completed reports are not modified and no unknown
was reinterpreted.

## Verdict

**A formulation bug, not a biological transport limitation.** The community builder kept each model's native
exchange bounds on the ports. iML1515 (BiGG convention) closes uptake (lower bound 0) on 307 of its 331
exchanges; only the 24 members of its own minimal medium can be imported. Every amino acid, nucleobase and
vitamin exchange therefore had `lb = 0`: an organism could secrete the metabolite into the shared pool but
**no organism could take it back up**. Cross-feeding was impossible by construction, so every rescue was
negative regardless of the biology, and the "growth only when both already grow" pattern follows directly.
The pool balance, flux signs and the matched monocultures were correct.

## Evidence (all bounded, no label used to choose anything)

Scripts for this check were scratch scripts; the evidence is reproduced by the regression tests below.

1. **Bounds.** `R_EX_trp__L_e`, `his__L_e`, `leu__L_e`, `arg__L_e`, `ser__L_e`, `ura_e` all have bounds
   `[0, 1000]`; `glc__D_e` is `[-10, 1000]`, `o2_e` `[-1000, 1000]`; 24 of 331 exchanges have `lb < 0`.
2. **Independent COBRApy 0.32.1 / GLPK, single strains.** Wild type grows (0.877 /h at the native medium).
   Gene-knockout (GPR evaluated by COBRApy) plus a 10 mmol/gDW/h supplement rescues the expected auxotrophs:
   trpA+Trp 1.25, hisB+His 0.89, argB+Arg 1.65, lysA+Lys 0.91, pyrF+uracil 0.89, metA+Met 0.90, proC+Pro 1.64,
   icd+Glu 1.50, glnA+Gln 1.54; trpA, hisB, argB, lysA, pyrF, metA, proC, icd, glnA do not grow without it. The
   model, gene (locus b-number to GPR) and amino-acid exchange IDs are as expected.
3. **Independent shared-pool pair (own construction, not `build_lp`).** Two prefixed copies of iML1515 joined
   by a pool metabolite per exchange (each exchange reaction gets +0.5 pool coefficient), supply and removal
   reactions, and `mu` coupled to both biomass reactions. With native bounds every pair gives mu = 0; with
   uptake opened at the ceiling:

   | pair | independent COBRApy | production engine after fix | transfer-forbidden control |
   |---|---|---|---|
   | trpA + argB | 0.8213 | 0.8213 | 0 |
   | hisB + lysA | 0.8212 | 0.8212 | 0 |
   | trpA + hisB | 0.8218 | 0.8218 | 0 |
   | pyrF + metA | 0.8202 | 0.8202 | 0 |

   Agreement to four digits between two independently written formulations; the rescue disappears when
   uptake of non-medium metabolites is forbidden, so it is attributable to the transfer. These pairs were
   chosen as complementary auxotrophs in unrelated pathways (mechanistic), not from labels; their measured
   labels were read afterwards and are *not* the point of this check (trpA+argB synergy, the other three
   intermediate). One viable-monoculture pair (serA+trpA) is correctly negative (`one_viable_monoculture`).

No literature pair is claimed here: the per-pair measured labels are already in the case file, and a
literature comparison against the paired model requires the primary-source figures that were not re-acquired.

## The six strains that grow alone are a separate, real model-versus-biology gap

Not a code bug: COBRApy reproduces them (KO growth 0.859 to 0.877 /h without supplement). Flux inspection of the
bypass (pFBA, no tuning):

| strain | bypass in iML1515 |
|---|---|
| serA, serC | serine from glycine by reverse `GHMT2r` (b2551); the deleted 3-phosphoglycerate route is not required |
| guaB | `XPPT` (xanthine phosphoribosyltransferase, b0238) supplies XMP, so `IMPD` is bypassable |
| pyrD | `DHORDfum`, a fumarate-dependent dihydroorotate oxidation with **no gene assignment**, so no deletion can block it |
| ppc | oxaloacetate regenerated from malate (`MDH`), the anaplerotic reaction is not essential in silico |
| leuB | isozyme `IPMD` b1800 (`b0073 or b1800`); the deletion blocks no reaction |

These strains are measured auxotrophs, so their pairs are expected model misses. They stay in every
denominator. Any change here (e.g. assigning a gene to `DHORDfum`) would be tuning against known phenotypes and is not made.

## Fix

* `scripts/community_metabolism.py`: an optional per-port `open_uptake` (finite, > 0). When present, `build_lp`
  opens the exchange's import direction to that rate, **before** closed ports and uptake closures are applied,
  so closed ports and the control still win. The pool balance, not the bound, limits what can be imported.
  Default (key absent) keeps the previous behaviour, so all earlier manifests are unchanged.
* `benchmark/v3/run_community_benchmark.py`: every port sets `open_uptake` to the protocol ceiling; PROTOCOL
  gets `port_uptake` and `protocol_revision: 2`; the cache key includes it. New protocol SHA-256
  `d62dc7ec...77dc` (pinned in the test).
* Tests: `TestNativeClosedUptake` (community module: native bounds make the transfer impossible, open uptake
  allows it, no nutrient from nothing, uptake closure overrides, invalid values rejected) and
  `test_native_closed_organic_uptake_does_not_prevent_cross_feeding` (runner, BiGG-style closed organics). The
  earlier fixtures gave every exchange `lb = -1000`, which is why the suite missed this.

Risk: this raises the number of open import routes for **every** exchange, including metabolites the model's
authors left closed, so false positives (rescues through a transfer that has no real transporter) are now
possible; they were previously impossible. The control and the per-pair `control_comparison` record whether
a rescue needs the transfer, but nothing checks that a transporter exists. Treat the new run as the first
interpretable result of this model, not as a tuned one.

## Why nothing else needs fixing, and what remains biological

* **Steady state versus 4-day fold change.** The model is static balanced growth at equal fractions; the data are
  day-4 fold growth in sealed microplates. A positive means "mutual growth is feasible if the transfers exist"
  and says nothing about lag, density-dependent release, ratio drift or the thresholds (synergy > 50-fold). The
  intermediate class (41%) has no counterpart in the binary prediction.
* **Export limits.** Free secretion (`removal` and unbounded `ub = 1000` on every export) makes any product
  exportable and the model has no leakage cost, so an optimistic rate is expected; real amino-acid release is slow
  and the paper reports growth only after days.
* **Transport.** After the fix, transport is as constrained as the model's internal transporters (periplasmic
  and inner-membrane reactions keep their bounds); the exchange boundary is no longer the limiter.
* **Model source.** The authors used iAF1260; iML1515 is its successor with different gene content (e.g. the
  GPR of the six strains above). The iAF1260 file was **not** acquired in this session (network access needs
  approval); do not infer that the bypasses above exist in iAF1260.

## Frozen follow-up design (to run by the root, not here)

1. Rerun all six arms with the corrected protocol on the same cases file and model; report as `protocol_revision 2`
   next to, not over, the first run. No condition (oxygen, trace, glucose) is changed after seeing results.
2. Acquire the original iAF1260 (BiGG / BioModels, hash recorded) and a cases file whose model hash and
   gene mapping match; the runner takes any model whose hash equals the case hash. Compare the six single-strain
   bypasses and the pair confusion side by side before choosing anything.
3. Hold-out rule: this set (1,035 pairs) is the development set for model comparison. Any later change motivated by
   a failure pattern (e.g. exclude ports without transporter evidence, per-strain growth gaps) must be evaluated
   on pairs not used to motivate it, declared before it is run.
4. Primary-source check for iML1515 versus iAF1260 (Monk 2017; Feist 2007) of the exchange bounds and gene
   assignments, before any claim about transportability.
