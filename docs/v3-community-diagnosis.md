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

The evidence below was first obtained with scratch scripts. It is now reproduced by a committed, hashed CLI,
`benchmark/v3/check_community_cobra.py` (section *Reproducible independent oracle*), and by the regression tests.

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
  opens the exchange's import direction to at least that rate (a more permissive native bound is kept), **before** closed ports and uptake closures are applied,
  so closed ports and the control still win. The pool balance, not the bound, limits what can be imported.
  Default (key absent) keeps the previous behaviour, so all earlier manifests are unchanged.
* `benchmark/v3/run_community_benchmark.py`: every port sets `open_uptake` to the protocol ceiling; PROTOCOL
  gets `port_uptake` and `protocol_revision: 2`; the cache key includes it. New protocol SHA-256
  `d62dc7ec...77dc` (pinned in the test).
* Tests: `TestNativeClosedUptake` (community module: native bounds make the transfer impossible, open uptake
  allows it, no nutrient from nothing, uptake closure overrides, invalid values rejected) and
  `test_native_closed_organic_uptake_does_not_prevent_cross_feeding` (runner, BiGG-style closed organics). The
  earlier fixtures gave every exchange `lb = -1000`, which is why the suite missed this.

`open_uptake` semantics, stated precisely: the import bound becomes `min(native lower bound, -rate)`. It opens *at
least* `rate` and never caps, so a native bound more permissive than the rate (e.g. -2000 against a ceiling of 1000) is
kept (pinned by `test_open_uptake_opens_at_least_the_rate_and_never_caps_a_wider_native_bound`). iML1515 has no
such exchange (all native import bounds are 0 or >= -1000), so the benchmark is unaffected; a general per-port use
that needs a cap must close or limit the import by other means. Ports are not opened by default: omitting
`open_uptake` keeps the native-bound semantics on purpose (earlier manifests are unchanged), and an application
manifest must register `open_uptake` per receiving port (`docs/scenarios/community-metabolism.md`). Nothing opens
ports on its own.

Risk: this raises the number of open import routes for **every** exchange, including metabolites the model's
authors left closed, so rescues are now possible that have no real transporter behind them; they were previously
impossible. Only the *exchange boundary* is opened. The model's internal transport reactions (outer-membrane
porins and inner-membrane transporters, with their bounds and GPRs) are retained, so opening an exchange does
not bypass membrane transport: a metabolite still has to cross through reactions the model contains. What this does
not establish is that those reactions correspond to real, expressed transporters in the strain or that the
exchange occurs; a model transport step with a broad GPR is as weak a claim as any other model-conditional one.
The control and the per-pair `control_comparison` record whether a rescue needs the transfer; the independent
oracle (below) additionally records, per rescue, which transferred species the growth depends on and the
consumer's own transport reactions for them. The transporter-evidence experiment (curated or predicted transporter
per port, transporter-knockout phenotypes) has not been done and its outcome is unknown. Treat the new run as the
first interpretable result of this model, not as a tuned one.

## Remaining limits, scoped

What this diagnosis covers is the all-negative pattern of the first run, which is explained by the closed import
bounds. It does not show that nothing else is wrong; the limits below are the ones identified, not a complete list.

* **Steady state versus 4-day fold change.** The model is static balanced growth at equal fractions; the data are
  day-4 fold growth in sealed microplates. A positive means "mutual growth is feasible if the transfers exist"
  and says nothing about lag, density-dependent release, ratio drift or the thresholds (synergy > 50-fold). The
  intermediate class (41%) has no counterpart in the binary prediction.
* **Export limits.** Free secretion (`removal` and unbounded `ub = 1000` on every export) makes any product
  exportable and the model has no leakage cost, so an optimistic rate is expected; real amino-acid release is slow
  and the paper reports growth only after days. Not investigated further here.
* **Transport.** Only the exchange boundary is opened; the model's internal transporters (periplasmic and
  inner-membrane reactions) keep their bounds and GPRs, so the exchange is no longer the limiter but membrane
  transport is still the model's own. Whether those transporters exist and are expressed in the strains is
  **unknown**: the transporter-evidence experiment has not been run, and the oracle shows model dependence on
  transport reactions, not transporter existence.
* **Model source.** The authors used iAF1260; iML1515 is its successor with different gene content (e.g. the
  GPR of the six strains above). The iAF1260 file was **not** acquired in this session (network access needs
  approval); do not infer that the bypasses above exist in iAF1260.
* **Not examined.** Condition choices (oxygen, trace metals, glucose rate) were fixed in advance and not varied
  here beyond the six declared arms; no assessment of the other 1,030 pairs was made in this diagnosis.

## Reproducible independent oracle

`benchmark/v3/check_community_cobra.py` re-solves a fixed, declared selection with COBRApy 0.32.1 and GLPK. It builds
the shared-pool formulation itself (prefixed organism copies, one pool metabolite per extracellular exchange with
coefficient `-c * fraction`, supply and removal reactions, `mu` coupled to every biomass reaction) from COBRApy
objects only, preserves the SBML identifiers (`f_replace={}`), and takes gene deletions from COBRApy's own GPR
evaluation. It imports nothing from the production parser, LP builder, GPR evaluator or solver. Measured labels are
never read (a test shows identical output with labels flipped or removed); the selection manifest
(`benchmark/v3/selections/community_oracle_selection.json`) declares the pairs, arms, protocol values and
numeric tolerances before any run, role `development`, `labels_used: false`.

Selection (declared on 2026-10-26, mechanistic, not held out): four complementary-auxotroph pairs trpA+argB,
hisB+lysA, trpA+hisB, pyrF+metA; one viable-monoculture control serA+trpA; single-strain controls (wild type,
deletion, deletion plus 10 mmol/gDW/h of the missing metabolite) for trpA, hisB, argB, lysA, pyrF, metA; all six
oxygen x trace arms. Agreement tolerances: `mu` relative 1e-6 and absolute 1e-7, identical growth class at the
1e-6 threshold. The report records the SHA-256 of the model file, cases file, selection manifest, production report
and per-arm files, oracle script, plus the production runner, engine and community-module hashes and git commit, the
Python/cobra/optlang/GLPK versions and the solver parameters (library defaults, nothing tuned).

Reproduce (about 6 min on one core; no remote job, no full 1,035-pair loop):

```bash
python benchmark/v3/check_community_cobra.py --cases wintermute2010_cases.jsonl \
    --selection benchmark/v3/selections/community_oracle_selection.json --write-subset-cases subset.jsonl
python benchmark/v3/run_community_benchmark.py --cases subset.jsonl --model iML1515.xml --out prod_fixture
python benchmark/v3/check_community_cobra.py --model iML1515.xml --cases subset.jsonl \
    --selection benchmark/v3/selections/community_oracle_selection.json \
    --production-report prod_fixture/report.json --out oracle_report.json
```

A full revision-2 production run can be validated the same way (`--production-report` of that run, or
`--arm-cases cases_<arm>.jsonl ...` for per-arm files, which leaves provenance marked unverified); only the selected
pairs are re-solved. Committed results:
- Bounded fixture (`benchmark/v3/results/community-oracle-revision2-bounded/`, production at `65ec9bd`, protocol revision 2, status `complete`, 5 pairs x 6 arms).
- Full 1,035-pair x 6-arm revision-2 production run (`benchmark/v3/results/development-3dfe6a6/report.json`, Slurm job 9746) and its oracle verification (`benchmark/v3/results/development-3dfe6a6/cobra_oracle_report.json`, 30/30 pair x arm checks agreeing, max $|\Delta \mu| = 1.65 \times 10^{-13}$). See `docs/v3-community-development.md` for the full 6-arm confusion matrices and strain-resampling intervals.

Result on this selection: all 30 pair x arm comparisons (four simulations each: both monocultures, community,
transfer-forbidden control) agree in growth class and `mu`, largest |delta mu| 1.7e-13. In the arms with trace metals
present, the four diagnostic pairs grow (0.1571 to 0.1574 /h at 0 oxygen, 0.2355 to 0.2360 at 2, 0.8202 to 0.8218 at 20 mmol/gDW/h) while both
monocultures do not, the transfer-forbidden control does not, and the same community with the model's native
exchange bounds (no `open_uptake`) does not; serA+trpA has one viable monoculture and is not a rescue. In the
trace-absent arms nothing grows, not even the wild type. In the trace-present arms the single-strain controls
reproduce the earlier statements (wild type grows; each deletion does not grow without its supplement and does with it).

Which part carries the rescue (one optimal solution per LP, one species closed at a time, given the other
transfers stay open): closing the consumer's own reactions that touch the extracellular species (its porins and
transporters, with the model's bounds and GPRs) removes growth for histidine (hisB) in hisB+lysA and hisB+trpA,
lysine (lysA) in hisB+lysA, and methionine (metA) in pyrF+metA. For trpA+argB no single listed species is
individually necessary, i.e. another transfer can substitute for each; the joint dependence on transfer is shown only
by the control. This shows that, in the model, growth depends on the *model's transport reactions* for those
metabolites, not only on the opened exchange bound. It does not show that these transporters exist or are expressed;
the transporter-evidence experiment is still open and its outcome unknown.

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
