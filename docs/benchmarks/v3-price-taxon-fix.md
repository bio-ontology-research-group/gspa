# Taxon consistency for enzyme claims: matched Price run

Status: complete (2026-10-31): biology checks and Price scoring for all six predictors (scoring job 15730518).

## What changed in GSPA

The biology checks of the first release candidate (`docs/benchmarks/v3-biology-checks.md`) found enzyme calls that cannot
occur in bacteria surviving into GSPA's output, for three reasons: consistency enforcement was off, the consistency
checker only looked at GO claims, and 21 of the 22 strain taxa were missing from the bundled taxonomy. Release
candidate 2 (JAR sha256 ade3ec60…, staged with a source snapshot in `/project/k10018/gspa-v3/release-candidate-2/`)
checks enzyme claims through their mapped GO terms, accepts a strain lineage (`--taxon-lineage`, manifest column
`taxon_lineage`), and records "taxon not in hierarchy: constraints not applied" instead of silently checking nothing.
Unmapped enzymes are reported as not checkable. Defaults are unchanged. The full JVM gate passes (1,583 tests,
0 failures, 6 optional skips).

## Run

Shaheen collection `/scratch/project/k10018/gspa-v3/out/price-taxon-rc2`: 132 cells (6 predictors × 22 Price genomes),
arm stage 0 (no priors, as before) and arm "with fix" (default priors, `--enforce-stage2 --enforce-consistency-mode
remove --taxon-lineage`, lineages from NCBI Taxonomy 2026-10-30 in `benchmark/v3/taxonfix/price22_taxon_lineages.tsv`).
All cells verified by the collection finalizer (job 15729830). Code: `benchmark/v3/taxonfix/`,
`benchmark/v3/shaheen/price-taxon-cells.sbatch`. Results: `benchmark/v3/results/price-taxon-rc2/` (scores) and
`benchmark/v3/results/biology-taxon-rc2/` (biology checks).

## Results

### Enzyme calls that cannot occur in bacteria (call threshold 0.5, all 22 genomes)

| Predictor | Stage 0: violating calls | Stage 0: genomes affected | With fix: violating calls | With fix: violating retained predictions |
|---|---:|---:|---:|---:|
| TopEC | 411 | 22 | 0 | 0 |
| MAPred | 44 | 22 | 0 | 0 |
| CLEAN | 6 | 6 | 0 | 115,535 |
| EnzBERT | 2 | 2 | 0 | 0 |
| DeepProZyme | 1 | 1 | 0 | 0 |
| GraphEC | 1 | 1 | 0 | 0 |

Stage 0 reproduces the counts of the first release candidate exactly. The CLEAN rows left after the fix are two
enzyme numbers, 1.14.99.51 and its renumbered successor 1.21.3.10, which CLEAN scores at its floor for every protein;
GSPA's enzyme-to-GO table has no entry for either number, so GSPA reports them as not checkable, while the biology
check resolves 1.14.99.51 to 1.21.3.10 and finds the forbidden GO term. None is above the call threshold.

### Agreement with the Price-148 experimental functions

| Predictor | Top-1, stage 0 → with fix | Top-3 | Top calls changed | Corrected | Harmed |
|---|---|---|---:|---:|---:|
| TopEC | 5 → 5 | 5 → 5 | 0 | 0 | 0 |
| MAPred | 50 → 50 | 95 → 95 | 0 | 0 | 0 |
| EnzBERT | 53 → 53 | 76 → 76 | 0 | 0 | 0 |
| DeepProZyme | 46 → 46 | 90 → 90 | 0 | 0 | 0 |
| GraphEC | 18 → 18 | 37 → 37 | 0 | 0 | 0 |
| CLEAN | 60 → 60 | 87 → 87 | 0 | 0 | 0 |

## What this shows

The fix removes every enzyme call above threshold that GSPA's own taxon rules forbid in bacteria, and it
removes no correct Price call for any of the six predictors (CAFA F-max unchanged for all six, `summary_exact.tsv`). It is a correctness improvement, not an accuracy gain on the
Price enzymes, because none of the removed calls was a Price reference function.

## Limits

* Renumbered enzymes absent from the enzyme-to-GO table stay unchecked; resolving renumbered entries before the lookup
  would close the remaining CLEAN gap.
* Strain lineages must be supplied; without them GSPA now says that constraints were not applied.
* An earlier run of these biology checks reported zero violations in both arms because its code snapshot lacked the
  constraint files; the check now stops with an error in that case, and the reported numbers come from the rerun.
