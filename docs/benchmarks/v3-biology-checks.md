# Biology checks of the release-candidate Price predictions

Status: complete for the release-candidate Price collection (2026-10-30). Code: `benchmark/v3/biology/`
(`extract.py`, `checks.py`, `report.py`, 9 tests), job `benchmark/v3/shaheen/biology-extract.sbatch`, Shaheen job 15716330.
Results: `benchmark/v3/results/biology-rc1/` (`report.md` has every genome × predictor × arm row, `evaluation.json`,
`receipts.json` with the sha256 of every input file, `gatcab_sensitivity.tsv`).

## What is checked

The checks take the real GSPA output for the 22 Price bacteria and the six enzyme predictors (132 cells) and test it
against biology that holds for essentially all free-living bacteria. Each rule and its source are listed in
`benchmark/v3/biology/checks.py`; the rules were written before the outputs were read. None of the 22 genomes is an
obligate symbiont, so none is exempt.

| Check | Rule | Source |
|---|---|---|
| Translation | a charging route for each of the 20 amino acids (aminoacyl-tRNA synthetases; GatCAB/GatDE accepted for Gln and Asn) | Ibba & Söll 2000; Sheppard et al. 2008 |
| Core enzymes | ten near-universal activities (RNA and DNA polymerase, topoisomerases, DNA ligase, triose-phosphate isomerase, phosphoglycerate kinase, enolase, PRPP synthetase, F-type ATP synthase) | Gil et al. 2004; Hutchison et al. 2016 |
| Taxon consistency | no enzyme whose GO terms are forbidden in Bacteria by GSPA's own bundled taxon constraints | GSPA constraint data (Deegan et al. 2010) |
| Truth agreement | top enzyme call of each Price-148 protein against its experimental EC, by metabolic stratum; other proteins are unconfirmed, never wrong | Price et al. 2018 |

"Call" means a score at or above the engine's call threshold (0.5). "Retained" means any retained score; for the dense
predictors (CLEAN, DeepProZyme, GraphEC), which score every enzyme for every protein, the retained level passes by
construction and carries no information. Only the call level is used below.

## Results

Stage 0 (raw predictor evidence) and stage 1 (full integration with default priors) are byte-identical in all 132 cells,
so every result below is the same for both arms. The integration's priors act on GO claims only, these cells hold enzyme
claims only, no gap or operon inputs were given, and symbolic consistency enforcement was off in this configuration. The
checks therefore measure the predictors, not a GSPA change.

### Translation and core enzymes (call level, genomes passing out of 22)

| Predictor | Translation, rule as registered | Translation, GatCAB counted for both Gln and Asn | Core enzymes |
|---|---:|---:|---:|
| CLEAN | 4 | 19 | 21 |
| DeepProZyme | 9 | 19 | 21 |
| EnzBERT | 4 | 19 | 21 |
| GraphEC | 0 | 0 | 9 |
| MAPred | 4 | 20 | 21 |
| TopEC | 0 | 0 | 2 |

Most failures of the registered rule are asparagine. The GatCAB amidotransferase performs both the Glu-tRNA(Gln) and the
Asp-tRNA(Asn) transamidation (Sheppard et al. 2008), but predictors label it with the glutamine activity (6.3.5.7) only, so
a genome that charges Asn through GatCAB fails the registered rule. The second column counts 6.3.5.7 for asparagine as
well; it is a sensitivity reading added after the results were seen, reported next to, not instead of, the registered
rule. Remaining failures name specific amino acids (Tyr, Arg, Glu) in individual genomes; per-genome detail is in
`report.md`. GraphEC and TopEC miss charging and core enzymes broadly at the call level.

### Taxon consistency (call level)

| Predictor | Genomes with any violating call | Violating call rows (all genomes) |
|---|---:|---:|
| CLEAN | 6 | 6 |
| DeepProZyme | 1 | 1 |
| EnzBERT | 2 | 2 |
| GraphEC | 1 | 1 |
| MAPred | 22 | 44 |
| TopEC | 22 | 411 |

These calls survive into the integrated output because enforcement was off and because the consistency checker only
inspects GO claims; 21 of the 22 strain taxa are also absent from the bundled taxon hierarchy, so the checker could not
have placed them under Bacteria. All three are configuration and coverage gaps of the current release candidate.

### Agreement with the Price-148 experimental functions (top call, four-level canonical EC)

| Predictor | Core biomass (91) | Condition-specific catabolism (55) | Peripheral (2) | All (148), retained top | All (148), at call threshold |
|---|---:|---:|---:|---:|---:|
| CLEAN | 24 | 10 | 0 | 60 | 34 |
| DeepProZyme | 7 | 7 | 0 | 49 | 14 |
| EnzBERT | 30 | 11 | 0 | 53 | 41 |
| GraphEC | 2 | 1 | 0 | 18 | 3 |
| MAPred | 32 | 12 | 1 | 50 | 45 |
| TopEC | 1 | 1 | 0 | 5 | 2 |

The stratum columns count calls at the threshold. CLEAN's 60 of 148 by retained top rank matches the earlier development
run (README, `docs/benchmarks/v3-price-development.md`), an independent consistency check of this collection.

## What this means for the acceptance record

* The suite is the requested biology test of real predictions, and it runs on every collection with one job.
* It shows no GSPA effect on the enzyme predictions in this configuration (identical arms), and it exposes three
  shortcomings to fix or state: consistency enforcement off by default here, EC claims ignored by the checker, and
  strain taxa missing from the bundled hierarchy.
* The Stage 2B collection (`price-stage2b-rc1`, flux-certified refinement switched on) is the arm that can change enzyme
  calls; the same checks will be run on it.
