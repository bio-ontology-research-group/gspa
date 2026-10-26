# Price development regression: full target candidates

At engine revision `e2303a9`, the current pathway-prior configuration produces
**no EC top-1 gain** on these inputs. All seven configurations executed on all
22 genomes. This is development reuse of Price data examined in previous audits,
with limited genomic context; it is not a blind v3 test or the completed
all-module benchmark.

## Matched comparisons

Each configuration uses identical claims in stage0 (priors disabled) and stage1
(default priors, with the same EC-to-GO/pathway/module reference files). The
candidate export retains all nonzero scores for the 148 target proteins. One
of the original 149 entries is excluded because its mapped identifier is a
replicon rather than a scored protein. There are 151 complete-EC labels. The
[input protocol](v3-price148.md) gives all identifiers, exclusions and limitations.

| Exported predictor evidence | Stage0 correct / 148 | Stage1 correct / 148 | Accuracy | All-case report |
| --- | ---: | ---: | ---: | --- |
| CLEAN | 60 | 60 | 40.5% | [JSON](../../benchmark/v3/results/development-e2303a9/price-clean-report.json) |
| DeepProZyme | 49 | 49 | 33.1% | [JSON](../../benchmark/v3/results/development-e2303a9/price-dpz-report.json) |
| EnzBERT | 53 | 53 | 35.8% | [JSON](../../benchmark/v3/results/development-e2303a9/price-enzbert-report.json) |
| MAPred | 50 | 50 | 33.8% | [JSON](../../benchmark/v3/results/development-e2303a9/price-mapred-report.json) |
| GraphEC | 18 | 18 | 12.2% | [JSON](../../benchmark/v3/results/development-e2303a9/price-graphec-report.json) |
| TopEC | 5 | 5 | 3.4% | [JSON](../../benchmark/v3/results/development-e2303a9/price-topec-report.json) |
| CLEAN + DeepProZyme + EnzBERT + MAPred | 60 | 60 | 40.5% | [JSON](../../benchmark/v3/results/development-e2303a9/price-ensemble_pool4-report.json) |

Every matched stage1 comparison has zero top-1 corrections and zero harms.
Reports include genome-level bootstrap intervals (2,000 resamples), fixed-score
metrics, descriptive Fmax, calibration diagnostics, per-genome metrics, execution
records, commands and hashes. Fixed-threshold values do not establish usable
confidence: the current default calibration is provisional and not fitted on
independent calibration data. Calibration summaries across different predictor
configurations use different candidate universes and should not be ranked as if
they were directly comparable.

These numbers describe the **acquired predictor outputs**, not a rerun of each
predictor's published benchmark protocol. Original endpoint training overlap,
reference databases, and some architecture/evidence-type declarations still need
verification. The single-predictor raw rankings give the same target top-1
accuracies as stage0, checked separately below.

## Why this comparison is limited

The export includes target proteins and truncated nearby evidence but has no
gene coordinates, complete proteomes, operon calls, gap input, or organism model.
Reaction-local, cohort, FBA and agent stages were not enabled. Much of the current
pathway/essentiality/context code acts on GO candidates; these inputs are EC
candidates. This result therefore exposes the limited connection between the
current pathway priors and direct EC predictions. It does not test the new
reaction/cohort/learned components, which still need production wiring, real
training and held-out evaluation. It also cannot measure enzyme detection:
every labelled target is already a known enzyme.

## Raw ranking and numerical precision

The [raw ranking audit](../../benchmark/v3/results/development-e2303a9/price-raw-rank-audit.json)
compares each single predictor's original scores with stage0 output and verifies
claims, labels and predictions against the recorded run hashes. All six preserve
the measured target top-1 outcomes. Five preserve the complete exported ranking.
CLEAN has 685,560 adjacent raw-distinct pairs that become tied after the current
sigmoid/log-odds transformations; none changes a unique raw top-1. Full-precision
TSV output fixed decimal-rounding loss, but cannot reverse differences already
lost to floating-point arithmetic. Consequently stage0 is not an exact raw-score
ranking baseline for every CLEAN candidate. The independently scored raw
baseline is retained in the audit rather than inferred from the transformed
posterior.

```bash
python benchmark/v3/audit_price_runs.py --inputs /path/to/price-inputs \
  --run-root /path/to/first-attempt/results \
  --run-root /path/to/retry/results --out /path/to/new-raw-audit.json
```

Slurm array 9717 produced DeepProZyme, TopEC and pool4 results. Other tasks failed
or were cancelled before producing results; node002 retries 9728 and 9731
produced EnzBERT, MAPred, GraphEC and CLEAN results in a separate directory.
Verified ranking audit job 9735 completed afterward. Failed attempts were not
counted as successful execution and their directories were retained.
