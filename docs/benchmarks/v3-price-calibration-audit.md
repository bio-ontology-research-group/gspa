# Price calibration metrics and the zero-score baseline

The historical Price reports do not establish useful enzyme probabilities merely
because the learned arms reduce Brier score or ECE relative to stage0. An
unfitted baseline that assigns zero to every candidate has a lower Brier score
than **13 of the 14 learned predictor/model combinations**. This baseline offers
no useful ranking: every candidate is tied. Its purpose is to show how heavily
these metrics depend on the very large number of unlabeled candidate pairs.

This is an arithmetic audit of all seven stored reports and all six arms, not
a rerun of current GSPA or a new biological benchmark. Each report evaluates
148 labelled proteins and 151 reference-positive protein/EC pairs in 22 genomes.
The denominator includes other emitted protein/EC pairs under the historical
closed-world convention. Those other pairs are **unconfirmed**, not experimentally
proven negatives. The pair universe differs substantially between predictors;
absolute Brier/ECE values across those universes are not matched comparisons.

| Predictor evidence | Evaluated pairs | Zero Brier | Historical logistic Brier | Historical tree Brier |
|---|---:|---:|---:|---:|
| CLEAN | 774,954 | 0.000195 | 0.001656 | 0.002749 |
| DeepProZyme | 401,840 | 0.000376 | 0.001798 | 0.002902 |
| EnzBERT | 14,852 | 0.010167 | 0.010387 | 0.011793 |
| GraphEC | 750,508 | 0.000201 | 0.001794 | 0.002741 |
| MAPred | 3,997 | 0.037778 | 0.040820 | 0.037724 |
| TopEC | 294 | 0.513605 | 0.603806 | 0.513833 |
| Four-predictor pool | 784,148 | 0.000193 | 0.001663 | 0.002744 |

For a binary reference label and a predicted probability of zero, squared error
is one only on a positive pair. Brier therefore equals `151 / evaluated_pairs`.
All predictions occupy the zero-score bin, so equal-width-bin ECE has the same
value. Neither quantity requires fitting a parameter or choosing a threshold
from the test results. The MAPred tree is the sole learned combination with a
lower Brier score than this baseline; that observation alone does not establish
accurate ranking or biological calibration.

The historical learned arms can reduce Brier/ECE while losing the correct EC
ranking. For example, CLEAN top-1 falls from 60/148 to 2/148 for logistic and
0/148 for the tree. The [complete Price report](v3-price-development.md) retains
these outcomes. Current predictor identity and input-range gates were introduced
after these historical runs; their efficacy must be evaluated separately.

Future comparisons should retain these full-universe metrics alongside raw-score
and production baselines, fixed-calibration-threshold precision/recall, top-k
ranking, candidate coverage, corrections and harms, and uncertainty across
genomes. Experimentally characterized alternatives and null cases are still
needed to validate biological probability claims. The zero-score baseline must
not be turned into evidence that a protein lacks a function.

## Reproduction and evidence

The [complete audit](../../benchmark/v3/results/price-calibration-baseline-audit-v1/report.json)
records all 42 arm summaries, report hashes, pair counts, ECE, top-1 results and
differences from the zero baseline. The script cross-checks each universe count
against every arm's reliability bins. It does not certify the original raw inputs
or inference from those summaries.

```bash
python3 benchmark/v3/results/price-calibration-baseline-audit-v1/audit.py \
  --reports benchmark/v3/results/price-context-ablation \
  --out "$TMPDIR/price-calibration-baseline-new.json"
```

The command refuses to overwrite an existing report. No annotations, model
parameters, thresholds or original benchmark artifacts were changed.
