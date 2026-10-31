# Historical Price report calibration audit

`audit.py` checks all seven historical Price reports and all six inference arms.
`report.json` binds their hashes and records comparison with an unfitted all-zero
probability baseline. This is report arithmetic, not a new inference run or
experimental-negative validation. The zero baseline has no meaningful ranking.

See [interpretation and reproduction](../../../../docs/benchmarks/v3-price-calibration-audit.md).
