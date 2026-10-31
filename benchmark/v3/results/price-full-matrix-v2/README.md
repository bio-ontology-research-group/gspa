# Complete Price export collection audit

On 2026-10-27 all 132 parts in the frozen v2 export plan had verified publication
receipts, totaling 1,526,919,081 nonzero complete/preliminary EC scores. Counts
match every planned matrix. The final two array tasks completed with exit 0.

This archive records the independent post-publication collection audit, job 10109.
`verify_price_export_collection.py` binds every matrix comparison and publication
receipt to the plan, checks complete predictor/assembly coverage and rehashes
every published artifact and original input (851 files). Its seven integrity
regressions pass. The source and job script are frozen in `launch-inputs.json`.

The audit completed and `verification.json` reports all 851 files verified,
all 132 parts, and 1,526,919,081 claims in 3,359.51 seconds. The complete progress
log is archived. Slurm had already expired job 10109 when terminal state was
requested, so no terminal Slurm record is claimed; the successful verifier receipt
and final log record establish this audit's result.
Remote output: `/data/hohndor/gspa-v3/price-full-collection-audit-v1/verification.json`.
Original data and frozen plan remain in `/data/hohndor/gspa-v3/price-full-matrix-v2/`.
The independent per-part matrix comparisons already ran during export; this
collection audit verifies their coverage, provenance and unchanged published
bytes. Neither step establishes biological annotation accuracy.
