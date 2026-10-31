# Nutrient phenotype benchmark evidence

This archive contains the complete registered production run, its biological
scores and independent numerical checks. It is not a whole-item acceptance
certificate.

- 1,144 assays, six fixed arms, 6,864 conditions, 13,200 production LPs.
- Predictions: 1,889 growth, 4,974 no growth, one solver unknown.
- Final independent replay: all 13,200 GPR sets agree; 13,198 flux comparisons
  agree. One production unknown and one finite baseline objective discrepancy
  remain unverified; strict exit is 1.
- Original replay: three unverified flux comparisons. Fresh-model diagnosis
  resolves two finite numerical differences caused by GLPK basis history.
- No arm beats the 80% constant-positive supplementation-rescue comparator.
  Known-genotype phenotype feasibility is not independent gene-function truth.
- Twelve focused tests and 638 benchmark tests pass (three optional skips).

## Contents

`data/plan/` contains the frozen protocol and complete outcome-free execution
plans. `data/remote/jobs/10144/` retains the models, original source, production
results, package versions, failed original oracle and exit status.
`data/remote/oracle-v2/10148/` contains the revised oracle source and complete
output. `data/numerical-verification.json` summarizes the final comparison.
`data/evaluation.json` and `data/constant-baselines.json` retain all strata,
unknowns and simple comparators. Large files have a `.gz` storage suffix.

`observations/` binds the raw-derived observations and their source audit/model
registration. The full original source acquisition is in the sibling
`nutrient-environments-v1` archive. `verification/` includes exact commands,
Slurm records, numerical diagnoses, test logs and failed/cancelled attempts.
Its `source-access/` subdirectory contains preliminary acquisition notes for
future datasets, not additional validation labels or a complete source archive.
`source/` snapshots the final benchmark implementation and documentation.
Original and revised oracle sources are intentionally different and separately
hashed; production predictions and their frozen implementation are unchanged.

## Integrity and reproduction

From this directory, run `sha256sum -c SHA256SUMS` to check stored files.
`storage-index.json` maps each logical payload to its stored path and gives its
original byte length and SHA-256. For a `.gz` payload, decompress before comparing
with that original hash. Both stored and decompressed payloads were independently
verified. Integrity verification is separate from scientific acceptance.

Follow `source/docs/benchmarks/v3-nutrient-phenotypes.md` for the production,
oracle and evaluation commands. Use the archived source versions and models
when replaying the frozen plan. A strict oracle exit of 1 is expected for this
registered prediction set because the production unknown and finite discrepancy are retained.
