# Full Price reference evaluation launch

Job 10432 is pending after successful collection audit 10410. This is a launch
checkpoint, not an accuracy report. Six predictors, 22 genomes, 132 cells,
1,526,919,081 acquired claims; 148 eligible proteins and 151 EC references per
predictor. Raw exported scores, no-prior integration and default-prior integration
are compared. No native-call, learned/context or full-system claim is made.

All 15 evaluator tests pass. Sonnet reviewed prerequisites.py and run.sbatch
without blocking findings. Its deployment caveats are checked: all 11 package
hashes match remotely, slurm/ exists, Python is 3.12.3, and prior jobs establish
writable node004 scratch. The scheduler receipt confirms afterok:10410,
KillOInInvalidDependent=Yes and no requeue.

The metadata preflight checks 73 completed cells at its recorded checkpoint.
It verifies small receipts, matched commands and file presence, not large-file
contents or biological accuracy. It deliberately does not score partial panels.
The dependent full evaluator will rehash its full inputs and output artifacts;
check_retained.py then independently checks reduced-table arithmetic.

The five source files exactly match the already tested price-reference-v1
package. Production inference and frozen collection files are unchanged.

Remote package: /data/hohndor/gspa-v3/price-reference-full-v1
Expected outputs: jobs/10432/{evaluation,independent-arithmetic.json,exit-code.txt}

Verify this portable launch archive with:

```bash
sha256sum -c archive.sha256
```

All twelve acceptance requirements remain open; no release has been made.
