# Observed-gap context exports and training-resume checkpoint

All eight arms completed in job 10163: 144 genome/arm runs and 3,026,272 candidate
vectors. Training stopped before model fitting because the remote package omitted
the training registry. The original complete exports remain on Unimatrix and are
hash-bound in export.json; their full 20 GB payload is not duplicated here.

This checkpoint retains the completed-export attestation, failure logs, exact
seven-input supervision supplement and before/after-stage validation wrapper.
Six supervision inputs match prior table-manifest hashes. The registry was not
hashed in those manifests; its current bytes match the pre-existing committed
file. The limitation is explicit. Job 10176 is live at this checkpoint and resumes
at training with the original source, runtime, recipes and completed exports.

No model result, biological benefit, promotion or whole-item acceptance is claimed.
See [scope and status](../../../../docs/benchmarks/v3-observed-gap-context-modes.md).
