# Context comparison audit v1

Complete descriptive audit of the finished context-mode run, without refitting.
Read [the result and limitations](../../../../docs/benchmarks/v3-context-comparison-results.md).

Run the portable check from any directory:

```sh
python3 benchmark/v3/results/context-comparison-audit-v1/verify_archive.py
```

`analysis/` retains the original 6,706 scored rows, thresholds and all primary
summaries/contrasts. `job-10223/audit/` contains the new 6,981-file hash ledger and
expanded EC report with every reference locus for all 36 arm/policy combinations.
`references/` contains the exact labels, splits, registry and parent registration;
the large GO ontology remains remote and hash-bound. `frozen-models/` holds all
12 fitted models and six table manifests. `metadata/` holds production/parity
receipts and all twelve per-genome production audit summaries. `original-job/`
and `logs/` retain execution status, wrappers and original logs. `frozen-source/`
contains the collected frozen Python sources; `retained-source-bindings.json`
explicitly identifies the 71 retained files checked against original bindings.
`audit-source/` preserves the new audit and five focused tests. Sonnet's read-only
review and its adjudication are retained; the review preceded expanded EC results.
The original supplement's registry wording is preserved as historical evidence;
the independent audit corrects that interpretation using the parent registration.

The portable verifier checks every retained file, the original dependency-chain
receipt, reference membership, all 18 thresholds, 36 original metric summaries,
EC/GO disaggregation and 36 expanded EC panels from the retained per-locus records.
It does not import either the original evaluator or the new cluster auditor.
The cluster auditor separately reconstructed scores from the full production
TSVs and all 48 paired statistics. Its terminal receipt records 337 consumed
bound inputs rechecked after computation.

Raw production TSVs, complete feature matrices, JVM runtime/JAR, GO ontology and
the 831 MB immutable package are not copied into this archive. Hashes alone are
not a replacement for those files. Full remote replay uses the immutable package
and original export base under
`/data/hohndor/gspa-v3/observed-gap-context-modes-v1-relocation/`, the resumed delta
under `/data/hohndor/gspa-v3/context-resume-v1/jobs/10176/`, and the supplement
under `/data/hohndor/gspa-v3/context-resume-v1/`. Audit job 10223 and its wrapper
are retained under `/data/hohndor/gspa-v3/context-comparison-audit-v1/`.
The per-locus records retain calls, reference matches, top five and top-score
ties, not complete score arrays. Ranking reconstruction requires the bound
remote TSVs. This archive proves numerical/report consistency, not independent
biological accuracy, causal context benefit, or superiority over CLEAN.
