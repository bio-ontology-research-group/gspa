# Protein identity and dense essentiality-prior execution

The 2026-10-26 investigation repaired an identifier crash and removed repeated
candidate scans without changing the essentiality prior's scoring rules.
These are engineering results; they do not validate the prior biologically or
replace the outstanding learned/null model and held-out evaluation requirements.

## Identifier crash

`ClaimKey.parse` split a key at its first two pipes, whereas `IntegrationState`
uses the last two delimiters to preserve protein identifiers containing pipes.
For `sp|P00956|SYI_ECOLI|GO|GO:0006412`, the old typed parser returned null.
Targeted essentiality rescue then crashed while accessing the pinned score.

The typed parser now uses the same delimiter rule as the integration state.
Three new round-trip cases fail before the repair and pass afterwards, including
an identifier containing the literal segment `GO`. A four-run production CLI
replay verifies that the old pipe-ID case fails, the repaired case succeeds,
ordinary-ID predictions remain unchanged, and both repaired predictions agree
after removing only the protein-ID column. This uses one synthetic GO claim;
the UniProt-shaped ID is a syntax fixture, not biological truth.

Evidence: `benchmark/v3/results/claim-key-v1/`. The first repaired build is
`ee72ddd579f606ded6f86a6a7b3af2c91877b6881ac315b5e86714000ae5b665`;
the replay also passes again on the final indexed build below.

## Dense-prior bottleneck

Slurm job 10058 runs the complete 9,790,188-candidate DeepProZyme matrix using
the original frozen capacity build. A live JVM thread dump, obtained through
`java -m jdk.jcmd/sun.tools.jcmd.JCmd`, samples its inference thread inside
`EssentialityPrior.candidates`, parsing keys through Groovy's list-to-array
conversion. The stack includes exception construction during that conversion.
This is a sampled stack, not a complete CPU-time profile.

Previously, each uncovered GO category scanned the entire candidate map. The
32-category EC-only regression performs 64,000 key visits for 2,000 candidates.
It now visits that universe once. Eligible borderline GO candidates are indexed
by term for the current iteration. Category lookup restores their original
encounter order before ranking, preserving legacy ties; targeted ties retain
their existing explicit key ordering. The index is rebuilt each iteration, so
same-size key replacement and likelihood changes cannot leave stale candidates.
Memory for the index scales with eligible GO candidates, not irrelevant EC rows.

`splitFunctionKey` also constructs a `String[]` directly instead of converting
a Groovy list. Scores, thresholds, specificity penalties, pinned-score rules,
rescue logging, ontology closure and failure handling are otherwise unchanged.

The synthetic replay uses both rescue modes, mixed GO/EC and malformed keys,
pipe-containing protein IDs, overlapping category closures, likelihood fallback,
empty claim counts, changed keys/likelihoods, and repeated iterations. All **240
recorded states** from 80 scenarios match exactly for boosts, pinned scores and
rescue logs; 90 recorded states contain rescues. Dedicated tests additionally
check scan count, ordering and same-size replacement. The two scan-count cases
fail before optimization; the two semantic cases already pass before it.

For 100,000 EC candidates and 32 essential categories, the median of three
post-warm-up `beginIteration` measurements falls from **4.493 s to 0.072 s**
(62.4-fold in this synthetic check). This excludes input construction, JVM startup,
other priors, output writing and validation. It is not an end-to-end speedup or
proof of full-matrix memory capacity.

Final JAR SHA-256:
`86805e4019bbe239b43b07664bd7a5d432a60b2a145907a08ddc807c5a6086cf`.
All 1,010 JVM tests complete without failures/errors (two skips). The benchmark
suite completes 525 tests without failures (three optional skips), including
all eight live Nextflow checks against that exact JAR.

Evidence, sources, the live stack, before/after replay and frozen capacity launch
are in `benchmark/v3/results/essentiality-index-v1/`. The launched capacity
checks use complete published EnzBERT and DeepProZyme matrices with exact
candidate-universe verification. **EnzBERT job 10079 completed successfully**:
all 349,400 candidates retained, Java wall time 224.56 s, peak RSS 3,531,416 KiB
(3.37 GiB), 8 GiB heap and 12 GiB allocation. Total measured time including
verification is 243.40 s. Independent job 10082 rechecks every candidate and
recorded file hash, and the entire prediction TSV is byte-identical to the
earlier complete stage1 run. Both Slurm jobs report COMPLETED and exit 0.

The older run used a different heap/allocation and had different concurrent
load; these resource measurements are not a controlled end-to-end speedup test.
The unchanged DeepProZyme reference job 10058 completed on 2026-10-27 with
all 9,790,188 candidates retained: 2:44:34 Java wall time and 52,331,280 KiB
(49.91 GiB) maximum RSS. Its terminal state, receipt and resource records are in
`benchmark/v3/results/price-capacity-v4/`. Optimized job 10080 also completed with exit 0, retaining all candidates:
Java wall time 1:09:50, peak RSS 63,324,216 KiB (60.39 GiB), and total measured
time including verification 4,550.47 s. The frozen engine used 64 GiB heap and
80 GiB allocation. Its records and terminal state are in the same archive.
Independent comparison 10108 completed with exit 0. It independently compares
all 9,790,188 identities in each arm, verifies byte-identical complete predictions,
and rehashes all 37 recorded files twice. The TSV SHA-256 is
`64701c3dd55cd6dc753139df7b93eb10e49f1efad361ea943f19f42f7aa6136f`.
Different concurrent loads prevent interpreting the wall-time difference as a
controlled speedup. Larger-matrix and cohort capacity remain unproven.

## Reproduction

Use new output paths and set `TMPDIR` and `JAVA_TOOL_OPTIONS` to disk-backed
scratch locations. Run the same Groovy replay through the pre-index and indexed
production JARs, then compare the complete `cases` arrays and inspect timings:

```bash
java -Xmx2g -cp "$BEFORE_JAR" groovy.ui.GroovyMain \
  benchmark/v3/check_essentiality_index.groovy "$BEFORE_JSON"
java -Xmx2g -cp "$AFTER_JAR" groovy.ui.GroovyMain \
  benchmark/v3/check_essentiality_index.groovy "$AFTER_JSON"
python benchmark/v3/check_claim_key.py --before "$CRASHING_JAR" \
  --after "$AFTER_JAR" --out "$IDENTITY_REPLAY"
```
