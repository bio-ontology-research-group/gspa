# Price full-matrix streaming export

`benchmark/v3/prepare_price_context.py` writes the Price inference input from the
132 full-proteome matrices (6 predictors x 22 assemblies, about 1.5 billion nonzero
complete-EC scores: CLEAN 601M, DeepProZyme 322M, GraphEC 580M) without holding
claims in memory. The older `prepare_price.py export --raw-matrix-map` path keeps
every claims dictionary in RAM and cannot run on these inputs. It is unchanged.
No production core code is touched, and no new feature contract or scoring formula
is introduced.

The historical Price context comparison did **not** consume these full exports;
it reused compact/truncated claims. Its completion claim is superseded by
[the current acceptance record](../../tasks/v3-acceptance.md). Full export,
production execution and biological evaluation are separate gates.

## What is held in memory

One matrix (its pandas pickle or `.npz` must be deserialised whole) plus the row being
formatted. Rows are streamed straight to plain JSONL, so output size does not count.
The largest matrix, not the cohort, sets the memory need. A test bounds the peak of a
4000 x 1500 stream to a fraction of the dense matrix.

## Workflow

Run on the machine that holds the matrices (the root's remote source directory).

```bash
V=benchmark/v3/prepare_price_context.py
# 1. label-free target file: aliases, targets and compact raw stats (verifies the artifact SHA-256 BEFORE unpickling)
python $V targets --artifact price149_all6_compact.pkl.gz --artifact-sha256 <sha> \
    --layout <layout dir> --statistics-projection <frozen enzymeobsolete.json> --out targets.json
# 2. one job per (predictor, assembly): schedule them independently
python $V plan --verification price-source-verification.json > jobs.jsonl
python $V export-part --predictor clean --assembly GCA_000006965.1 \
    --verification price-source-verification.json --trust-pickle \
    --layout <layout dir> --targets targets.json --out parts/
# 3. only after every selected part is complete
python $V merge --parts parts/ --targets targets.json --manifest suite6_price149_manifest.tsv --out merged/
```

`export-part` options per job: `--context-min-score`, `--context-top-k` (context retention),
`--evidence-type`, `--unplaced-declarations`. Use the same retention for every job of
one merge; `merge` refuses mixed retention.

## Outputs

* `parts/<predictor>/<assembly>/claims.jsonl`: one claim per line in the same serialisation as
  `prepare_price._claim` (tested byte for byte), rows in matrix order, ECs sorted within a protein.
* `join_map.tsv`: every matrix id, its canonical id, join status and layout status.
* `part.json`: written last; hashes, counts before and after retention, target coverage, retention.
* `merged/`: pooled ensemble claims (streaming concatenation, uniqueness of
  source/protein/EC validated), `labels/`, `manifests/v3_price_<variant>.json`,
  `unplaced_proteins.json`, `merge_report.json`. A single-predictor variant references its
  part file and is not copied. A claim line costs roughly 200 bytes, so the full cohort is on the
  order of 300 GB before pooling; plan disk accordingly.

## Frozen remote execution and independent verification

The new `price-full-matrix-v1` plan selects all **132** verified original matrices
on Unimatrix, covering all six predictors and 22 assemblies. It contains
**1,526,919,081** nonzero complete **and preliminary** EC scores (the older
1,517,016,382 figure excluded preliminary ECs). No context score cutoff or top-k
is used. Jobs are ordered by nonzero count for capacity measurement; that order
uses no labels and does not change the selected set.

`export_price_matrix_job.py` runs the existing exporter on node-local scratch,
then invokes `verify_price_part.py`. The verifier independently derives EC-column
selection and identifier joins, parses every emitted JSON record and compares it
to its matrix score and identity. It catches changed scores or missing/extra
claims even when a claims checksum has been recomputed. Weak, preliminary,
unplaced and unresolved candidates remain present; a different accession version
is never substituted. The verifier accepts only complete all-nonzero exports.

The worker checks the frozen source/code/layout/target hashes before and after
execution, records package versions, and compares its count with the original
matrix audit. It copies a verified part to an unpublished destination, verifies
the copy's hashes, and renames it into place only when complete. Failure receipts
and Slurm logs remain visible; completed destinations cannot be overwritten.
The per-part claim verification retains the scratch paths used during execution;
`job-receipt.json` binds the published files by relative name and SHA-256.

The current plan and outputs are under
`/data/hohndor/gspa-v3/price-full-matrix-v2/`. Slurm array **9965** verified the two
previously failing regression cases; **9967** processes the remaining 130 jobs,
initially two at a time on node003 (2 CPUs and 8 GB per job). After checking node
capacity, its throttle was increased to four; source/code/input hashes are unchanged.
Bulk output goes to its node-local
scratch before publication to `/data`. At launch the persistent store had 7.2 TB
free and node003 scratch had 3.4 TB free. This removes the local 255 GB storage
constraint without introducing lossy retention. Inspect live Slurm state and
completed receipts; an old job number or log alone is not proof of liveness or
completion.

`run_price_matrix_capacity.py` runs the built production JAR on a published part,
with matched stage-0/stage-1 inputs, verified execution receipts, resource
measurements and an exact candidate-universe comparison. The JVM currently
materializes claims; successful sparse-input execution must not be extrapolated
into a claim that dense matrices or the full cohort fit in memory.

The capacity script's `stage0` means production integration with
`--enable-priors none`; `stage1` enables the default priors. This is a resource
measurement, and its `stage0` name does not mean the raw-predictor maximum-score
baseline used in the accuracy benchmark.

The complete EnzBERT GCF_000236665.1 matrix now passes both arms: **349,400**
candidates from 3,494 rows survive exactly. Baseline integration took 94.47 seconds
with 7,782,496 KiB peak Java RSS; the priors arm took 403.07 seconds with
6,492,324 KiB. Both prediction TSVs are byte-identical. There is no measured
improvement from priors and no biological accuracy estimate in this test.

The subsequent [essentiality-index build](v3-essentiality-index.md) repeats the
complete EnzBERT priors arm in job 10079: 224.56 s Java wall time, 3,531,416 KiB
(3.37 GiB) peak RSS, with an 8 GiB heap/12 GiB allocation. Independent job 10082
checks all 349,400 candidate identities and rehashes 36 recorded files twice;
the full prediction TSV is byte-identical to the earlier stage1 run. Heap and
concurrent-load differences prevent a controlled application-speedup claim.
Original dense-prior job 10058 completed on 2026-10-27, preserving all
9,790,188 DeepProZyme candidates: Java wall time 2:44:34, maximum RSS
52,331,280 KiB (49.91 GiB). Records are in
`benchmark/v3/results/price-capacity-v4/`. Indexed job 10080 completed with exit 0, retaining all candidates: Java wall
time 1:09:50 and peak RSS 63,324,216 KiB (60.39 GiB). Independent comparison
10108 completed with exit 0: every candidate identity matches and the full
prediction TSV is byte-identical. All 37 bound files were rehashed twice.
These runs do not constitute a controlled speedup measurement.

`verify_candidate_universe.py` now compares externally sorted identity keys
instead of materializing the input set and every output row in Python. GNU sort
uses one thread and a 64 MiB buffer, with all temporary files on the caller's job
scratch. Every key is compared; duplicate, missing, extra and changed identities
fail. This bounds the candidate comparison, not the JVM or `verified_batch.py`'s
separate ledger loading. A fresh independent recheck matches both real EnzBERT
outputs and rehashes 33 input/output/code files. Capacity failures preserve a
machine-readable failed record and any completed arm; process kills still require
the Slurm log and must never be interpreted as success.

Frozen capacity job **10037**, under `/data/hohndor/gspa-v3/price-capacity-v2/`,
measures the 9,790,188-candidate DeepProZyme matrix from the same assembly. It uses
a 64 GiB Java heap in an 80 GiB allocation and requests only the no-priors arm.
It has now completed: the no-priors arm preserves all **9,790,188** candidates,
takes **930.31 seconds**, and peaks at **56,327,604 KiB (53.72 GiB)** Java RSS.
The final exact candidate comparison passes. A separate post-run job rehashes
25 recorded inputs/outputs and confirms that their bytes still match the receipt.
This measures one dense baseline, not the default-prior or joint-cohort path.
The output TSV is about 859 MiB and its full stage ledger about 2.5 GiB; these
remain at the recorded remote paths. Execution used the original frozen
`c469273...` JAR. The later candidate-count repair has its own build and tests;
the capacity run was not changed or relabelled mid-execution. The completed
record is in `benchmark/v3/results/price-capacity-v3/`.

The earlier export snapshot at
2026-10-26 17:44 UTC records 70/132 published parts and 55,334,600 claims with no
recorded failures; those progress counts are receipt-based, not a new bulk-byte
rehash. See `benchmark/v3/results/price-capacity-v2/` for the recheck, commands,
resource records, code snapshot and tests. Eight new focused tests pass; the
broader suite reports 523 tests with 11 skips (eight live Nextflow opt-ins and
three other optional checks). The earlier live Nextflow tests apply to the
unchanged route and production JAR.

Eleven new regressions cover independent verification, input mutation, immutable
publication, corrupted copies and compact-statistics projection. The focused exporter/verifier suite passes
37 tests; the complete benchmark suite passes 515 tests with three optional
skips, including the live Nextflow regressions. No JVM source changed here.

These records establish export integrity, not independent accuracy or exact
sequence identity to the original predictor FASTAs. Current acquired NCBI
proteomes and the original matrix inventories differ for some identifiers;
those discrepancies remain explicit in the join maps and existing source audit.

### Real failure found and repaired

Array 9890's first plan was stopped after real failures exposed an incorrect
target-statistics comparison. The original compact extractor first converts
scores to float32 and maps obsolete EC columns into its vocabulary, merging
collisions by maximum. The full exporter compared those statistics to raw values
restricted to that vocabulary. For example, TopEC's raw `2.4.1.129` call for
WP_012434194.1 was stored as `2.4.99.28` by the compact extractor; excluding the
raw identifier from the compact vocabulary falsely reported a missing score.
The same omission also invalidated EnzBERT comparisons. Float32 underflow can
likewise change nonzero counts for very small raw scores.

The retrieved extraction reference `enzymeobsolete.pkl` has SHA-256
`3214ffd31cdea7f19a84a6a94715a754f58d7f82ee591e118fca907a7934f5b8`, exactly matching
the archived PSB code snapshot. Its vocabulary also matches the compact artifact
exactly. This reference reproduces a historical transformation; it does not
establish current EC nomenclature for biological evaluation.

The corrected target metadata declares that projection and includes a SHA-256 of
each complete compact float32 target vector. Verification now reproduces the
projection, compares the entire vector, and separately compares every emitted
raw claim against the original matrix. **Raw scores and EC identifiers are not
changed, rounded or merged.** A regression reproduces obsolete-ID mapping,
many-to-one maximum reduction and float32 underflow while checking unchanged raw
claims; it also detects a nonmaximal score change that counts/maxima would miss.

The corrected real TopEC GCF_000020125.1 part verifies all 7,163 claims across
7,194 matrix rows. EnzBERT GCF_000236665.1 verifies all 349,400 claims across
3,494 rows. Both previously failed; their old failure receipts remain under v1.
All 132 v2 parts now have verified publication receipts, with exactly
1,526,919,081 claims and no per-job count mismatches against the frozen plan.
The last array tasks completed on 2026-10-27. Independent collection audit 10109
has rehashed all 851 bound input/output files and verified complete plan coverage;
its final receipt reports success in 3,359.51 seconds. The seven collection-verifier regression
tests detect changed bulk files/sources, missing or unexpected parts, duplicate
plan entries, incorrect job/count records and disconnected matrix comparisons.
Launch evidence is in `benchmark/v3/results/price-full-matrix-v2/`.
A future biological comparison must handle obsolete/current EC identity consistently;
passing this export check alone does not solve that evaluation requirement.

A part directory is built under a temporary name and renamed only when complete, so an interrupted job
leaves no part and a finished part is never overwritten. Merge output is likewise immutable. A missing
matrix or a SHA-256 mismatch stops the job; the compact artifact is never used as a fallback.

## Retention

* Target proteins keep every nonzero numeric complete EC and every preliminary EC (`1.1.1.n11`).
* Hierarchical partial columns (`1`, `1.1`, `1.1.1`, `1.1.1.-`) are excluded and counted separately.
* Non-target context keeps all nonzero scores by default. A score threshold and/or top-k is explicit and lossy:
  the part then reports `context.pruned: true` and `full_candidate_coverage: false`, and the runner manifest
  says the context is pruned. Pruned context is never called full candidate coverage.
* No rule reads a label. `export-part` has no truth parameter, and a test shows claim hashes are unchanged
  when the manifest truth is mutated while labels change.
* Target candidate coverage is checked per target against the compact artifact's `raw_target_stats`
  (`n_nonzero_vocab`, `max`). A mismatch or an absent target marks the part unverified and `merge` refuses it.

## Protein id join

Exact version first. A bare id (no `.version`) maps only when the NCBI FASTA has exactly one version
of it. A versioned id absent from the FASTA is never moved to another version. An unresolved id keeps its
original string in the claims; `join_map.tsv` records why (`unresolved_absent`,
`unresolved_version_mismatch`, `unresolved_ambiguous_bare`). Duplicate axis ids and two rows resolving to
one protein stop the job. The audit expects 80 versioned matrix proteins absent from the current NCBI FASTAs
(the root is recovering a genomic mapping separately) and 21 bare missing ids.

Proteins acquired without a coordinate are listed from the layout report's `unavailable_loci`. Ids absent
from the FASTA count as unplaced only if the explicit `--unplaced-declarations` file names them; the rest are
reported in `merge_report.json` under `undeclared_unresolved_claim_ids`, which the cohort stage would reject.
No sequence or coordinate is guessed.

## Limitations

* Complete-matrix exports and the capacity measurements above have executed;
  they do not establish inference throughput on the largest 45M-claim matrices
  or the complete cohort. Dense default-prior capacity remains under measurement.
* Sequence identity between the NCBI FASTA and the predictors' original inputs remains unverified.
* `merge` re-reads every claims file once to hash and validate it, which is I/O-bound at this scale.
* A bare predictor id (DeepProZyme) whose versioned form is absent from the FASTA stays unresolved
  and will not match the versioned target id used by labels; the target is then reported absent
  and the part unverified.
* Evidence types are the provisional ones in `prepare_price.PREDICTORS`.

The subsequent [production collection](v3-price-production-collection.md) registers
all 132 parts for both prior settings and repairs multi-GB ledger verification.
Its initial wave and separate dense capacity diagnostic remain in progress;
export completeness must not be reported as completed inference or EC accuracy.

The [complete EC identity audit](v3-ec-identity.md) now checks the old projection
against a pinned official release. Its 101 split-branch choices and six deleted
outcomes must not be reused as unconditional aliases in new scoring. Raw exports
and the historical statistics projection remain unchanged.
