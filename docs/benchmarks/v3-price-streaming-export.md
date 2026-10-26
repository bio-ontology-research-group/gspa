# Price full-matrix streaming export

`benchmark/v3/prepare_price_context.py` writes the Price inference input from the
132 full-proteome matrices (6 predictors x 22 assemblies, about 1.5 billion nonzero
complete-EC scores: CLEAN 601M, DeepProZyme 322M, GraphEC 580M) without holding
claims in memory. The older `prepare_price.py export --raw-matrix-map` path keeps
every claims dictionary in RAM and cannot run on these inputs. It is unchanged.
No production core code is touched, and no new feature contract or scoring formula
is introduced.

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
    --layout <layout dir> --out targets.json
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

* Nothing was run on the full matrices here; tests use synthetic fixtures. Throughput on 45M-claim
  matrices is untested (claims are formatted by string template, not `json.dumps`).
* Sequence identity between the NCBI FASTA and the predictors' original inputs remains unverified.
* `merge` re-reads every claims file once to hash and validate it, which is I/O-bound at this scale.
* A bare predictor id (DeepProZyme) whose versioned form is absent from the FASTA stays unresolved
  and will not match the versioned target id used by labels; the target is then reported absent
  and the part unverified.
* Evidence types are the provisional ones in `prepare_price.PREDICTORS`.
