# Dense Price baseline capacity: DeepProZyme

Slurm job 10037 completed production baseline integration and exact retention
verification for the complete GCF_000236665.1 DeepProZyme matrix: 9,790,188
candidates from 3,494 rows. Java takes 930.31 seconds with 56,327,604 KiB peak
RSS (53.72 GiB). These resources cover Java, excluding Python input/output
verification. They include writing the output artifacts to persistent storage.

The job requests only the production no-priors combiner (`stage0`), not the
raw-predictor maximum-score accuracy baseline. It uses the original immutable
c469273224723bff080efd9440d1c1ae07fc0c943c2cacffa04a98003481f801 JAR;
the subsequent candidate-count fix was not substituted into this live job.
Frozen runner/helper code and launch parameters are in the preceding
`price-capacity-v2/launch/` archive. The output TSV (about 859 MiB) and full stage
ledger (about 2.5 GiB) stay at the remote paths bound in `verified_run.json`.

The capacity script checks every candidate identity using bounded-memory
external sorting. An independent post-run job (10053) rehashes 25 recorded
inputs/outputs after measurement. The record and runnable script are included.
No label, pruning, sampling or top-k substitutes for full candidate retention.
This is neither an accuracy result nor proof of default-prior, larger-matrix or
joint-cohort capacity.

`export-progress.json` is a receipt-only snapshot of the ongoing 132-part export.
Its array throttle was increased from two to four after checking available
node CPU/memory. Frozen source/code/input hashes are unchanged. This snapshot
counts published worker receipts; it does not rehash all bulk export bytes.
See docs/benchmarks/v3-price-streaming-export.md and tasks/v3-acceptance.md.
