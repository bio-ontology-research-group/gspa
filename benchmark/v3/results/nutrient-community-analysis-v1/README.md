# Coculture analysis checkpoint

This archive freezes the dependence-sensitive analysis for all 1,320 coculture
cases and 7,920 execution plans. Full simulation/independent checking is still
running in Slurm array 10150; no empirical accuracy or whole-item acceptance is
claimed. Seven focused tests and the 652-test benchmark gate pass (13 optional
skips). The initial mutable-alias failure and successful repair are retained.

Each payload is gzip-compressed with its original and stored SHA-256 in
`index.json`. `SHA256SUMS` binds every stored file except itself. The unchanged
execution package/models/plans are in the preceding
`nutrient-community-checkpoint-v1` archive. The new freeze binds its exact
execution protocol and source hashes. Results of the full run will be archived
separately.
