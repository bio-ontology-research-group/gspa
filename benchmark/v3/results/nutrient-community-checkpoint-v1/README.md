# Nutrient coculture execution checkpoint

This is a verified implementation/input checkpoint, not a completed empirical
evaluation. Slurm array 10150 is running the full production and independent
GLPK replay. Submission, checksum and running-state receipts are retained.

The outcome-blind plan contains all 1,320 cocultures (1,200 primary), six
medium/oxygen arms and five LPs per case/arm: 7,920 plans, 39,600 production calls.
Separate genotypes and organism models preserve full GPR effects and native
maintenance. The three community policies cover permissive uptake, native
unsupplied uptake and closure of unsupplied uptake. Bulk endpoint OD does not
prove both members grew or identify exchanged metabolites/transporters.

Seven focused tests pass, including independent synthetic reciprocal-rescue
checks and shard/serial identity. The 645-test benchmark gate passes with 13
optional skips. Sixty deterministic preliminary real-model comparisons agree
with GLPK; no claim is made for the full catalogue yet. Initial resource-warning
logs are retained; the final focused gate is warning-free.

`package/` is the exact immutable remote package: source, models and final
plan-v3 bindings. No observation file is in that package. Development plan
versions have identical plan rows but older implementation hashes; their hash
comparison is recorded separately. Use only the final frozen package to run.

Current periodic production messages report original interleaved plan positions,
not per-shard completion counts; some shards emit no periodic messages. The
checkpoint read actual files through a read-only Slurm step after worker SSH
was unavailable: 243 and 234 written records, still awaiting full verification.

Sources, commands, package versions, preliminary results and assumptions are
retained. Paths in receipts preserve the original environment and need remapping
elsewhere. The full source report documents the protocol and limitations. Large
files use deterministic gzip; `storage-index.json` records uncompressed sizes
and hashes. Every payload was expanded and verified against its original file.
Run `sha256sum -c SHA256SUMS` to check stored files.

No full numerical or biological accuracy result, transporter validation,
gene-attribution result or whole-item acceptance is claimed.
