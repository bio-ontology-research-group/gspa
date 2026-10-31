# Price reference evaluation checkpoint

Job 10306: complete TopEC and MAPred panels from the original frozen production
collection, 44 cells, 2,129,370 whole-proteome claims, 148 eligible reference proteins
and 151 EC functions. Both default-prior comparisons show zero top-1 gains or
losses. Details: `docs/benchmarks/v3-price-reference-evaluation.md`.

Run `python3 benchmark/v3/results/price-reference-v1/verify_archive.py` from the
repository root. The verifier checks this archive, source hashes and independent
arithmetic on all retained reference score tables. Large original artifacts stay
remote; their 792 input/artifact bindings and completed rehash are recorded in
the evaluation registration and report. This portable verifier does not reread
those remote bytes or independently certify the extraction.

Fifteen remote tests pass. `review/` retains Sonnet 5.5's read-only findings and
their dispositions. No inference, model or source prediction changed. This is
development reuse, not full v3 acceptance or a release.
