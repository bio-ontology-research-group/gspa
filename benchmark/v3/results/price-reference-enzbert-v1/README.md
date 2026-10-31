# Complete EnzBERT Price reference panel

Job 10313 completed with exit 0 after all 22 original-production cells finished.
It reads all 11,553,500 acquired whole-proteome claims and compares raw exported
scores, no-prior integration and default-prior integration for all 148 eligible
reference proteins / 151 original EC assertions. All three arms match 53/148 at
top 1 and contain 96/151 reference functions among their candidates. No top-1
identity changes, gains or losses occur. Exact and unambiguous-current EC results
agree. These are acquired-input results, not native published EnzBERT accuracy,
enzyme detection or a full GSPA comparison.

The frozen evaluator rehashes all 506 bound files, including complete ledgers,
and passes all 15 tests. Runtime: 8:45.86 wall time, peak RSS 34,400 KiB. No native
inference, predictions, fitted model or threshold is changed.

Run from any directory:

```sh
python3 benchmark/v3/results/price-reference-enzbert-v1/verify_archive.py
```

The verifier checks every retained file, source identity, tests, exit status and
independent reduced-table arithmetic on 888 protein/arm/identity-policy records.
It also checks the lack of top-1 identity changes directly from the retained score
tables. It does not independently extract those tables from the remote full files;
that extraction and the 506-file recheck belong to the recorded evaluator run.
The original evaluator's Sonnet review remains in `../price-reference-v1/review/`.

See `docs/benchmarks/v3-price-reference-evaluation.md` for population policy,
identifier handling, uncertainty and limits. All twelve whole v3 requirements
remain open; this completes one additional scoped predictor panel.
