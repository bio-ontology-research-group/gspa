# Role-truth registry (`role_truth_cases.jsonl`)

One JSON object per line describes a **reviewed, complete candidate locus**, not a sufficient
null-model training set. Candidate membership is limited by `universe_note`; absence of other
annotations is never evidence for a negative.

## Fields

- `case_id`: unique, stable slug; also used as `locus_group_id` by `build-loci`.
- `status`: `accepted`, `pending_review`, or `rejected`.
- `genome_id`, `taxid`, `organism`: exact genome/taxon identity. Evaluation-only genomes may use a
  benchmark identifier rather than an assembly accession.
- `in_training_cohort`: boolean. `false` means evaluation-only and excludes the case from
  `build-loci`; the evaluator ignores this flag.
- `function_type`: `GO` or `EC`; `function_id`: full primary identifier (`GO:…` or `EC:…`).
- `role`: candidate role, normally `UNSPECIFIED` for GO and `CATALYTIC` for EC.
- `members`: at least two objects, with unique `protein_id` within the case. Each contains
  `protein_id`, `swissprot_accession`, `label`, free-text `role_truth`, and non-empty `basis`.
  `label` is integer 0/1 or `null`. A null label requires `role_truth: "unknown"` and a basis
  explaining why the member is neither a clean positive nor a clean negative. Unknown members never
  count as positives or negatives. Empty `swissprot_accession` is accepted only for a null-label
  member. Every case needs at least one positive and at least one negative. Exactly one positive is
  required unless optional boolean `allow_multi_positive` is true. A positive may optionally declare
  `evidence_code`; `build-loci` uses it when adding a missing declared-positive label row.
- `references`: non-empty list of source identifiers; at least one must start with `PMID:`, `DOI:`,
  or `PMC`.
- `source_statement`: free text locating the experiment/table/figure supporting the labels.
- `reference_files`: list of archived source objects with relative `path` (resolved below
  `--references`) and lowercase SHA-256 `sha256`. Every listed digest is checked during validation
  and locus building. UniProt `<swissprot_accession>.txt` records are optional: when present,
  `build-loci` also checks accession, taxon, versioned protein cross-reference, and full sequence;
  when absent it records a skipped sequence check.
- `universe_note`: precise scope and completeness limit of the candidate set.
- `added`: ISO date; `reviewer`: review provenance (pending cases must still say who supplied the
  initial review).

## Optional GAF cross-check

`gaf_assertion` is used only for cases derived from the GAF negative inventory. It contains
`inventory` (path relative to `--references`, or absolute), `accession`, `relation`, `evidence`,
`reference`, and `candidate_status`. The builder verifies the inventory schema/status and candidate
file digest, exactly one matching assertion, and the expected genome/protein candidate join. The
inventory snapshot date also bounds the UniProt sequence-version check.

`build-loci` uses only `in_training_cohort: true` cases whose status passes `--statuses`; it never
silently omits a selected case whose exported candidate universe is incomplete. Unknown members get
no new label rows, and any case containing one bypasses locus certification with reason
`candidate_universe_incomplete_by_declaration`: its candidate universe is not fully labelled by
declaration. Known positive and negative members are still emitted.

When a label-1 member has no original positive row for the case function, the builder adds the
declared positive with `evidence` from optional `evidence_code` (default `IDA`), `tier` =
`experimental`, the case's primary `reference`, and the member's Swiss-Prot accession. Every such
row is listed in `verification.json` under `declared_positive_additions`. Original rows remain
immutable; conflicting original label-0 truth for a positive, or an original positive for a
negative, is refused.

The evaluator measures within-locus ranking and threshold calls only, not genome-wide performance.
It reports unknown members under `unknown_scores`, excludes them from negatives, ranks, and call
counts, and reports `n_unknown_above_active`. `active_rank` and MRR use strict competition rank
(tied negatives count as above); Top-1 is reported as both `top1_strict` and `top1_lenient`. An
unscored active has `active_rank: null`.

## Negative-only observations (`role_truth_negative_only.jsonl`)

This file uses the same schema except that no positive member is required. These are
**negative-only observations: usable for false-positive counting, never as locus groups**. Validate
them with `validate --negative-only`. The evaluator's `--negative-only FILE` option reports
`negative_only_called_at_threshold / negative_only_total` separately and excludes unknown members.
