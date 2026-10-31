# Role-truth benchmark, evaluation archive v1 (2026-10-28)

Registry: `benchmark/v3/role_truth_cases.jsonl` (sha256 81bdd93b…, 8 cases: 5 accepted, 3 pending review)
and `benchmark/v3/role_truth_negative_only.jsonl` (sha256 5ed517c0…, 2 negative-only loci).
Evaluator: `benchmark/v3/evaluate_role_truth.py --format integrated_tsv --score-field score` with
`--negative-only`. Strict top-1 requires the active member strictly above every evidence-labelled
negative; lenient top-1 accepts ties. Unknown members are never ranked or counted.

| folder | scores evaluated | sha256 of the M. tuberculosis table |
|---|---|---|
| `legacy-bench9` | April-2026 legacy integration `/data/hohndor/gspa/proteomes/bench9/final/mtb_integrated.tsv` (9-source claims incl. DIAMOND+Pfam, additive priors) | 49b60586… |
| `job10615-jar10526-priors-none` | legacy bench9 claims re-integrated with JAR `/data/hohndor/gspa-v3/supplied-gff-v1/jobs/10526/gspa-3.0.0.jar`, batch manifest, priors `none` (`/data/hohndor/gspa-v3/role-truth-v1/jobs/10615`) | 0d668670… |
| `job10619-jar10526-layout` | as 10615 plus `--genome-layout` built from the RefSeq GFF (`role-truth-v1/layouts`) | 2e796f69… |
| `simple-baseline-legacy-claims-max` | best single claim score per candidate from `bench9/mtb_claims.jsonl` (evaluator `--format claims_jsonl`), no integration | — |
| `job10622-jar10526-operons` | as 10619 plus a legacy operons file and the default prior set (genomic_context fired on 4 rows) | 5e6f66aa… |

Only the four M. tuberculosis cases are scorable on these tables (the Saccharolobus cases live in the
training cohort; the Sulfolobus acidocaldarius and S. elongatus cases are pending). Results and
interpretation: `docs/benchmarks/v3-role-truth-benchmark.md`.

## evaluations-registry-v2 (2026-10-28, after merging three drafted cases and one null group)

Same arms as above plus `claims-max` (simple method), each on the concatenated M. tuberculosis,
E. coli K-12 and B. subtilis 168 tables of that arm (`legacy`: bench9 `final/*_integrated.tsv`;
`rt1061x`: `role-truth-v1/jobs/1061x` outputs). Registry sha256 recorded in each `summary.json`.
