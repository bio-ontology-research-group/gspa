# Cohort pathway replay and reference-overlap audit

Latest status, 2026-10-27: **complete engineering replay, 18/18 genomes**. All
72 production CLI calls pass with the context-mode engine
`7f9d1d477098087469f90564d79a923ee28e85ba0a010c33840f2db44131671b`.
The complete audit verifies 175,272 reaction rows, 19,424 pathway rows and 56,613
raw candidates. There are 497 eligible reaction groups, 521 gap records (464
with GO), 7,448 RLC rows and 373 suggestions. Automatic and explicit gap routes
match across all 484,202 homology candidates; functional TSV bytes are unchanged.
The [complete replay archive](../../benchmark/v3/results/cohort-pathway-replay-v2/README.md)
retains all outputs and 338 bound inputs. This extends the frozen-output replay
below to the ten completed archaeal searches; it does not run the new absence
comparison or establish reference independence. No whole item is accepted.

Earlier checkpoint, retained below: **partial engineering verification, 8/18 genomes**. No whole
audit item is accepted. All 32 production CLI calls completed successfully, and
the final benchmark gate passes 600 tests with twelve optional skips. This work
does not change JVM code or the engine used by the preceding 1,221-test JVM gate
and nine live Nextflow checks.

The [evidence archive](../../benchmark/v3/results/cohort-pathway-replay-v1/README.md)
contains complete replay outputs, input/output hashes, exact commands, source
snapshots, tests, overlap findings and the failed initial overlap job. The
[acquisition checkpoint](v3-cohort-pathway-acquisition.md) supplies the actual
search provenance. At that checkpoint ten archaeal searches awaited reference
publication; they have since completed and are included in the replay above.

## Production behavior checked

`benchmark/v3/replay_cohort_pathways.py` first independently audits the search
attempts. It rejects an incomplete cohort by default; `--allow-incomplete` is an
explicit inspection mode whose receipt says `partial_verified` and retains all
ten missing genomes. No function labels are loaded for replay or verification.

For each completed genome, the annotation CLI receives a wrapper that validates
the complete query sequences, taxonomy and thread budget, then returns the two
unaltered observed tables. This is a **frozen-output replay**, not another search.
The independent Python oracle checks all raw candidate identities, roles,
individual hit rows and identity-derived scores. Complete raw tables and their
hashes must survive the production annotation route. Gap eligibility is checked
against complete reaction groups and the actual pathway-prediction table.

Three matched integration calls then use the same homology claims, real reaction
graph, EC aliases, EC-to-GO map and genome layout: no gap input, automatically
constructed gaps, and explicit gaps. Each call uses local reaction context with
no peer genomes. All candidate identities survive. Automatic and explicit routes
have identical full feature vectors, gap records, RLC feature TSVs and suggestion
TSVs. Relative to no gaps, only the three declared `ctx.gap.*` columns change.
Functional prediction TSVs remain byte-identical across all three arms.
Stage ledgers distinguish skipped no-gap RLC from executed RLC with gaps.

| Genome | Raw candidates | Eligible reaction groups | Gap records | Homology candidates | RLC rows | Suggestions |
|---|---:|---:|---:|---:|---:|---:|
| GCA_000008485.1 | 5,120 | 61 | 63 | 30,077 | 1,060 | 53 |
| GCA_000008765.1 | 4,758 | 26 | 31 | 36,403 | 420 | 21 |
| GCA_000009085.1 | 2,579 | 18 | 18 | 17,150 | 140 | 7 |
| GCA_000009985.1 | 8,059 | 70 | 73 | 42,570 | 1,039 | 52 |
| GCA_000012525.1 | 3,471 | 46 | 50 | 29,079 | 820 | 41 |
| GCA_000027305.1 | 2,595 | 34 | 35 | 21,721 | 560 | 28 |
| GCA_000046845.1 | 8,519 | 63 | 63 | 37,621 | 940 | 47 |
| GCA_000147055.1 | 7,694 | 54 | 56 | 50,066 | 660 | 33 |
| **Total** | **42,795** | **372** | **389** | **264,687** | **5,639** | **282** |

Of 389 gap records, 346 carry GO mappings. Records can represent alternatives
within one reaction group; these counts are not independent missing enzymes.
The full real-table oracle checks eligible groups, namespace and weak-hit flags;
it does not independently establish the biological correctness of each EC-to-GO
mapping. Existing JVM regressions cover multiple-association expansion.

Engine SHA-256:
`1f0a7da89894896b05626d386879be494fd1e8d1ecbacba3b048e600086befd6`.
Five replay regressions cover raw candidate conversion, complete-group gap
eligibility, missing/invented/wrongly flagged groups, malformed tables and
rejection of partial acquisition without explicit opt-in.

## Confirmed reference overlap

The eight reaction tables contain 1,295 rows reporting 100% identity and query
coverage, involving 411 genome-scoped target proteins. Alignment summaries alone
do not establish full-sequence identity. A frozen plan names 482 unique
genome/protein/reference pairs in 245 reference FASTA files. Slurm job **10133**
compares actual full sequences and rehashes query, reference and reported-hit
inputs before and after the scan. It completes with exit 0 on node004.

The first attempt, **10132**, failed immediately on node007 with scheduler
`RaisedSignal:53`, exit `0:53`, and no application log. Its scheduler record is
preserved; the cause is not established. The same frozen plan and code succeeded
when explicitly scheduled on node004. The failed attempt is not counted as a scan.

The scan confirms **466 exact pairs involving 395 target proteins**. Whitespace
is removed and sequence case normalized; terminal stop characters are retained.
The other sixteen proteins demonstrate why a rounded identity/coverage summary
must not be treated as proof of full-sequence identity.

| Genome | Frozen split | Proteins with an exact reference sequence |
|---|---|---:|
| GCA_000008485.1 | train | 28 |
| GCA_000008765.1 | test | 52 |
| GCA_000009085.1 | calibration | 57 |
| GCA_000009985.1 | test | 24 |
| GCA_000012525.1 | train | 92 |
| GCA_000027305.1 | calibration | 80 |
| GCA_000046845.1 | train | 45 |
| GCA_000147055.1 | test | 17 |

Thus exact reference overlap includes **93 test proteins**, 137 calibration
proteins and 165 training proteins. Splits are joined only after the sequence
audit; they do not select which pairs are inspected. This is a **lower bound**
from named reported hits. It does not scan every reference sequence, exclude
related families, inspect pretrained predictors or establish independence of
functional labels. These data cannot certify unseen-sequence or unseen-family
performance. Preserve them for execution checks; any independent accuracy
evaluation must separately enforce and audit its reference-exclusion protocol.

Four overlap regressions cover genome-scoped sequence identity, changed
references, changed reported-hit sources and conflicting duplicate reference IDs.

## Repeated search discrepancy

The earlier GCA_000012525.1 search (job 10111) and the cohort search (job 10125)
have identical proteome, tool, DIAMOND and all 41,274 reference-file hashes. They
ran on node004 with the same taxonomy, search arguments and four-thread budget.
Their recorder source and staging paths differ. The complete pathway tables are
byte-identical, but reaction rows increase from 13,147 to 13,149.

The preserved row-multiset comparison finds two added weak hits and two weak-hit
rows with changed reference/alignment values. This adds two raw EC candidates
(`ABB56303.1 / EC:2.1.1.133` and `ABB56352.1 / EC:1.1.1.135`), removes none, and
changes evidence for two existing candidates. Candidate counts are 3,469 versus
3,471. All eligible gap groups remain identical; the production gap outputs
contain the same 50 records after removing only the route-derived genome ID.

The cause remains unresolved. Inspection of the local upstream source finds
reference concatenation iterating Rust `HashMap` values, so query ordering is a
candidate explanation. This checkout is not proven to match the installed binary,
and ordering has not been isolated experimentally as the cause. Do not claim
search determinism, attribute the change to the GSPA adapter, or silently replace
the earlier observation. The archive preserves both the diagnostic and inspected
source. Frozen observed-output replay is deterministic for the checked routes.

## Remaining acceptance work

The archaeal acquisition and all eighteen-genome replays now complete. Matched
pathway, universal-network and organism-specific comparisons are recorded in
[the context-mode report](v3-context-modes.md), without evidence of a reliable
calibration benefit. Resolve or bound search reproducibility and enforce
reference independence for any generalization claim. These local-context replays are not
joint cohort inference, a trained-model comparison, calibrated locus/role/null
prediction or biological validation. RLC scores and gap assignments remain
heuristic hypotheses. Every whole audit item remains open in
[the acceptance record](../../tasks/v3-acceptance.md).
