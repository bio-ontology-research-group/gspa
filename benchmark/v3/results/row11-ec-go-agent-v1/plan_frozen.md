# DRAFT plan: EC-aware agent comparison on the Price genomes (acceptance row 11)

Status: FROZEN 2026-11-02. Drafted the same day after a feasibility round (cheap checks only, no new compute).
The SHA-256 of this file as frozen is recorded in `benchmark/v3/results/row11-ec-go-agent-v1/PLAN_SHA256`
before any score or label is read; later changes are appended under "Deviations", dated; the body above that
heading is not edited.

Approvals (Robert Hoehndorf, 2026-11-02, relayed by the coordinator): (1) the locus-anchored Universe B with
2:1 decoys, disclosed as written; (2) the floor 0.20 and max-10-candidates CLI options, defaults unchanged;
(3) the `ec6+contract` set is included once the other session's DIAMOND contract claims are frozen; `ec6` runs
first if they are not ready, recorded as a dated deviation; (4) the Universe A cap is left to this plan's
author (50 per genome if agent wall-clock needs it). Implementation by agy-auto from specs, reviewed here.

Row 11 bar: "Public evidence separated from truth before selection; opaque IDs; no invented scores; identical
inputs/guardrails/budgets for simple and learned methods; unresolved-only calls", measured on "unseen genomes /
families; real unresolved cases; calibration-frozen decisions; correctness / harm / abstention and complete
inference / tool cost". The missing element is correctness and harm: no real unresolved case so far has
experimental truth (`docs/benchmarks/v3-stage3-production-comparison.md`, `v3-roles-agent-closure.md`).

Robert's direction (2026-11-02): "EC and GO should go together; many EC numbers are a part of GO anyway and can
be aligned by mappings in go.obo. Otherwise, the agent should get EC numbers, not GO (or rather: in addition)."

## 1. What the feasibility round measured (2026-11-02, ws, read-only)

Sources: the Price stage-1 integrations of TopEC, MAPred and EnzBERT from the roles closure
(`/mnt/data2/leechuck/gspa-v3-data/roles-agent-closure-v1/agent/<predictor>/stage1/`, 22 genomes, EC rows only),
the retained release-candidate scores of all six predictors at the 148 Price loci
(`price-calibration-context-v1/inputs/rc1-evaluation/retained/`, 132 cells), the Price reference files
(`ec2go.txt` 2026-03-25, `kegg_pathways.tsv` 5,262 steps, `kegg_modules.tsv` 2,822 steps), and
`~/.cache/gspa-v3-data/training-cohort/go-basic.obo` (2026-07-26). Scripts (feasibility only, not pipeline code):
`~/tmp/gspa-v3-agy/row11-feasibility/{simulate_ec_cases.py,locus_status.py}`.

### 1.1 EC to GO alignment

| quantity | value |
|---|---:|
| GO terms in go-basic.obo (2026-07-26) | 48,340 |
| GO terms carrying at least one `xref: EC:` | 4,662 (4 obsolete) |
| EC xref lines (distinct ECs) | 4,835 (4,578 complete, 257 partial `-`) |
| ECs in the pipeline's `ec2go.txt` (2026-03-25; one GO per EC, no duplicates) | 4,822 (4,565 complete, 257 partial) |
| Price truth: distinct ECs / (protein, EC) pairs | 56 (all complete) / 152 over 148 loci |
| truth ECs found in ec2go / in go.obo xrefs | 56 / 56 (pairs: 152 / 152) |
| EC rows in the Price integrations with a partial EC | 0 (all six predictors emit four-level ECs) |

`ec2go` is generated from the go.obo xrefs, so the two agree on every truth EC; the pipeline already loads
ec2go into `IntegrationState.ecToGoAlternatives` (`IntegrationReferenceLoader`), uses it for the pathway
database (`PathwayGraph.getRequiredGoTerms`, `ecToGo[ec]`), for the gap-filling prior and for taxon
consistency. What it does not do: the integrated table keeps EC claims as `EC` rows (no GO row is derived), and
`AnomalyExtractor.extract` builds its annotated set from `GO` rows only and skips non-GO keys for borderline
conflicts. That is the whole reason the Price collections yield 0 cases.

### 1.2 Where the truth loci sit in the production outputs

Three local predictors, stage 1, max-merged over predictors, per (protein, EC) truth pair
(calibration = the 7 genomes of the frozen genus-blocked split, seed 20261101; test = 15 genomes):

| truth pair status | calibration (72 pairs) | test (80 pairs) | Swiss-Prot experimental (35, all test) |
|---|---:|---:|---:|
| called (posterior >= 0.5) | 42 | 35 | 25 |
| candidate 0.2–0.5 | 9 | 15 | 0 |
| candidate 0.05–0.2 | 5 | 8 | 1 |
| absent | 16 | 22 | 9 |

All six predictors (retained release-candidate scores at the 148 loci, max over predictors, stage 1), rank of
the truth EC among the protein's own EC candidates:

| truth EC rank on its protein | calibration | test | total |
|---|---:|---:|---:|
| rank 1 | 19 | 13 | 32 |
| rank 2–8 | 35 | 43 | 78 |
| rank 9–20 | 5 | 9 | 14 |
| below 20 | 13 | 14 | 27 |

Rank-2–8 loci have a median top posterior of 0.60 and a median truth posterior of 0.55: these are proteins with
several complete ECs called at once, that is, genuine EC conflicts, not sub-threshold gaps. Only one locus has
the truth at rank 1 below 0.5.

### 1.3 Dry run of an EC-aware extractor (truth-blind case universe)

Pathway gaps keyed through ec2go (present = annotated GO or a descendant at >= 0.5; candidates = proteins with
posterior in [floor, 0.5) on an EC of the missing step; top 8 kept), all 22 genomes, before the 150-case cap:

| collection | floor | gap cases (KEGG maps) | cases with a labelled candidate | positive (protein, EC) pairs | gap cases (modules) | positives (modules) |
|---|---:|---:|---:|---:|---:|---:|
| TopEC | 0.05 | 43 | 0 | 0 | 416 | 0 |
| MAPred | 0.05 | 6,332 | 174 | 4 | 5,681 | 3 |
| MAPred | 0.2 | 2,522 | 49 | 4 | 2,389 | 3 |
| EnzBERT | 0.05 | 3,019 | 75 | 7 | 4,748 | 10 |
| EnzBERT | 0.2 | 233 | 1 | 0 | 399 | 0 |
| three merged | 0.05 | 8,397 | 273 | 3 | 5,715 | 4 |
| three merged | 0.2 | 3,163 | 56 | 2 | 2,297 | 4 |

EC conflicts over the whole proteomes (three merged, 115,535 proteins with a candidate): two complete ECs both
>= 0.5 on one protein: 20,428 proteins (6,615 of them sibling ECs sharing the first three digits); called top
with a second within 0.10: 23,286; borderline top (0.35–0.65) with a close second: 30,392. On the 148 labelled
proteins: 128 conflicts, 58 of which involve the truth EC.

Two facts follow. First, the engine's provisional calibrator sigmoid(4·raw − 2) maps raw 0 to posterior 0.12,
so the extractor's floor 0.05 admits every exported candidate (EnzBERT exports 100 ECs per protein, CLEAN and
GraphEC about 5,400); the floor has to sit above the calibrator's zero point. Second, and decisive: the base
rate of labelled loci is 148 in about 115,000 candidate proteins (0.13 %). A truth-blind, budgeted case sample
of a few hundred cases touches at most a handful of truth loci (3–10 positive pairs in all 22 genomes in the
dry run, and that is before any budget). A truth-blind universe therefore cannot produce the correctness and
harm numbers; only a locus-anchored universe can, and it has to be declared as such.

### 1.4 Agent backend and cost of record

`evaluate_stage3_agent.py agent` calls `gspa-core/src/main/resources/agent/run_gspa_agent.py`: an
OpenAI-compatible chat endpoint (BORG LLM gateway, `GSPA_LLM_URL`, default
`http://unimatrix01.kaust.edu.sa:8000/v1`, key from `BORG_LLM_API_KEY` or `~/.config/borg-llm/api-key`; the key
file exists on ws and the endpoint answers, with an authentication error without the key), model `qwen3.8-27b`,
temperature 0, seed 42, thinking disabled, `max_tokens` 450 per turn, two turns, four workers, five local
inspection tools. Cost of record (196 cases): 392 requests, 730,626 prompt and 71,169 completion tokens,
589.9 s wall; per case 2.0 requests, about 3,700 prompt and 360 completion tokens, 3.0 s wall at four workers.
No billing applies (self-hosted); the token counts are the cost.

## 2. Design

### 2.1 EC to GO alignment (one mapping, used everywhere)

- Mapping of record: the pipeline's `ec2go.txt` (the Price reference copy, sha256 recorded), which is generated
  from the go.obo xrefs. A check script verifies at run time that every EC in ec2go is an `xref: EC:` of the same
  GO term in the pinned `go-basic.obo` and lists the differences (13 ECs newer in the 2026-07 obo); the truth ECs
  must be in both (currently 56/56).
- Direction: EC to GO is one-to-one in ec2go; GO to EC is one-to-many (one GO term can carry several ECs, and
  partial ECs map to broader terms). The extractor therefore keys cases on the GO term of record and carries the
  EC list of that term; a case built from GO-bearing claims (Step G) and a case built from EC claims land on the
  same key.
- No GO rows are derived into the integrated table (no invented scores); alignment happens inside the extractor
  and in the truth join only.

### 2.2 EC-aware extractor (JVM change, `gspa.agent.AnomalyExtractor`)

1. Annotated set: `annotatedGo` is built from GO rows and from EC rows mapped through
   `state.ecToGoAlternatives` (already loaded). Pathway completeness then sees EC-only collections.
2. Candidates: for a missing GO term, candidates come from GO provenance of that term and from EC provenance of
   every EC mapping to it; each candidate row records `function_id` as claimed (EC or GO), the mapped GO term,
   and the usual posterior, log-odds, supporting claims, operon and neighbour fields.
3. New case types, all keyed on one protein ("locus cases"), target type `EC`:
   - `EC_CONFLICT`: two or more distinct complete ECs at or above the annotated threshold on one protein;
     candidates are the called ECs plus any further EC above the floor; context: shared first three digits,
     pathway membership of each EC (KEGG maps and modules), operon neighbours' called ECs/GO.
   - `EC_BORDERLINE`: top EC in [0.35, 0.62] with a second EC within 0.15 or with multi-source support or a fired
     prior (the existing borderline rule without the GO-only filter).
   - `EC_PARTIAL`: a called partial EC (`x.y.z.-`) with complete children above the floor (none occur in the
     six predictors' outputs; kept for GO-bearing/InterPro claims).
   - `PATHWAY_GAP` as today, now reachable through EC rows (2.2.1–2).
4. Floor: `subThresholdFloor` becomes a posterior floor of 0.20 (raw about 0.15 under the provisional
   calibrator) with the rationale in 1.3; `maxCandidatesPerAnomaly` 10 for locus cases. Both are CLI options
   (`--anomaly-floor`, `--anomaly-max-candidates`) so Step G's frozen runs are unaffected (defaults unchanged
   unless given).
5. Budget: `maxAnomalies` stays 150 per genome; within the cap, locus cases are ordered by a truth-free priority
   (closeness of the two best posteriors, then multi-source support, then pathway membership), pathway gaps by
   completeness. The order is recorded in the case.
6. Operons: pass derived operons (production distance rule from the layouts, as `roles_closure/operon_context.py`
   does) so `in_pathway_operon` and neighbour functions are populated; the roles-closure Price runs had none.
7. Tests (Spock): EC-only collection yields a pathway gap through ec2go; an EC conflict case; partial EC case;
   floor and cap options; GO-only behaviour unchanged with default options (byte-identical anomalies on an
   existing fixture).

The JAR is built to a private copy (`row11-ec-go-agent-v1/jar/`, sha256 and git head recorded); the worktree
JAR used by the other session is not rebuilt.

### 2.3 Inputs and integration cells (Shaheen)

- Claims: the six predictors' whole-proteome Price claims (`/data/hohndor/gspa-v3/price-full-matrix-v2/parts/`
  on unimatrix: CLEAN 116 GB, GraphEC 113 GB, DeepProZyme 62 GB, EnzBERT 2.3 GB, MAPred 397 MB, TopEC 28 MB; the
  Shaheen /project mirror at the `/data` paths), merged per genome into one claims file (the ensemble_pool4
  convention), every part digest re-verified against `part.json`.
- Two claim sets: `ec6` (six predictors) and `ec6+contract` (six predictors plus the contract-bearing
  Swiss-Prot DIAMOND claims produced by the other session's Step C, read from its published, hashed claim files
  once they exist; never regenerated or edited here). The second set is what gives the learned models anything
  to score (they abstain on contract-less claims by design).
- Arms per claim set: `stage1` (default priors), `learned_logistic`, `learned_gbdt`; same layouts, taxon,
  essential profile bacteria, Price reference files, derived operons, `--anomalies-out`, EC-aware extractor with
  floor 0.20; receipts and schema-2 ledgers through `v3_runner.run_arm`. 22 genomes × 2 sets × 3 arms = 132
  cells, one Slurm array, `--exclude=nid00022,nid00497,nid01033`.
- Full-matrix merges are large (about 1.5 M rows per genome); the roles-closure engine handled 0.86 M rows in
  9 s with 40 GB heap, so 64 GB heap, one genome per task.

### 2.4 Two case universes, truth separated first

Universe A (truth-blind, production): the extractor's output as is, capped at 150 per genome. Expected: a few
thousand cases before the cap, 3,300 after it, with 0–10 touching a truth pair (1.3). It is run and scored
because it is the production path, and its (near-)zero truth intersection is reported as the measured base rate.

Universe B (locus-anchored, disclosed): the case universe is the set of proteins with experimental EC truth
(148 Price loci; 35 Swiss-Prot experimental, 15 of them circular for the DIAMOND arm and reported separately),
together with truth-free matched decoys. Selection uses only the locus identity, never the EC label; the locus
list is public (Price et al. 2018 supplementary; Swiss-Prot). Construction:

1. Freeze the locus list (sha256) and the decoy rule before any integration output is read.
2. For every listed protein, the EC-aware extractor's locus case is taken from the production output if the
   protein is unresolved under the frozen definitions (EC_CONFLICT, EC_BORDERLINE, or no EC at or above 0.5 with
   at least one candidate above the floor). Resolved loci (one EC called, no competitor) and no-candidate loci are
   not cases; they are counted and reported (a call cannot be scored there; "never correct").
3. Decoys: for each labelled case, two proteins of the same genome in the same case type, nearest in top
   posterior, not in the locus list, not sharing an orthogroup with a labelled protein (cohort orthogroups).
   Decoys carry no truth; they are scored as unknown.
4. The public bundle (`evaluate_stage3_agent.py prepare` format: opaque case IDs, candidates, evidence, learned
   scores where present, no truth-bearing field) is written and hashed. Only after that hash is recorded is any
   label file opened (`row11_truth` join, 2.6).
5. Leak check, pre-registered: per method, the call rate on labelled cases and on decoys must be
   indistinguishable (two-sided Fisher, p > 0.05, and absolute difference < 10 points). A method that calls
   more on labelled cases has an information channel; its results are reported as invalid.

Universe B conditions on "this protein has an experimentally established enzyme activity", which enriches for
enzymes that earlier annotation got wrong. It does not condition on the answer. The report states this in the
first paragraph and in the acceptance row.

Expected scorable cases from 1.2: truth within the top 10 candidates on 110 of 151 pairs (54 calibration, 56
test), of which about 78 are conflicts with the truth at rank 2–8 and 32 have the truth on top; the single
three-predictor merge gives 98 of 152 with the truth as a candidate at >= 0.2. The merged six-predictor
integration will move these numbers (log-odds combination raises agreeing calls); the measured counts replace
the estimates. Swiss-Prot adds about 24 test cases (all in test genomes), 15 of them circular for the DIAMOND arm.

### 2.5 Methods under identical inputs, guardrails and budgets

Same public bundle, same guardrail (`SymbolicGuardrail`: GO ontology, taxon constraints; EC targets pass the EC
identity check and the same final-candidate constraints), same decision vocabulary per candidate:
`PROMOTE` (assert the EC), `DEMOTE` (remove a called EC), `ABSTAIN`; `RELAX` only toward a GO ancestor, as today
(an EC is relaxed through its ec2go term).

| arm | source of the decision | learned scores available |
|---|---|---|
| `abstain_only` | none | n/a |
| `rank_by_raw_seq` | best raw predictor score | n/a |
| `rank_by_stage1_posterior` | stage-1 posterior | n/a |
| `rank_by_context_support` | `in_pathway_operon` | n/a |
| `learned_logistic`, `learned_gbdt` | production learned probability, where the candidate carries a contract | only on `ec6+contract` candidates with a DIAMOND claim; otherwise abstain (measured) |
| `qwen_agent` | `run_gspa_agent.py`, two turns, temperature 0, seed 42, four workers, 450 tokens per turn, five tools over the bundle only | same bundle |

Budgets: the agent's turn, token and tool limits are those of the cost-of-record run; simple and learned
methods are deterministic over the same fields. No method sees labels, decoy flags or locus-list membership.

Policies: the existing `calibrate` rule, fixed before any test outcome: at least ten known calibration
decisions and at least 90 % demonstrated precision; maximise correct, then minimise harm, then the highest tied
threshold; otherwise abstain. Change for Universe B, fixed now: precision is computed over known decisions
(labelled cases); calls on decoys are reported as the unknown-call rate, not counted against precision, because
decoys are two thirds of the bundle by construction. Universe A keeps the original rule.

### 2.6 Truth join and scoring

Truth: Price-148 complete ECs (`<assembly>.ec.tsv`, exact protein ID) and the Swiss-Prot experimental Price
ECs (exact sequence), joined after the public hash. EC targets match exactly; GO targets are positive when
ec2go maps a labelled EC of that protein to the target term (the other session's rule, reused). No negative is
derived: a decision naming another EC on a truth protein is `unconfirmed`. Harm is defined only as demoting or
removing an experimentally supported EC that the production output had called or carried as a candidate, or
promoting a different EC while demoting the supported one.

Per method and split: correct, unconfirmed, harmful, abstained, unknown (decoys), failed; no-candidate loci
separately; Wilson 95 % intervals; paired bootstrap over loci for agent versus each simple method and versus
learned; recall of the truth EC among unresolved loci; cost (requests, tokens, decision wall time, shared
production time).

Pass criteria for the row (stated now, reported as measured):

- C1 extractor: EC-only collections yield cases; GO-only outputs unchanged under default options (tests pass;
  anomalies byte-identical on the fixture).
- C2 separation: public hash recorded before any label file is opened; leak check passed by every reported
  method.
- C3 power: at least 40 known unresolved test loci and at least 10 known calibration decisions available to at
  least one method; otherwise the row stays Not met with the measured counts and the comparison is reported as
  underpowered (with n = 40, a Wilson interval has half-width about 0.15 at p = 0.5, and a paired contrast can
  show only differences above about 25 points).
- C4 evidence for the bar element: at least one method makes at least one known test decision under its frozen
  policy; correctness, harm and abstention are then reported for all seven methods with intervals.
- Row 11 moves to Met only if C1–C4 hold on Universe B and Universe A's truth intersection is reported; a
  negative result (every method abstains, or the agent is not better) still satisfies the bar if measured.

### 2.7 Cost and wall-clock

- Agent: Universe A 3,300 cases (two sets × stage1 only, agent runs once per public bundle, the three arms share
  the agent run where candidates are identical) and Universe B about 390 cases (130 labelled + 260 decoys) per
  set, so about 4,100 cases, 8,200 requests, about 15 M prompt tokens and 1.5 M completion tokens, about 3.5 h
  wall at four workers on the unimatrix endpoint (3.0 s per case). If the endpoint is slower than in the run of
  record, Universe A is capped at 50 cases per genome (1,100 cases, 1 h); Universe B is never cut.
- Shaheen: 132 cells, about 10–20 min each with full-matrix merges, one node per task; well inside the budget.
- ws: merges, bundle preparation, truth join, simple methods, evaluation: light Python.
- No monetary cost: self-hosted model. Token totals are reported as the cost.

### 2.8 How this combines with the other session's Step G (`tasks/v3-price-contract-learned-2026-11-02.md`)

Step G runs the unchanged extractor on GO-bearing DIAMOND contract claims and joins truth with the ec2go rule.
This plan reads Step G's published claims (not its results) as the `+contract` component, so that the learned
models have contract-bearing candidates inside EC-aware cases, and reuses its truth rule and `row11_truth`
format. Nothing in Step G's folders is edited. The row-11 report cites both: Step G for GO-bearing claims with
the production extractor, this plan for EC-aware cases and the locus-anchored measurement. If Step G already
finds known test decisions, its numbers are the production-universe evidence and Universe B is the
locus-anchored complement; if it finds none, Universe B is the only source of correctness and harm.

### 2.9 Risks

- Base rate: Universe A will almost certainly intersect truth on fewer than ten pairs; it is reported, not
  relied on.
- Disclosure: Universe B conditions on label existence; a reviewer may call this selection. Mitigations are
  decoys, the leak check, opaque IDs, and the explicit statement; the acceptance row records the limitation.
- Calibrator: posteriors are compressed (raw 0.8 maps to 0.77, raw 0 to 0.12); the EC-conflict and borderline
  definitions are stated on posteriors of the production engine, not retuned on Price.
- Learned arm: can decide only where DIAMOND contract candidates exist; most EC-conflict candidates have none,
  so the learned methods will largely abstain. That is a true property of the shipped models and is reported.
- Pretraining overlap of the six predictors with Price proteins is unassessed (stated in `v3-price148.md`).
- Dependence: several loci per genome and per pathway; intervals are per locus, and genome-blocked bootstrap
  is reported alongside.
- Circular Swiss-Prot labels (15, S. meliloti) for the DIAMOND arm: reported separately, as in Step G.
- Endpoint availability: the agent run is retried once; an unavailable endpoint is reported as not applicable,
  as in Step G's G3.

## 3. Work plan

| step | where | who | output |
|---|---|---|---|
| S0 freeze this plan, locus list, decoy rule, split | ws | me | `PLAN.sha256`, `LOCI.sha256` |
| S1 ec2go versus go.obo check | ws | agy-auto (spec `~/tmp/gspa-v3-agy/row11-ec2go-check/spec.md`) | `ec2go_obo_check.json` |
| S2 EC-aware extractor + Spock tests + CLI options, private JAR | ws build | agy-auto (spec `row11-extractor/spec.md`), I review | JAR sha256, test log |
| S3 claim merges (`ec6`, `ec6+contract`), derived operons | Shaheen | agy-auto (spec `row11-cells/spec.md`) | parts, digests |
| S4 132 integration cells | Shaheen | sbatch `benchmark/v3/shaheen/row11-ec-go-agent.sbatch` | receipts, anomalies |
| S5 bundles for Universe A and B, decoys, public hash | ws | agy-auto (spec `row11-bundle/spec.md`) | `public_*.json`, `PUBLIC.sha256` |
| S6 truth join (after S5 hash) | ws | reuse of `row11_truth.py` logic in a copy under `benchmark/v3/row11_ec_go/` | `truth_*.json` |
| S7 six simple methods, learned, agent; calibrate; evaluate; leak check | ws + endpoint | agy-auto (spec `row11-methods/spec.md`) | proposals, policies, reports |
| S8 report `docs/benchmarks/v3-row11-ec-go-agent.md`, acceptance row text | ws | me | report, `SHA256SUMS` |

Results folder: `benchmark/v3/results/row11-ec-go-agent-v1/`; raw outputs under
`/mnt/data2/leechuck/gspa-v3-data/row11-ec-go-agent-v1/` (ws) and Shaheen scratch (`LARGE_FILES.tsv`).
No commits in this round. `tasks/v3-acceptance.md` is edited only in S8.

## Deviations (append only, dated)

(none at freezing)
