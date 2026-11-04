# EC-aware agent comparison on the Price genomes (acceptance row 11)

DRAFT: numbers are filled in from `benchmark/v3/results/row11-ec-go-agent-v1/tables.md` when the collection run
completes; this file is not final until the "Status" line below says so.

Status: in progress (2026-11-03).

Frozen plan: `tasks/v3-row11-ec-go-agent-2026-11-02.md` (sha256 of the frozen text `72556ac7…` in
`benchmark/v3/results/row11-ec-go-agent-v1/PLAN_SHA256`; the plan with its dated deviations `c6abff8a…` in
`PLAN-with-deviations.sha256`). Results folder: `benchmark/v3/results/row11-ec-go-agent-v1/`. Raw runs, bundles'
private role maps, truth files and integrated tables: `/mnt/data2/leechuck/gspa-v3-data/row11-ec-go-agent-v1/` on
ws and `/scratch/project/k10018/gspa-v3/out/row11-ec-go-agent-v1/` on Shaheen. Code: `benchmark/v3/row11_ec_go/`,
`benchmark/v3/shaheen/row11-ec-go-agent.sbatch`, tests `benchmark/v3/tests/test_row11_ec_go*.py`; the JVM change
(EC-aware anomaly extractor) lives in the worktree branch `agent/row11-ec-extractor` (patch and JAR with sha256 under
the data folder's `jar/`).

## 1. Question and bar

Row 11 asks for an agent comparison on real unresolved cases with experimental truth, under identical inputs,
guardrails and budgets for simple, learned and agent methods, with calibration frozen before test scoring, and with
correctness, harm, abstention and complete cost reported. Until now no real unresolved case had experimental truth,
and the only collection with experimental enzyme truth on new genomes (the 22 Price genomes, 148 reannotated loci)
produced no case, because its outputs are EC numbers and the production extractor keyed every case on GO terms.

Robert's direction was to align EC and GO through the go.obo cross-references and to give the extractor and the agent
EC numbers in addition to GO. This report measures the result.

## 2. What was built

### 2.1 EC to GO alignment

The pipeline's `ec2go.txt` is the mapping of record; it is generated from the go.obo `xref: EC:` lines, and a check
script compares the two (`ec2go_obo_check.json`).

| quantity | value |
|---|---:|
| ECs in ec2go (2026-03-25) | 4,822 |
| GO terms with an EC xref in go-basic.obo (2026-07-26) | 4,662 |
| ECs in both files / mapped to the same term | 4,819 / 4,819 |
| ECs only in ec2go / only in the obo | 3 / 16 |
| Price truth ECs (distinct) in ec2go / in the obo | 56 / 56 |

### 2.2 EC-aware extractor

`gspa.agent.AnomalyExtractor` gains an `ecAware` mode (CLI `--anomaly-ec-aware`, off by default; default output is
unchanged, verified by a test that compares the default output with and without EC rows). In that mode EC rows at or
above the annotated threshold count toward the annotated GO set through ec2go, candidates for a missing GO term come
from GO rows and from EC rows of the ECs mapping to it, and a second file `*.anomalies.locus.jsonl` holds locus cases:
EC_CONFLICT (two or more complete ECs called on one protein), EC_PARTIAL (a called partial EC with complete children)
and EC_BORDERLINE (top EC in [0.35, 0.62] with a close second, multi-source support or a fired prior). The candidate
floor and the candidate and locus caps are CLI options (`--anomaly-floor 0.20`, `--anomaly-max-candidates 10`,
`--anomaly-max-locus`).

### 2.3 Inputs

The six enzyme predictors' whole-proteome Price claims (CLEAN, DeepProZyme, EnzBERT, GraphEC, MAPred, TopEC) merged
per genome after a uniform raw-score floor (`raw_score >= 0.05`, see the deviation of 2026-11-02 and the gate in 3.1),
as set `ec6`; set `ec6c` adds the other session's frozen Swiss-Prot DIAMOND contract claims. Derived operons
(production distance rule), the Price layouts, taxon, essential profile bacteria and the Price reference files are the
same as in every Price analysis. Arms stage1, learned_logistic, learned_gbdt through `v3_runner.run_arm` with receipts.

### 2.4 Two case universes

Universe A is the production extractor's target-centred output (truth-blind, capped at 150 per genome by the
extractor and at 50 per genome in the bundle for agent cost). Universe B is locus-anchored: the frozen list of 180
proteins with experimental EC truth (148 Price loci, 32 Swiss-Prot experimental; identities only, `loci.tsv`, sha256
frozen before any run) is joined to the locus file, every unresolved listed protein becomes a case, and each labelled
case is accompanied by two truth-free matched decoys (same genome and case type, nearest top posterior, no shared
orthogroup with a listed protein). Opaque case IDs; the public bundles carry no label, role or source field; their
sha256 is recorded in `PUBLIC_SHA256` before any label file is opened. Universe B conditions on the existence of an
experimental label, not on its value; the leak check (call rate on labelled versus decoy cases) tests whether any
method could tell the two apart.

## 3. Gates and deviations

### 3.1 Input-reduction gate

TABLE_FULLCHECK

### 3.2 Learned stage

With the merged claims the engine's learned stage skipped every genome in both sets ("no candidates supported by the
model input-coverage manifest"), the DIAMOND contract claims notwithstanding: the models admit only candidates whose
source and availability profile was seen in training, and six additional EC sources make every profile unseen. The
learned methods therefore abstain by construction on both sets; all methods decide on stage-1 scores. The learned
arms' receipts record the skip.

### 3.3 Locus cap

The first run hit the 150-case locus cap in every genome with EC_CONFLICT cases alone; the cap was lifted for the
second run (`--anomaly-max-locus 1000000`), the first run's receipts are archived and unused.

## 4. Results

### 4.1 Cases

TABLE_CASES

### 4.2 Universe B coverage of the labelled loci and truth among candidates

TABLE_COVERAGE

### 4.3 Methods under frozen policies (Universe B)

TABLE_METHODS_B

### 4.4 Universe A

TABLE_A

### 4.5 Contrasts

TABLE_CONTRASTS

### 4.6 Cost

TABLE_COST

## 5. Verdict against the frozen criteria

TABLE_VERDICT

## 6. Limitations

LIMITATIONS

## Files

FILES
