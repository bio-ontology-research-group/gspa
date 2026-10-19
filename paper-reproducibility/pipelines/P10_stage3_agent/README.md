# P10 — Stage 3 (guarded language-model agent): benchmark v2, LRO-50, trigger rates

Regenerates the Stage 3 benchmark from **real integrator anomalies** with committed code,
re-runs every agent configuration with full transcripts, recomputes all baselines, the
LRO-50 re-ranking and the Stage 3 trigger rate over the 88 multikingdom proteomes.
Benchmark v1 (`results/pillar4_agent_benchmark_120.jsonl`) is superseded: it had no
generating code and leaked labels (see "Why v2").

Compute: `gspa integrate` runs on unimatrix01 (SLURM `debug`) and, on 2026-10-19, on `ws`
(leechuck-office, 56 cores) against a byte-identical mirror of the inputs (`ws_env.sh`),
because every debug node was allocated by other jobs and node007 has no `/data` mount.
LLM calls run from any host that reaches the vLLM gateway. Remote run directory:
`/data/hohndor/gspa/proteomes/pillar4_repro_2026-10-19/` (unimatrix01).

## Software and model

| Item | Value |
| --- | --- |
| Jar | `/data/hohndor/gspa/bin/gspa-unified.jar`, `gspa 2.1.0`, sha256 `e415d9c6e13fbf96382ba7d798074b2ca4ec2991499de6714a5e4da2b19f72a7` |
| Driver | `p10_integrate.groovy` (runs the jar's `IntegrateCommand` unmodified; only `AnomalyExtractor.maxAnomalies` 150 -> 100000; production output = first 150 records) |
| JDK | `/data/hohndor/envs/jdk-21.0.10+7` (same build copied to the ws mirror) |
| GO | `reference/go.owl` sha256 `5ba93c63…82d`, `reference/go.obo` (releases/2026-03-25) sha256 `a77e3567…4bd`; guardrail OBO `benchmark/neural/cafa6_recon/go-basic.obo` (releases/2025-06-01) sha256 `ca48860e…00b` |
| Pathways | `kegg_pathways.tsv` `f112964c…82d`, `kegg_modules.tsv` `27184381…8ea`, `ec2go.txt` `98ca99f8…abd` |
| Sidecar | `benchmark/agent/run_gspa_agent.py` sha256 `250d878a1b60e556891dcaaeb7c94c7c7c2a05a0953edf5dce57ece27bcb00ae` (imported unmodified by `run_pillar4_v2_agent.py`; the sha256 is also in every transcript header) |
| Model | `Qwen/Qwen3.8-27B-FP8` on the vLLM gateway `http://unimatrix01.kaust.edu.sa:8000/v1` (alias `qwen3.8-27b`, max_model_len 131072), `chat_template_kwargs.enable_thinking=false`, `max_tokens=450` per turn |
| Configurations | c0 zero-shot (1 turn, only `propose_resolution`, no guardrail); c2 fast (2 turns, guardrail); c3 guarded (6 turns, guardrail) seeds 42 (T=0.0), 43 (T=0.25), 44 (T=0.25); c4 unguarded (6 turns, no guardrail, seed 42, T=0.0) |
| API key | Read by the sidecar's `load_api_key()`: env `BORG_LLM_API_KEY` or `OPENAI_API_KEY`, else `~/.config/borg-llm/api-key` (mode 600). Never printed, never written to transcripts (the recorder stores the request body only) |

## Genomes (`genomes_v2.tsv`)

| Genome | Taxon asserted | Claims | Truth |
| --- | --- | --- | --- |
| ecoli, bsubtilis, paeruginosa, mtb | NCBITaxon_562 / 1423 / 287 / 1773 (as in the 88-proteome manifest) | `bench9/with_iprscan/{tag}_claims_merged.jsonl` (DIAMOND + InterProScan + Pfam); operons `operons/{tag}_operons.tsv`; gaps `bench9/gapseq/{tag}_gaps.jsonl` | `truth_dual/{tag}_truth_all.tsv` |
| scerevisiae, celegans, athaliana | NCBITaxon_4932 / 6239 / 3702 | `multikingdom_batch/blind_claims/{sp}_claims.jsonl` (DIAMOND + Pfam) | `multikingdom_batch/blind_truth/{sp}_truth_all.tsv` |

Integrator configuration for every run (`integrate_one.sh`) = the 88-proteome `iter3_relax`
configuration: theta `max_iter=3`, all priors, `--go-owl --lite`, KEGG pathways + modules +
ec2go, essential profile, asserted taxon, Stage 2 consistency enforcement (`relax`), operons,
gaps; when a gap table is present the dark-matter suggester is on (its suggestions are
extractor candidates). Truth support uses the is_a/part_of closure of the truth annotations.

## Benchmark v2 design (leave-annotation-out; same extractor path and schema for positives and negatives)

* **PATHWAY_GAP** — pathways/modules whose required GO steps are all annotated (posterior
  ≥ 0.50, closure) and truth-supported; a step T with a single covering protein p is held out
  by deleting p's claims for T and every descendant of T (`tp`; "delete that protein's claims
  for the target function"), the held-out step is appended to the gap table, `gspa integrate`
  is re-run, and the extractor's own PATHWAY_GAP record for (pathway, T) is taken. **TP** = p
  is among the extractor's candidates (it can only re-enter through the dark-matter
  suggester); expected PROMOTE p. **Decoy** = same case, p's claims deleted entirely and p
  removed from the operon table (the gene is absent); expected ABSTAIN/DEMOTE.
  Eligible TP cases (50 of 597 hold-outs) were split by sha256 order into a TP half and a decoy
  half *before* any decoy was run (`plan`), so TPs and decoys come from the same population.
  Decoys were matched 1:1 to TPs on maximum candidate posterior (±0.05) and candidate count.
* **ONTOLOGY_CONFLICT** — per genome 40 real proteins receive a GO term the guardrail forbids
  for the organism (pool: GO terms of the bundled GO-taxon and GO-anatomy constraint tables;
  seed 42) as two claims from the genome's evidence sources (raw scores on a 0.40–0.85 grid).
  **Genuine** = the extractor's BORDERLINE_CONFLICT record for an injected (protein, term)
  (posterior in [0.35, 0.50)); expected RELAX to a valid ancestor or DEMOTE. **Trap** = a
  BORDERLINE_CONFLICT record from the same runs whose term is truth-supported for the
  candidate and valid for the organism; expected ABSTAIN (or PROMOTE = keep). Traps matched
  1:1 on candidate count, supporting-claim count, posterior (±0.05), same genome preferred.
  These runs use **Stage 1 only**: with Stage 2 on, all 280 injected violations were relaxed
  or removed before the extractor ran (0 reached Stage 3; `holdouts/*/*_CONFSTAGE2_r0`).
* **INCOMPLETE_COMPLEX / BORDERLINE_ESSENTIAL — not constructible.** The extractor builds
  candidates only from sub-threshold provenance of the *exact* complex / essential term.
  Over the 7 genomes the selection found 7 complex and 32 essential hold-outs; under the `tp`,
  `tpw` (weakest claim kept) and `tph` (only DIAMOND claims deleted) constructions 1, 0 and 0
  of them produced an extractor record containing the held-out protein (most produced no
  record at all). The 88-proteome run shows the same scarcity in natural data (6 complex and
  83 essential anomalies, uncapped, against 31,839 anomalies in total). Both categories are
  therefore absent from v2 rather than filled with constructed records.
* The `tpw` construction was run for PATHWAY_GAP (82 eligible) but rejected: it leaves the
  held-out protein with exactly one supporting claim while decoy candidates have none
  (leak-check AUC 0.98 for "a candidate has a single claim"). `assemble` output of that draft
  is not used.
* **De-identification** — every protein id the model can see (candidates, operon neighbours,
  anchors) is replaced by an opaque token `P####`, operons by `O####`, anomaly ids by `A###`;
  tokens and anomaly ids are assigned in sha256 order, candidate order is sha256-shuffled per
  anomaly. No gene symbols or accessions appear (leak check: 0 original ids, 0 UniProt-pattern
  hits). Mapping: `results/pillar4_v2_id_mapping.tsv`.
* **Labels** — the jar's GO label lookup returns Reactome cross-reference labels for some
  terms (e.g. GO:0004069 "GOT2 dimer transaminates L-Asp"; an `rdfs:label` on an
  `owl:Axiom` annotating a `hasDbXref`). All target labels and `GO:x(label@p)` annotation
  strings are therefore replaced by the `name:` of the same GO release (go.obo 2026-03-25):
  9/80 target labels and 890/4,077 annotation strings changed. The 20-label check in the leak
  report lists jar label vs OBO name.
* Targets absent from the guardrail's GO release (go-basic 2025-06-01) were excluded before
  matching (the guardrail fails closed on unknown ids): 7 PATHWAY_GAP and 1 conflict record.

Composition: 80 items = PATHWAY_GAP 20 TP + 20 decoy (4 bacteria), ONTOLOGY_CONFLICT 20
genuine (all 7 genomes) + 20 traps (5 genomes; no eukaryote trap with two claims in the
window existed for athaliana/celegans). Leak checks (`results/pillar4_v2_leak_checks.tsv`):
model-visible key paths and context field sets identical between positives and negatives in
both categories; PATHWAY_GAP anomaly-level AUCs TP-vs-decoy 0.47 (max posterior), 0.53 (max
suggestion score), 0.49 (fraction of candidates in the pathway operon). ONTOLOGY_CONFLICT:
posterior AUC 0.51, but every genuine item carries the integrator's own `consistency` prior
(-0.68) and no trap does (AUC 1.00): the integrator already flags these terms, so the
category measures whether a method uses that symbolic signal, not hidden-label inference.

## Scoring (v2, `evaluate_pillar4_v2.py`)

TP correct = PROMOTE of an acceptable held-out protein with the target term; decoy correct =
ABSTAIN/DEMOTE; genuine conflict correct = DEMOTE or RELAX to a guardrail-valid proper
ancestor (ABSTAIN counted wrong; the lenient variant is in `pillar4_v2_per_category.tsv`);
trap correct = ABSTAIN or PROMOTE of the same term on the same protein. Precision = TP
promotions / all PROMOTE actions on TP, decoy and genuine items; recall over TP items;
accuracy over all items; violations = final PROMOTE/RELAX whose term fails the guardrail.
Baselines: posterior threshold (≥0.35) + guardrail; LR / RF / GBDT / MLP on 3 pre-tool scalars
(posterior, max claim probability, suggestion score) and 8 post-tool features, candidate-level,
5-fold StratifiedKFold over anomalies (category × label strata, seed 42), each with and
without the guardrail rule; single features (posterior, final log-odds, suggestion score),
argmax candidate with a threshold tuned on the training folds. Post-processing never reads
labels. McNemar: exact two-sided binomial on per-item correctness, c3 seed 42 vs every method.
LRO-50: `results/pillar4_lro_50_cases.jsonl` unchanged (sha256 `059d3f02…336`), protein ids
de-identified by the harness for the model and mapped back for scoring; single-feature rows
from `evaluate_pillar4.py` (whose `--lro-stats-only` path is unchanged).

## Trigger rate (88 proteomes)

`trigger_rates_v2.py` re-runs every proteome of the 88-proteome study with its own
`iter3_relax` manifest row (batch configuration: no dark-matter suggester) and the extractor,
no LLM. `n_triggered` = distinct proteins referenced by the extractor's (capped) anomalies.

## Run

See `run.sh` (steps `stage`, `baseline`, `select`, `holdouts`, `plan`, `decoys`, `assemble`,
`agent`, `evaluate`, `trigger`). Round scripts actually executed on 2026-10-19:
`r1_ws.sh` (r1tp, r1tpw, Stage-2 conflict pass), `r1b_ws.sh` (r1tph, Stage-1 conflicts),
`r2a_ws.sh` (r2tpw), `r2b_ws.sh` (decoys of the superseded tpw plan), `r3_ws.sh` (r3tp),
`r4_ws.sh` (decoys of the tp plan), `r5_solo_ws.sh` (check: 2 TP + 1 decoy benchmark cases
re-run alone; candidates and context identical to their packed runs). Include lists:
`make_includes.py`.

## Jobs (2026-10-19)

| Job | Host | What |
| --- | --- | --- |
| 9239 | unimatrix01 node004 | pilot (ecoli, maxIter 6, taxon 83333; not used) |
| 9264 | node004 | eukaryote baselines |
| 9263 | node007 | failed silently: node007 has no `/data` |
| 9277 | node007 (`/scratch/hohndor_p10` mirror via `node007.sh`) | prokaryote baselines |
| 9296 | — | cancelled (node007 fully allocated) |
| ws r1, r1b, r2a, r2b, r3, r4 | ws | hold-out and conflict runs (logs in `logs/`) |
| 9308 | node003 | 88-proteome trigger runs (88/88 rc 0) |
| ws r5 | ws | batching check (identical records) |

## Results (2026-10-19)

Benchmark v2 (80 items; c3 = guarded multi-turn agent):

| Method | P | R | F1 | Acc | Violations |
| --- | --- | --- | --- | --- | --- |
| Deterministic Stage 1+2 (abstain) | 0 | 0 | 0 | 0.500 | 0 |
| c3 seed 42 (T=0) | 0.706 | 0.600 | 0.649 | 0.775 | 0 |
| c3 3 seeds, mean ± SD | 0.618 ± 0.087 | 0.550 ± 0.041 | 0.581 ± 0.061 | 0.717 ± 0.046 | 0 |
| c2 fast (2 turns) | 1.000 | 0.150 | 0.261 | 0.625 | 0 |
| c4 unguarded | 0.545 | 0.600 | 0.571 | 0.600 | 1 |
| c0 zero-shot | 0 | 0 | 0 | 0.275 | 1 |
| Posterior ≥ 0.35 + guardrail | 0 | 0 | 0 | 0.500 | 0 |
| Best supervised (LR/RF/MLP, either feature set) + guardrail | 0 | 0 | 0 | 0.750 | 0 |
| Same without guardrail | 0 | 0 | 0 | 0.500 | 0–7 |

PATHWAY_GAP true positives (correct held-out protein chosen, 20 items): c3 12/11/10,
c4 12, c2 3, c0 0; argmax posterior 2%, final log-odds 13%, suggestion score 15%, uniform
random 12.5%, operon-context heuristic 44% (`pillar4_v2_pg_tp_hit1.tsv`). Per category (c3
seed 42): TP 12/20, decoys 16/20, genuine conflicts 20/20, traps 14/20. Exact McNemar (c3 seed
42 vs): deterministic p = 9.4e-4, c2 p = 0.029, c4 p = 0.0066, c0 p = 2.3e-8, threshold
p = 2.0e-4, supervised + guardrail p = 0.42–0.83, supervised without guardrail p ≤ 9.4e-4.

LRO-50 (unchanged cases, de-identified prompts): RLGC order Hit@1 34%, MRR 0.570; c3 Hit@1
84 / 86 / 88% (MRR 0.892 / 0.899 / 0.913); c2 36%, c0 40%, c4 94% (MRR 0.956); single
features suggestion score 34%, final log-odds 72%, ESM-2 centroid 96% (MRR 0.977); supervised
84–96%. McNemar c3 seed 42 vs ESM-2 centroid p = 0.070 (1 vs 7 discordant).

Trigger rate (88 proteomes, 111,253 proteins with claims, production cap 150): 8,893
anomalies referencing 5,975 proteins = 5.37% pooled (median 5.61%; bacteria 5.18%, archaea
6.75%, eukaryotes 6.53%, viruses 4.95%); the cap is reached in 52/88 proteomes. Uncapped:
31,839 anomalies, 16,388 proteins (14.7%); PATHWAY_GAP 4,294, INCOMPLETE_COMPLEX 6,
BORDERLINE_ESSENTIAL 83, BORDERLINE_CONFLICT 27,456. The previously quoted 1.84% is not
reproduced by any definition computed here.

Latency is wall-clock per anomaly on a shared gateway with three configurations running
concurrently (c3 ≈ 264–273 s, c4 215 s, c2 165 s, c0 86 s) and is not comparable to v1.
