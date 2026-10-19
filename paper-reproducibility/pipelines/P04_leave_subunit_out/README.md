# P04 — Leave-subunit-out (LSO) heteromeric complex benchmark

Manuscript elements: Table 3 (`tab:complex_lso`), Fig. 4a–b.
Compute host: unimatrix01 (SLURM, partition `debug`, node007 excluded; never the login node).
Entry point: `benchmark/unified/pillar3/run_pillar3_complex_lso.py` (regenerated 2026-10-19),
driven by `run.sh`. Stage 3 runs through `p04_agent_wrapper.py` (this folder).

## How to run

```bash
bash run.sh deploy                                   # local checkout -> unimatrix01:$PIPE
ssh unimatrix01 'bash /data/hohndor/gspa/proteomes/pillar3_unified/repro_2026-10-19/P04_leave_subunit_out/run.sh submit'
bash run.sh fetch                                    # results/pillar3_complex_lso{,_mcnemar}.tsv
```

The job runs four resumable steps (`prepare`, `integrate`, `agent`, `evaluate`; outputs that
already exist are skipped). Any single step: `sbatch run.sh job <step>`.

## Inputs (all on unimatrix01)

| Input | Path |
| --- | --- |
| Frozen unrelated_loci cases (Oct 15 23:39) | `/data/hohndor/gspa/proteomes/pillar3_unified/{g}_lso_cases.json` |
| Frozen unrelated_loci LSO claims (Oct 15) | `/data/hohndor/gspa/proteomes/pillar3_unified/{g}_claims_lso.jsonl` |
| Source claims (DIAMOND + InterProScan + Pfam) | `/data/hohndor/gspa/proteomes/bench9/with_iprscan/{g}_claims_merged.jsonl` |
| Truth (GOA, dual) | `/data/hohndor/gspa/proteomes/truth_dual/{g}_truth_all.tsv` |
| Operons (UniProt space, used by the jar) | `/data/hohndor/gspa/proteomes/operons/{g}_operons.tsv` |
| Operons (RefSeq space, used for within_operon decoys) | `/data/hohndor/gspa/proteomes/operons/{g}_refseq_operons.tsv` |
| Gene order | `/data/hohndor/gspa/proteomes/{g}_genomic.gff` (RefSeq CDS) |
| RefSeq -> UniProt | `/data/hohndor/gspa/proteomes/{g}.refseq_to_uniprot.tsv`, `.../refseq_uniprot_maps/{g}.refseq_to_uniprot.tsv`, GFF `Dbxref` UniProtKB/GOA |
| GO | `/data/hohndor/gspa/reference/go.obo` |
| Bundled complex classes | `/data/hohndor/gspa/proteomes/pillar3_unified/resources/complex-terms/complex_terms.tsv` |
| Guardrail tables (taxon, anatomy) | `/data/hohndor/gspa/gspa-core/src/main/resources/{taxon-constraints,anatomical-constraints}` |
| Jar | `/data/hohndor/gspa/proteomes/pillar3_unified/gspa-pillar3.jar` (GSPA v2.1.0, built Oct 15; the jar that produced the Oct-15 `{g}_m*.tsv`) |
| Java | `/data/hohndor/envs/jdk-21.0.10+7/bin/java` |
| Stage 3 sidecar | `/data/hohndor/gspa/benchmark-unified/agent/run_gspa_agent.py` (identical to `benchmark/agent/run_gspa_agent.py`) |
| LLM | vLLM at `http://unimatrix01.kaust.edu.sa:8000/v1`, served name `qwen3.8-27b` (Qwen3.8-27B-FP8), temperature 0, seed 42, 5 turns; API key read by the sidecar from `$BORG_LLM_API_KEY` or `~/.config/borg-llm/api-key` |

`{g}` ∈ ecoli (calibration genome), bsubtilis, mtb, paeruginosa, synechocystis, ecolo157.
Out-of-sample cohort = the five non-calibration genomes.

## Seeds and frozen parameters

* unrelated_loci: `random.seed(42)` once, genomes processed in the order
  ecoli, bsubtilis, mtb, paeruginosa, synechocystis, ecolo157 (one RNG stream). 30 complexes per
  genome (25 for Synechocystis), held-out raw 0.50 (2-subunit complex) / 0.45, companions raised to
  ≥ 0.85, 5 decoys with U(0.42, 0.49) rounded to 3 decimals. The `prepare` step re-runs the original
  generator in memory and records in `regen_check.tsv` whether it reproduces the frozen case files.
* within_operon: `random.Random(42*1000 + genome_index)` (index in the order above), used only to
  sample 5 of ≥ 5 eligible operon members.
* Stage 3: seed 42, temperature 0, `--max-turns 5`, `--max-concurrency 8`.

## Decoy designs

* `unrelated_loci` — the original design, unchanged (frozen case/claim files reused). Decoys are
  random proteins without any truth annotation, i.e. from loci unrelated to the complex.
* `within_operon` — same complexes, held-out subunits, companion boosts and the **same five
  synthetic decoy score values per case** (so the decoy score distribution is identical across
  designs). Decoys are drawn from the held-out subunit's own operon (RefSeq operon file, mapped to
  UniProt) when it has ≥ 5 eligible members; otherwise all eligible operon members plus the nearest
  genes of the held-out subunit's co-directional run (directon: maximal run of adjacent same-strand
  CDS); if still < 5, the nearest CDS on the same contig regardless of strand ("fill"). Eligible =
  UniProt-mapped and not annotated to the complex term in the truth set. When the held-out subunit
  has no chromosome position (2 cases in ecolo157), the first placed companion subunit is used as
  the anchor. Per-genome counts: `within_operon_decoy_stats.tsv`; per-case: `n_decoys_operon`,
  `n_decoys_directon`, `n_decoys_fill`, `decoy_anchor_mode` in `within_operon/{g}_lso_cases.json`.
  The claims file is the frozen unrelated_loci file with the 5 synthetic unrelated decoy claims per
  case removed and the 5 within-operon decoy claims appended; real claims a within-operon decoy
  already has for the complex term are kept (count in `pre_existing_claims_on_new_decoys`).

## Methods

| method_id | What it ranks by |
| --- | --- |
| `input_attenuated` | max raw score over all claims for (protein, complex) in the LSO claims file; no calibration, no priors |
| `real_diamond_channel` | real DIAMOND claims (`source=diamond`) for the complex term from `{g}_claims_merged.jsonl`; the held-out subunit's real claim, and every decoy's, is capped at the case's attenuation level (0.50 / 0.45), i.e. the same attenuation as the LSO input; synthetic decoy claims are not used |
| `real_interproscan_channel` | same with `source ∈ {interproscan, pfam}` |
| `operon_only` | jar `--enable-priors=genomic_context` |
| `noisy_or` | jar `--enable-priors=none` (Stage 0 pooling) |
| `coherence` | jar `--enable-priors=coherence` |
| `coherence_operon` | jar `--enable-priors=coherence,genomic_context` |
| `coherence_operon_stage2` | + `--enforce-stage2` |
| `full_stage3` | Stage 2 output with Stage 3 resolutions applied (PROMOTE → max(0.58, confidence), DEMOTE → 0.35) |

Metrics per case among {held-out, 5 decoys}: Hit@1 (ties count 1/ties, as originally),
`hit_at_1_strict` (held-out uniquely top), recall/precision/F1 at 0.50; complex coherence is
genome-wide (fraction of heteromeric complex terms with ≥ 1 protein at ≥ 0.50 that have ≥ 2).
Cohort rows are case-weighted means; coherence is the unweighted mean over genomes.
McNemar: exact two-sided binomial test on discordant pairs of the strict Hit@1 indicator
(`scipy.stats.binomtest`).

## Audit of the Oct-17 script (fixed or documented here)

1. ~80 % of complexes are selected to be co-operonic and the held-out subunit is chosen to share an
   operon with a companion — kept for `unrelated_loci`; `within_operon` is the control.
2. Companions are forced to raw ≥ 0.85 — kept in both designs.
3. The rows labelled "DIAMOND-BLASTp" and "InterProScan / Pfam" ranked the synthetic LSO input
   (the 0.45/0.50 held-out claim and U(0.42,0.49) decoy claims all carry `source=diamond`); the
   InterProScan row also matched only `interpro|pfam|iprscan` and so missed `source=interproscan`.
   Replaced by `input_attenuated` and the two `real_*_channel` rows.
4. Hit@1 ties scored 1/ties — kept; strict Hit@1 added.
5. Stage 3 failed closed: the sidecar looks for `benchmark/neural/cafa6_recon/go-basic.obo` two
   directories above itself, which does not exist on unimatrix01, so `SymbolicGuardrail` reported
   missing resources and rejected every PROMOTE. Fixed by `p04_agent_wrapper.py`, which points the
   unmodified sidecar at `reference/go.obo` and aborts if the guardrail does not load.
6. The script's `load_operons` treats the first member of each operon line as an id and does not
   map it (the jar parses all members). Kept unchanged for comparability (it chose the frozen
   held-out subunits and feeds Stage 3 anomaly features); within_operon decoys use a correct parse.
7. The "p < 1e-15, exact McNemar" statement had no code behind it; now computed
   (`pillar3_complex_lso_mcnemar.tsv`).
8. Operon-vote Hit@1 0.93 (results/pillar3_report.md) vs 0.58: in commit 3310b97 (Oct 17) the
   operon row was a hand-written Python heuristic (candidates co-operonic with a ≥ 0.70 anchor got
   +0.12, all others ×0.90), which by construction favours the co-operonic held-out subunit. Commit
   c147097 (Oct 19) replaced it by the jar's `genomic_context` prior alone (no coherence prior, no
   complex-terms file), giving 0.58. The same 3310b97 commit also contained a hand-coded
   promote/demote rule standing in for Stage 3, removed in c147097.

## Reproduction checks (run of 2026-10-19, job 9245 on node004)

* `regen_check.tsv`: the original generator with seed 42 reproduces all six frozen case files
  and both complex term/module files exactly.
* All 30 unrelated_loci jar outputs (5 methods × 6 genomes) are byte-identical to the Oct-15/17
  files in `pillar3_unified/`.
* `--enforce-stage2` produces output byte-identical to Stage 1 coherence + operon in every genome
  and both designs (no Stage 2A change on this benchmark), so `coherence_operon_stage2` equals
  `coherence_operon` and Stage 3 is applied on top of an unchanged Stage 1 matrix.
* `sha_inputs.sh` records path/size/sha256 of every input in `$REPRO/P04_inputs.sha256.tsv`
  (also appended to `manifest/data_manifest.tsv`, `used_by_pipeline = P04`).
