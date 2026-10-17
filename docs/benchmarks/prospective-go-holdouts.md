# Prospective GO holdouts

## Question

Does pooling complementary per-protein predictors in Stage 0 give better Gene Ontology predictions than single predictors for proteins that received experimental annotations only after the training cutoff, and does genome-scale refinement change accuracy in the low-identity "twilight zone"?

Paper: Supplementary Tables S4 and S15 (Suites 1A and 1B in Supplementary Table S1).

## Datasets

### CAFA temporal holdout

Proteins that acquired experimental GO annotations between the CAFA3 and CAFA5 snapshots. Each (protein, aspect) pair with new annotations is one evaluation. The archived evaluation input `benchmark/manifests/reference_data/ibex_extracted/cafa_temporal_holdout_real.json` records the cohort per aspect:

| Aspect | Evaluations | \(<30\%\) identity | 30-50% | 50-70% | 70-100% |
|---|---|---|---|---|---|
| MF | 4,712 | 901 | 224 | 479 | 3,108 |
| BP | 5,444 | 1,085 | 137 | 332 | 3,890 |
| CC | 5,969 | 1,454 | 302 | 458 | 3,755 |

In total there are 16,125 protein-aspect evaluations, of which 3,440 share less than 30% sequence identity with any training sequence. The temporal-holdout proteins and their ground-truth terms are frozen in `benchmark/manifests/suite1_cafa_temporal_manifest.tsv`.

### Prospective CAFA6 targets

769 CAFA6 targets (412 no-knowledge, 357 limited-knowledge; 298 MF, 415 BP, 384 CC targets), evaluated as isolated sequences. Because these targets come without source genome coordinates, genome-scale Stage 1 priors cannot act; Stages 0, 2A, and 3 are active.

## Protocol

- **CAFA temporal holdout.** Stage 0 pools five channels (DIAMOND, Pfam, InterProScan, ESM-2, ProtT5), followed by Stage 1 and Stage 2A for 1, 2, or 3 iterations. Protein-centric \(F_{\max}\) is computed per aspect on ancestor-propagated annotations, for the full cohort and for each identity stratum. Macro \(F_{\max}\) is the mean of the MF, BP, and CC values. See [Metrics](index.md#metrics).
- **CAFA6.** Seven channels (including ESM-2, ProstT5 3Di, DIAMOND, InterProScan, DeepGO-SE) feed Stage 0, then Stage 2A and the Stage 3 agent. Reported metrics are per-aspect \(F_{\max}\) and information-accretion-weighted \(WF_{\max}\), their macro means, and macro \(F_{\max}\) separately for no-knowledge (NK) and limited-knowledge (LK) targets.

## Baselines

- **CAFA temporal holdout:** DIAMOND (BLAST), DeepGOCNN, DeepGOPlus, DeepGOZero, InterPro MLP, DeepGOZero+BLAST, ProteInfer.
- **CAFA6:** an MLP on ESM-2 and ProstT5 3Di embeddings, and an unguarded agentic baseline (Qwen3-14B with domain and EC tools).

## Reproduce

```bash
uv run --with numpy --with scipy --with pandas \
    python scripts/eval_scenario1_holdouts_real.py
```

This one script writes four tables: the CAFA temporal holdout, the Holdout-125 and 95-genome enzyme tables (see [Enzyme EC holdouts](enzyme-ec-holdouts.md)), and the CAFA6 summary. Note:

- For the CAFA temporal holdout the script re-tabulates per-method, per-aspect \(F_{\max}\) values that are already stored in `cafa_temporal_holdout_real.json`; it does not rerun the predictors. The `runtime_sec` and `peak_ram_mb` columns are constants written in the script.
- The CAFA6 rows are written as literal values inside the script (function `run_part_1d_cafa6_prospective`). `results/PROVENANCE_MAP.md` lists `results/suite1_cafa6_prospective_summary.tsv` as transcribed from cluster archives. The CAFA6 predictions and ground truth are not in the repository.
- The script sets `ROOT = "/home/leechuck/Public/software/gspa"` and writes into `results/`; see the warning in [How to reproduce](index.md#how-to-reproduce).

## Result tables

- `results/scenario1_cafa_temporal_holdout_real.tsv`: one row per method; columns `mf_all`, `mf_lt30`, `mf_30_50`, `mf_50_70`, `mf_70_100`, the same for `bp_*` and `cc_*`, `macro_fmax_all`, `macro_fmax_lt30`. ([GitHub](https://github.com/bio-ontology-research-group/gspa/blob/main/results/scenario1_cafa_temporal_holdout_real.tsv))
- `results/suite1_cafa6_prospective_summary.tsv`: per-aspect `*_fmax` and `*_wfmax`, `macro_fmax`, `macro_wfmax`, `nk_macro_fmax`, `lk_macro_fmax`. ([GitHub](https://github.com/bio-ontology-research-group/gspa/blob/main/results/suite1_cafa6_prospective_summary.tsv))

## Main finding

Stage 0 pooling reaches `macro_fmax_all` 0.6929 (row `Noisy-OR (Ensemble)`), above the strongest single predictor, ProteInfer (0.6138). Genome-scale refinement leaves full-cohort accuracy essentially unchanged and raises `macro_fmax_lt30` from 0.6723 to 0.6775 (row `GSPA (Iter 1)`); further iterations do not add to this. Most of the protein-level accuracy therefore comes from pooling; refinement contributes a small gain in the twilight zone.

On CAFA6, the guarded Stage 3 configuration has the highest macro \(F_{\max}\) of the three methods in the table. See [Guarded agentic completion](agentic-completion.md) for the agent itself.
