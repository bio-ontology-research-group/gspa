# Price development regression: full target candidates

At engine revision `e2303a9`, the current pathway-prior configuration produces
**no EC top-1 gain** on these inputs. All seven configurations executed on all
22 genomes. This is development reuse of Price data examined in previous audits,
with limited genomic context; it is not a blind v3 test or the completed
all-module benchmark.

## Matched comparisons

Each configuration uses identical claims in stage0 (priors disabled) and stage1
(default priors, with the same EC-to-GO/pathway/module reference files). The
candidate export retains all nonzero scores for the 148 target proteins. One
of the original 149 entries is excluded because its mapped identifier is a
replicon rather than a scored protein. There are 151 complete-EC labels. The
[input protocol](v3-price148.md) gives all identifiers, exclusions and limitations.

| Exported predictor evidence | Stage0 correct / 148 | Stage1 correct / 148 | Accuracy | All-case report |
| --- | ---: | ---: | ---: | --- |
| CLEAN | 60 | 60 | 40.5% | [JSON](../../benchmark/v3/results/development-e2303a9/price-clean-report.json) |
| DeepProZyme | 49 | 49 | 33.1% | [JSON](../../benchmark/v3/results/development-e2303a9/price-dpz-report.json) |
| EnzBERT | 53 | 53 | 35.8% | [JSON](../../benchmark/v3/results/development-e2303a9/price-enzbert-report.json) |
| MAPred | 50 | 50 | 33.8% | [JSON](../../benchmark/v3/results/development-e2303a9/price-mapred-report.json) |
| GraphEC | 18 | 18 | 12.2% | [JSON](../../benchmark/v3/results/development-e2303a9/price-graphec-report.json) |
| TopEC | 5 | 5 | 3.4% | [JSON](../../benchmark/v3/results/development-e2303a9/price-topec-report.json) |
| CLEAN + DeepProZyme + EnzBERT + MAPred | 60 | 60 | 40.5% | [JSON](../../benchmark/v3/results/development-e2303a9/price-ensemble_pool4-report.json) |

Every matched stage1 comparison has zero top-1 corrections and zero harms.
Reports include genome-level bootstrap intervals (2,000 resamples), fixed-score
metrics, descriptive Fmax, calibration diagnostics, per-genome metrics, execution
records, commands and hashes. Fixed-threshold values do not establish usable
confidence: the current default calibration is provisional and not fitted on
independent calibration data. Calibration summaries across different predictor
configurations use different candidate universes and should not be ranked as if
they were directly comparable.

These numbers describe the **acquired predictor outputs**, not a rerun of each
predictor's published benchmark protocol. Original endpoint training overlap,
reference databases, and some architecture/evidence-type declarations still need
verification. The single-predictor raw rankings give the same target top-1
accuracies as stage0, checked separately below.

## Why this comparison is limited

The export includes target proteins and truncated nearby evidence but has no
gene coordinates, complete proteomes, operon calls, gap input, or organism model.
Reaction-local, cohort, FBA and agent stages were not enabled. Much of the current
pathway/essentiality/context code acts on GO candidates; these inputs are EC
candidates. This result therefore exposes the limited connection between the
current pathway priors and direct EC predictions. It does not test the new
reaction/cohort/learned components, which still need production wiring, real
training and held-out evaluation. It also cannot measure enzyme detection:
every labelled target is already a known enzyme.

## Full-Proteome 7-Variant × 6-Arm Context, Cohort, Learned-Model, and Locus-Decoy Evaluation

After acquiring all 22 exact-version NCBI assemblies (`115,542` proteins, `115,263` unambiguous CDS coordinates, `279` explicitly declared `unplacedProteins` in `cohort_metadata.json`) and building the 35,298 MMseqs2 30%-identity/80%-coverage sequence clusters (`familySearchComplete: true`, 18 NCBI genera/species in [`docs/benchmarks/v3-price-context-data.md`](v3-price-context-data.md)), all 7 predictor variants (`clean`, `dpz`, `enzbert`, `graphec`, `mapred`, `topec`, `ensemble_pool4`) were executed through the production `gspa integrate` CLI across 6 matched inference arms (`benchmark/v3/results/price-context-ablation/`, `summary.json` SHA-256 `a6bb077ac048d4d79ff13a9ad859687e4728acf4dfe5584412e664b7b84396aa`, `locus_decoy_evaluation.json` SHA-256 `5fea2947077f245f4a7bd4f7ccf3bfd28a6d5879b59b0158bacc08c8be84ead6`):

1. `stage0`: sequence-only integration (`--enable-priors false`)
2. `stage1`: default pathway/module/operon priors (`--enable-priors true`)
3. `stage1_rlc`: Stage 1 + Reaction-Local Context (`--rlc-suggester`)
4. `stage1_rlc_cohort`: Stage 1 + RLC + clade-held-out Cross-Genome Cohort Recurrence (`--cohort --cohort-holdout CLADE`)
5. `learned_logistic`: Stage 1 + RLC + Cohort + `model_logistic.json` (trained on the genus/family-disjoint `swissprot_label_rich_v1` cohort)
6. `learned_gbdt`: Stage 1 + RLC + Cohort + `model_gbdt.json` (trained on `swissprot_label_rich_v1`)

| Variant | Arm | Protein EC Top-1 | Protein EC Top-3 | Descriptive $F_{\max}$ | Brier | 10-Bin ECE | Same-EC Decoy Locus Top-1 | All-148 Locus MRR | Report |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| `clean` | `stage0` / `stage1` / `stage1_rlc` / `stage1_rlc_cohort` | **0.4054** (60/148) | **0.5811** (86/148) | **0.3862** | 0.0061 | 0.0774 | **0.5938** (38/64) | **0.8998** | [JSON](../../benchmark/v3/results/price-context-ablation/price-clean-report.json) |
| `clean` | `learned_logistic` | 0.0135 (2/148) | 0.0135 (2/148) | 0.0219 | **0.0017** | **0.0373** | 0.1250 (8/64) | 0.7939 | [JSON](../../benchmark/v3/results/price-context-ablation/price-clean-report.json) |
| `clean` | `learned_gbdt` | 0.0000 (0/148) | 0.0068 (1/148) | 0.0007 | 0.0027 | 0.0447 | 0.2344 (15/64) | 0.8131 | [JSON](../../benchmark/v3/results/price-context-ablation/price-clean-report.json) |
| `dpz` | `stage0` / `stage1` / `stage1_rlc` / `stage1_rlc_cohort` | **0.3311** (49/148) | **0.6014** (89/148) | **0.3696** | 0.0063 | 0.0775 | **0.6667** (40/60) | **0.7917** | [JSON](../../benchmark/v3/results/price-context-ablation/price-dpz-report.json) |
| `dpz` | `learned_logistic` | 0.0608 (9/148) | 0.0878 (13/148) | 0.0693 | **0.0018** | **0.0362** | 0.2500 (15/60) | 0.6948 | [JSON](../../benchmark/v3/results/price-context-ablation/price-dpz-report.json) |
| `dpz` | `learned_gbdt` | 0.0000 (0/148) | 0.0000 (0/148) | 0.0011 | 0.0029 | 0.0443 | 0.4167 (25/60) | 0.7314 | [JSON](../../benchmark/v3/results/price-context-ablation/price-dpz-report.json) |
| `enzbert` | `stage0` / `stage1` / `stage1_rlc` / `stage1_rlc_cohort` | **0.3581** (53/148) | **0.5068** (75/148) | **0.4403** | 0.0139 | 0.0715 | **0.4167** (5/12) | **0.6216** | [JSON](../../benchmark/v3/results/price-context-ablation/price-enzbert-report.json) |
| `enzbert` | `learned_logistic` | 0.2027 (30/148) | 0.2095 (31/148) | 0.2488 | **0.0104** | **0.0076** | **0.4167** (5/12) | **0.6216** | [JSON](../../benchmark/v3/results/price-context-ablation/price-enzbert-report.json) |
| `enzbert` | `learned_gbdt` | 0.0000 (0/148) | 0.0000 (0/148) | 0.0138 | 0.0118 | 0.0299 | **0.4167** (5/12) | **0.6216** | [JSON](../../benchmark/v3/results/price-context-ablation/price-enzbert-report.json) |
| `graphec` | `stage0` / `stage1` / `stage1_rlc` / `stage1_rlc_cohort` | **0.1216** (18/148) | **0.2500** (37/148) | **0.2095** | 0.0089 | 0.0934 | **0.6562** (42/64) | **0.9223** | [JSON](../../benchmark/v3/results/price-context-ablation/price-graphec-report.json) |
| `graphec` | `learned_logistic` | 0.0068 (1/148) | 0.0473 (7/148) | 0.0305 | **0.0018** | **0.0388** | 0.2500 (16/64) | 0.8240 | [JSON](../../benchmark/v3/results/price-context-ablation/price-graphec-report.json) |
| `graphec` | `learned_gbdt` | 0.0000 (0/148) | 0.0068 (1/148) | 0.0006 | 0.0027 | 0.0445 | 0.3594 (23/64) | 0.8446 | [JSON](../../benchmark/v3/results/price-context-ablation/price-graphec-report.json) |
| `mapred` | `stage0` / `stage1` / `stage1_rlc` / `stage1_rlc_cohort` | **0.3378** (50/148) | **0.6351** (94/148) | **0.3551** | 0.0528 | 0.1134 | **0.5909** (13/22) | **0.7387** | [JSON](../../benchmark/v3/results/price-context-ablation/price-mapred-report.json) |
| `mapred` | `learned_logistic` | 0.1081 (16/148) | 0.3176 (47/148) | 0.1517 | 0.0408 | 0.0215 | 0.2727 (6/22) | 0.7128 | [JSON](../../benchmark/v3/results/price-context-ablation/price-mapred-report.json) |
| `mapred` | `learned_gbdt` | 0.0811 (12/148) | 0.1419 (21/148) | 0.0734 | **0.0377** | **0.0066** | 0.2273 (5/22) | 0.7083 | [JSON](../../benchmark/v3/results/price-context-ablation/price-mapred-report.json) |
| `topec` | `stage0` / `stage1` / `stage1_rlc` / `stage1_rlc_cohort` | **0.0338** (5/148) | **0.0338** (5/148) | 0.0344 | 0.6282 | 0.7230 | 0.0000 (0/1) | **0.0338** | [JSON](../../benchmark/v3/results/price-context-ablation/price-topec-report.json) |
| `topec` | `learned_logistic` | **0.0338** (5/148) | **0.0338** (5/148) | 0.0369 | 0.6038 | 0.6675 | 0.0000 (0/1) | **0.0338** | [JSON](../../benchmark/v3/results/price-context-ablation/price-topec-report.json) |
| `topec` | `learned_gbdt` | **0.0338** (5/148) | **0.0338** (5/148) | **0.0380** | **0.5138** | **0.4903** | 0.0000 (0/1) | **0.0338** | [JSON](../../benchmark/v3/results/price-context-ablation/price-topec-report.json) |
| `ensemble_pool4` | `stage0` / `stage1` / `stage1_rlc` / `stage1_rlc_cohort` | **0.4054** (60/148) | **0.7095** (105/148) | **0.4419** | 0.0063 | 0.0778 | **0.7121** (47/66) | **0.9324** | [JSON](../../benchmark/v3/results/price-context-ablation/price-ensemble_pool4-report.json) |
| `ensemble_pool4` | `learned_logistic` | 0.0541 (8/148) | 0.1824 (27/148) | 0.0903 | **0.0017** | **0.0372** | 0.2727 (18/66) | 0.8226 | [JSON](../../benchmark/v3/results/price-context-ablation/price-ensemble_pool4-report.json) |
| `ensemble_pool4` | `learned_gbdt` | 0.0000 (0/148) | 0.0000 (0/148) | 0.0010 | 0.0027 | 0.0447 | 0.3030 (20/66) | 0.8212 | [JSON](../../benchmark/v3/results/price-context-ablation/price-ensemble_pool4-report.json) |

### Mechanistic Findings from the Full-Proteome Price Ablation

1. **Candidate retention is 100% of predictor-emitted truth (`candidate_recall = 1.0000` on `clean`, `graphec`, `ensemble_pool4`)**: No candidate pruning removes true functions in `stage0`, `stage1`, `stage1_rlc`, `stage1_rlc_cohort`, `learned_logistic`, or `learned_gbdt`.
2. **Ensembling 4 predictors (`ensemble_pool4`) substantially improves top-3 EC recall and locus decoy discrimination over any single predictor**: While top-1 protein EC accuracy matches `clean` (`0.4054` = `60/148`), top-3 EC recall rises from `0.5811` (`86/148` in `clean`) and `0.6351` (`94/148` in `mapred`) to **`0.7095`** (`105/148`), descriptive $F_{\max}$ reaches **`0.4419`**, and same-EC chromosomal decoy locus top-1 rises to **`0.7121`** (`47/66`, MRR across all 148 loci = **`0.9324`**).
3. **Why `stage1_rlc` and `stage1_rlc_cohort` match `stage0` in default sidecar mode**: By design (`IntegrateCommand`), `--rlc-suggester` and `--cohort` compute structured reaction-local and clade-held-out cohort evidence and emit them to the `.cohort.json` and `--feature-export-dir` sidecars without overriding the Stage 1 score unless `--learned-model` is supplied.
4. **Historical transfer reduces descriptive Brier/ECE while degrading within-protein EC ranking; useful biological calibration is unproven**:
   - **Reference-label metrics**: The historical learned arms reduce Brier/ECE relative to stage0 (e.g., CLEAN Brier falls from `0.0061` to `0.0017`). However, an unfitted all-zero predictor has lower Brier than 13 of the 14 learned combinations, and offers no useful ranking. The large unlabeled-pair denominator makes these reductions insufficient evidence of accurate probabilities. See the [complete zero-baseline audit](v3-price-calibration-audit.md). Unlabeled pairs follow a closed-world convention and are not experimentally proven negatives.
   - **Within-protein EC ranking under covariate shift**: In `swissprot_label_rich_v1`, sequence claims come from DIAMOND homology (`seq.group.homology.max_raw` weight `+0.1739`), whereas `seq.group.ec_specialist.max_raw` is identically `0.0` in training (`weight = 0.0`). On Price, claims come from EC specialist neural networks (`clean`, `dpz`, `enzbert`, `graphec`, `mapred`, `topec`), so only `seq.max_raw_any` (`+0.1739`) fires. Furthermore, in dense full-vocabulary matrices (`clean`, `graphec`, `dpz` where every protein has up to 5,242 EC claims), every common EC is claimed by dozens of proteins in the same genome (`ctx.competitor.n_with_claim` weight `-0.2215`) and maps to multiple ModelSEED reactions (`cohort.w.n_reaction_alternatives` weight `-0.2560`). Consequently, the linear/tree combiner penalizes common ECs with many reaction mappings and competitors more than the `+0.1739` sequence weight rewards the top EC, causing rare ECs with zero ModelSEED reaction mappings (`n_reaction_alternatives = 0`) to outrank the true EC within a protein. On sparse predictors (`topec`, `enzbert`), where each protein has only 1–5 EC claims, `learned_logistic` retains much higher top-1 accuracy (`0.2027` on `enzbert`, `0.0338` on `topec`).

## Raw ranking and numerical precision

The [raw ranking audit](../../benchmark/v3/results/development-e2303a9/price-raw-rank-audit.json)
compares each single predictor's original scores with stage0 output and verifies
claims, labels and predictions against the recorded run hashes. All six preserve
the measured target top-1 outcomes. Five preserve the complete exported ranking.
CLEAN has 685,560 adjacent raw-distinct pairs that become tied after the current
sigmoid/log-odds transformations; none changes a unique raw top-1. Full-precision
TSV output fixed decimal-rounding loss, but cannot reverse differences already
lost to floating-point arithmetic. Consequently stage0 is not an exact raw-score
ranking baseline for every CLEAN candidate. The independently scored raw
baseline is retained in the audit rather than inferred from the transformed
posterior.

```bash
python benchmark/v3/audit_price_runs.py --inputs /path/to/price-inputs \
  --run-root /path/to/first-attempt/results \
  --run-root /path/to/retry/results --out /path/to/new-raw-audit.json
```

Slurm array 9717 produced DeepProZyme, TopEC and pool4 results. Other tasks failed
or were cancelled before producing results; node002 retries 9728 and 9731
produced EnzBERT, MAPred, GraphEC and CLEAN results in a separate directory.
Verified ranking audit job 9735 completed afterward. Failed attempts were not
counted as successful execution and their directories were retained.
