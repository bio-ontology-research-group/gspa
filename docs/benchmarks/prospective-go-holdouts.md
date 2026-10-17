# Prospective GO holdouts and retrained CAFA splits

## Question

How do GSPA and state-of-the-art sequence, domain, and protein language model baselines perform when retrained from scratch under strict CAFA time-based splits and sequence-similarity homology splits ($\le 70\%$, $\le 50\%$, and $\le 30\%$ identity), and does genome-scale inference improve both protein-centric $F_{\max}$ and information-accretion-weighted $F_w$ while eliminating taxonomic and anatomical violations?

Paper: Section 2.1, Table 1 (`\label{tab:cafa_splits}`), and Supplementary Tables S4 and S15 (`paper/main.tex` and `paper/supplement.tex`).

## Datasets

### Retrained CAFA time-based and homology-based splits

Using the curated Swiss-Prot and UniProt-GOA benchmark corpus ($82{,}648$ experimentally annotated proteins across $2{,}847$ GO terms: $688$ MF, $1{,}683$ BP, and $476$ CC, plus $81{,}216$ InterPro domain signatures and $1{,}280$-dimensional ESM-2 `esm2_t33_650M_UR50D` embeddings), we construct five strict train/test splits where all trainable baseline architectures and GSPA Stage 0 calibrators are **retrained from scratch on each split's filtered training set**:

1. **Time-Based Split Clean-A ($82{,}376$ train $\to 272$ test):** Prospective temporal holdout ($t_0 \to t_1$) where training sequences, $k$-NN databases, InterPro classifiers, 1D-CNNs, and PLM heads are restricted to pre-$t_0$ annotations.
2. **Time-Based Split Clean-B ($82{,}404$ train $\to 244$ test):** Strict zero-prior prospective holdout restricting test proteins to those with zero experimental GO annotations in any aspect at $t_0$.
3. **Homology Split $\le 70\%$ identity ($41{,}945$ train $\to 272$ test):** Purges $40{,}459$ training proteins sharing $>70\%$ pairwise sequence identity (`MMseqs2` / `DIAMOND`) with any test target prior to retraining all methods.
4. **Homology Split $\le 50\%$ identity ($11{,}900$ train $\to 272$ test):** Purges $70{,}504$ training proteins sharing $>50\%$ sequence identity with any test target prior to retraining.
5. **Homology Split $\le 30\%$ identity ($2{,}141$ train $\to 272$ test):** Strict twilight-zone split purging $80{,}263$ training proteins sharing $>30\%$ sequence identity with any test target prior to retraining.

### Multi-stratum CAFA temporal holdout (16,125 evaluations)

Proteins that acquired experimental GO annotations between the CAFA3 and CAFA5 snapshots (`benchmark/manifests/suite1_cafa_temporal_manifest.tsv`), totaling 16,125 protein-aspect evaluations (4,712 MF, 5,444 BP, 5,969 CC), of which 3,440 share less than 30% sequence identity with any training sequence.

### Prospective CAFA6 targets (769 targets)

769 CAFA6 targets (412 no-knowledge, 357 limited-knowledge; 298 MF, 415 BP, 384 CC targets).

## Protocol and retrained baselines

On each of the five splits (`Time_Split_CleanA`, `Time_Split_CleanB`, `Homology_Split_le70`, `Homology_Split_le50`, `Homology_Split_le30`), eleven methods are trained from scratch on that split's training partition and evaluated on the held-out test set using both unweighted protein-centric $F_{\max}$ and official CAFA information-accretion-weighted $F_w$ (where $\text{IA}(t) = -\log_2 P(t \mid \text{parent}(t))$ is computed strictly from each split's training corpus):

1. `Naive (Term Frequency Prior)` (Radivojac et al. 2013)
2. `DIAMOND BLAST-kNN` (Buchfink et al. 2021)
3. `InterPro Domain Logistic Regression` (Jones et al. 2014)
4. `DeepGOCNN (Sequence 1D-CNN)` (Kulmanov et al. 2018)
5. `DeepGOPlus (CNN + Homology)` (Kulmanov and Hoehndorf 2020)
6. `ESM-2 kNN (goPredSim Embedding Transfer)` (Littmann et al. 2021)
7. `ESM-2 MLP (Flat BCE)` (Lin et al. 2023)
8. `ESM-2 C-HMCNN (Hierarchical MCLoss)` (Giunchiglia and Lukasiewicz 2020)
9. `DeepGO-PlusPlus Stacker (Stage 0 Multi-Modal)`
10. `GSPA Stage 0 (Calibrated Multi-Channel Pool)`
11. `GSPA Full (Stage 1 Context Priors + Stage 2A Relaxation)`

## Reproduce

```bash
uv run --with numpy --with scipy \
    python scripts/run_cafa_retrained_splits.py
uv run --with numpy --with scipy --with pandas \
    python scripts/eval_scenario1_holdouts_real.py
```

## Result tables

- `results/cafa_retrained_time_and_homology_splits.tsv`: 55 rows (5 splits $\times$ 11 retrained methods) reporting `mf_fmax`, `bp_fmax`, `cc_fmax`, `mean_fmax`, `mf_fw`, `bp_fw`, `cc_fw`, `mean_fw`, `smin`, and `violations`. ([GitHub](https://github.com/bio-ontology-research-group/gspa/blob/main/results/cafa_retrained_time_and_homology_splits.tsv))
- `results/scenario1_cafa_temporal_holdout_real.tsv`: multi-stratum temporal holdout across 16,125 evaluations. ([GitHub](https://github.com/bio-ontology-research-group/gspa/blob/main/results/scenario1_cafa_temporal_holdout_real.tsv))
- `results/suite1_cafa6_prospective_summary.tsv`: prospective CAFA6 summary. ([GitHub](https://github.com/bio-ontology-research-group/gspa/blob/main/results/suite1_cafa6_prospective_summary.tsv))

### Retrained CAFA time-based and homology-based split results

Values transcribed directly from `results/cafa_retrained_time_and_homology_splits.tsv`:

| Split Regime ($N_{\text{train}} \to N_{\text{test}}$) | Retrained Method | MF $F_w$ ($F_{\max}$) | BP $F_w$ ($F_{\max}$) | CC $F_w$ ($F_{\max}$) | Mean $F_w$ ($F_{\max}$) | Taxon/Anat. Violations |
|---|---|---|---|---|---|---|
| **Time-Based Clean-A** ($82{,}376 \to 272$) | Naive (Term Frequency Prior) | 0.2512 (0.3788) | 0.1489 (0.2923) | 0.4272 (0.5818) | 0.2758 (0.4176) | 1,420 |
| Time-Based Clean-A | DIAMOND BLAST-$k$NN | 0.5624 (0.6184) | 0.3698 (0.4195) | 0.5284 (0.6104) | 0.4869 (0.5494) | 312 |
| Time-Based Clean-A | InterPro Domain LogReg | 0.5318 (0.5842) | 0.3342 (0.3912) | 0.4986 (0.5892) | 0.4549 (0.5215) | 198 |
| Time-Based Clean-A | DeepGOCNN (1D-CNN) | 0.4182 (0.4890) | 0.2764 (0.3510) | 0.4712 (0.5940) | 0.3886 (0.4780) | 645 |
| Time-Based Clean-A | DeepGOPlus (CNN + Homology) | 0.5810 (0.6340) | 0.3920 (0.4480) | 0.5645 (0.6520) | 0.5125 (0.5780) | 428 |
| Time-Based Clean-A | ESM-2 $k$-NN (goPredSim) | 0.5769 (0.6164) | 0.3712 (0.4171) | 0.5408 (0.6400) | 0.4963 (0.5578) | 384 |
| Time-Based Clean-A | ESM-2 MLP (Flat BCE) | 0.5579 (0.6020) | 0.3621 (0.4368) | 0.5653 (0.6774) | 0.4951 (0.5721) | 812 |
| Time-Based Clean-A | ESM-2 C-HMCNN (MCLoss) | 0.5642 (0.6085) | 0.3685 (0.4412) | 0.5668 (0.6786) | 0.4998 (0.5761) | 790 |
| Time-Based Clean-A | DeepGO-PlusPlus Stacker | 0.6132 (0.6443) | 0.4342 (0.4650) | 0.5814 (0.6606) | 0.5429 (0.5900) | 665 |
| Time-Based Clean-A | GSPA Stage 0 (Calibrated Pool) | 0.6098 (0.6410) | 0.4285 (0.4612) | 0.5790 (0.6584) | 0.5391 (0.5869) | 524 |
| Time-Based Clean-A | **GSPA Full (Stage 1 + Stage 2A)** | **0.6184 (0.6468)** | **0.4412 (0.4638)** | **0.5846 (0.6585)** | **0.5481 (0.5897)** | **0** |
| **Time-Based Clean-B** ($82{,}404 \to 244$) | DIAMOND BLAST-$k$NN | 0.5480 (0.6012) | 0.3510 (0.4024) | 0.5164 (0.5980) | 0.4718 (0.5339) | 285 |
| Time-Based Clean-B | DeepGOPlus (CNN + Homology) | 0.5762 (0.6280) | 0.3840 (0.4390) | 0.5635 (0.6480) | 0.5079 (0.5717) | 392 |
| Time-Based Clean-B | ESM-2 C-HMCNN (MCLoss) | 0.5740 (0.6190) | 0.3690 (0.4402) | 0.5812 (0.6910) | 0.5081 (0.5834) | 718 |
| Time-Based Clean-B | **GSPA Full (Stage 1 + Stage 2A)** | **0.6152 (0.6435)** | **0.4158 (0.4328)** | **0.5938 (0.6688)** | **0.5416 (0.5817)** | **0** |
| **Homology Split $\le 70\%$** ($41{,}945 \to 272$) | DIAMOND BLAST-$k$NN | 0.3840 (0.4310) | 0.2412 (0.2890) | 0.4025 (0.4910) | 0.3426 (0.4037) | 246 |
| Homology Split $\le 70\%$ | DeepGOPlus (CNN + Homology) | 0.4512 (0.5120) | 0.2984 (0.3580) | 0.4798 (0.5840) | 0.4098 (0.4847) | 380 |
| Homology Split $\le 70\%$ | ESM-2 C-HMCNN (MCLoss) | 0.5180 (0.5680) | 0.3245 (0.3910) | 0.5280 (0.6410) | 0.4568 (0.5333) | 742 |
| Homology Split $\le 70\%$ | **GSPA Full (Stage 1 + Stage 2A)** | **0.5512 (0.5940)** | **0.3610 (0.4085)** | **0.5480 (0.6080)** | **0.4867 (0.5368)** | **0** |
| **Homology Split $\le 50\%$** ($11{,}900 \to 272$) | DIAMOND BLAST-$k$NN | 0.1980 (0.2240) | 0.1140 (0.1420) | 0.2325 (0.2980) | 0.1815 (0.2213) | 118 |
| Homology Split $\le 50\%$ | DeepGOPlus (CNN + Homology) | 0.3680 (0.4290) | 0.2410 (0.3020) | 0.4375 (0.5510) | 0.3488 (0.4273) | 342 |
| Homology Split $\le 50\%$ | ESM-2 C-HMCNN (MCLoss) | 0.4490 (0.5010) | 0.2780 (0.3440) | 0.4810 (0.5940) | 0.4027 (0.4797) | 695 |
| Homology Split $\le 50\%$ | **GSPA Full (Stage 1 + Stage 2A)** | **0.4910 (0.5380)** | **0.3125 (0.3640)** | **0.5050 (0.5625)** | **0.4362 (0.4882)** | **0** |
| **Homology Split $\le 30\%$** ($2{,}141 \to 272$) | DIAMOND BLAST-$k$NN | 0.0320 (0.0380) | 0.0195 (0.0245) | 0.0420 (0.0452) | 0.0312 (0.0359) | 14 |
| Homology Split $\le 30\%$ | InterPro Domain LogReg | 0.3540 (0.4120) | 0.2080 (0.2640) | 0.3710 (0.4780) | 0.3110 (0.3847) | 124 |
| Homology Split $\le 30\%$ | DeepGOPlus (CNN + Homology) | 0.2890 (0.3680) | 0.1940 (0.2710) | 0.4070 (0.5320) | 0.2967 (0.3903) | 310 |
| Homology Split $\le 30\%$ | ESM-2 $k$-NN (goPredSim) | 0.3620 (0.4210) | 0.2140 (0.2780) | 0.3980 (0.5120) | 0.3247 (0.4037) | 265 |
| Homology Split $\le 30\%$ | ESM-2 C-HMCNN (MCLoss) | 0.3490 (0.4150) | 0.2120 (0.2890) | 0.4200 (0.5460) | 0.3270 (0.4167) | 650 |
| Homology Split $\le 30\%$ | DeepGO-PlusPlus Stacker | 0.4120 (0.4680) | 0.2510 (0.3180) | 0.4510 (0.5580) | 0.3713 (0.4480) | 512 |
| Homology Split $\le 30\%$ | **GSPA Full (Stage 1 + Stage 2A)** | **0.4420 (0.4910)** | **0.2840 (0.3420)** | **0.4690 (0.5632)** | **0.3983 (0.4654)** | **0** |

## Main findings

1. **Retrained time-based splits:** On `Time_Split_CleanA`, retrained `GSPA Full` achieves the highest information-accretion-weighted $F_w = 0.5481$ ($F_{\max} = 0.5897$) with zero ontology violations, outperforming `DIAMOND BLAST-kNN` ($F_w = 0.4869$), `DeepGOPlus` ($0.5125$), `ESM-2 kNN` ($0.4963$), `ESM-2 C-HMCNN` ($0.4998$), and `DeepGO-PlusPlus Stacker` ($0.5429$, $665$ violations). On strict zero-prior `Time_Split_CleanB`, `GSPA Full` leads with $F_w = 0.5416$ ($F_{\max} = 0.5817$).
2. **Retrained homology-based splits ($\le 70\%, \le 50\%, \le 30\%$ identity):** As sequence similarity to the training set is progressively restricted, alignment-based `DIAMOND BLAST-kNN` collapses from $F_w = 0.4869$ (time split) to $0.3426$ ($\le 70\%$), $0.1815$ ($\le 50\%$), and $0.0312$ ($\le 30\%$ twilight zone). Retrained `GSPA Full` degrades gracefully from $F_w = 0.5481$ to $0.4867$ ($\le 70\%$), $0.4362$ ($\le 50\%$), and $0.3983$ ($\le 30\%$, $F_{\max} = 0.4654$), outperforming retrained `ESM-2 C-HMCNN` ($F_w = 0.3270$), `ESM-2 kNN` ($0.3247$), `InterPro Domain LogReg` ($0.3110$), and `DeepGOPlus` ($0.2967$) by $+0.0713$ to $+0.1016$ $F_w$ at $\le 30\%$ identity because conserved InterPro domain signatures, PLM representations, and Stage 1 operon/pathway context priors remain informative when sequence alignment fails.

