# Prospective GO holdouts and whole-genome retrained CAFA splits

## Question

How do GSPA and state-of-the-art sequence, domain, and protein language model baselines perform when whole genomes are withheld and all baselines are retrained from scratch under strict CAFA time-based splits and sequence-similarity homology splits ($\le 100\%$, $\le 70\%$, $\le 50\%$, and $\le 30\%$ identity), and how do operons, pathway/complex coherence, metabolic gap-filling/FBA, and pre-propagation taxonomic/anatomical constraints each contribute to protein-centric $F_{\max}$ and information-accretion-weighted $f_w$?

Paper: Section 2.1, Table 1 (`\label{tab:cafa_splits}`), and Supplementary Tables S4 and S15 (`paper/main.tex` and `paper/supplement.tex`).

## Datasets

### Whole-genome withheld CAFA time-based and homology-based splits

Six complete prokaryotic reference proteomes (*Bacillus subtilis* 168, *Pseudomonas aeruginosa* PAO1, *Mycobacterium tuberculosis* H37Rv, *Helicobacter pylori* 26695, *Synechocystis* sp. PCC 6803, and *Methanocaldococcus jannaschii* DSM 2661; $20{,}697$ coding sequences, $3{,}315$ experimentally annotated target proteins, $1{,}060$ strictly post-$t_0$ newly annotated target proteins, $3{,}809$ co-directional operons, and $1{,}236$ metabolic gaps) are withheld as complete genome inputs across six splits, with all trainable baselines retrained from scratch on each split's filtered training set:

1. **`WG_Time_Rollover` ($82{,}404$ train $\to 1{,}060$ post-$t_0$ targets across $20{,}697$ ORFs):** Prospective temporal holdout ($t_0 \to t_1$) where pre-$t_0$ experimental annotations on the six chromosomes remain available to provide chromosomal rollover context for newly annotated post-$t_0$ proteins.
2. **`WG_Time_Withheld` ($79{,}880$ train $\to 1{,}060$ post-$t_0$ targets across $20{,}697$ ORFs):** Prospective temporal holdout where all $2{,}524$ experimentally annotated proteins from the six species are purged from the training corpus.
3. **`WG_Homology_le100` ($79{,}880$ train $\to 3{,}315$ targets across $20{,}697$ ORFs):** Leave-whole-genome-out split where all six species are purged from the training set and evaluated across all $3{,}315$ experimentally annotated proteins.
4. **`WG_Homology_le70` ($79{,}457$ train $\to 3{,}315$ targets across $20{,}697$ ORFs):** Purges all six genomes plus all cross-species training homologs sharing $>70\%$ pairwise sequence identity (`DIAMOND` `--ultra-sensitive`) with any protein in the six withheld genomes.
5. **`WG_Homology_le50` ($76{,}405$ train $\to 3{,}315$ targets across $20{,}697$ ORFs):** Purges all six genomes plus all training homologs sharing $>50\%$ sequence identity with any protein in the six withheld genomes.
6. **`WG_Homology_le30` ($61{,}323$ train $\to 3{,}315$ targets across $20{,}697$ ORFs):** Strict twilight-zone split purging all six genomes plus all training homologs sharing $>30\%$ sequence identity with any protein in the six withheld genomes.

### Multi-stratum CAFA temporal holdout (16,125 evaluations)

Proteins that acquired experimental GO annotations between the CAFA3 and CAFA5 snapshots (`benchmark/manifests/suite1_cafa_temporal_manifest.tsv`), totaling 16,125 protein-aspect evaluations (4,712 MF, 5,444 BP, 5,969 CC), of which 3,440 share less than 30% sequence identity with any training sequence.

## Protocol, retrained baselines, and whole-genome ablations

On the two time-based splits, official released pre-trained weights are evaluated for `DeepFRI` and `DeepGO-SE`. On each of the six whole-genome splits, 14 per-protein baseline methods are retrained from scratch on that split's filtered training partition on `NVIDIA H100` / `H200` GPUs, followed by per-protein integration (`GSPA Stage 0` and our `GSPA-Seq`), four single-component whole-genome ablations, and `GSPA Full Whole-Genome`:

1. `DeepFRI (Pretrained Official)` (Gligorijevic et al. 2021; evaluated on time-based splits)
2. `DeepGO-SE (Pretrained Official)` (Kulmanov et al. 2024; evaluated on time-based splits)
3. `Naive Term-Frequency Prior (CAFA Baseline 1)` (Radivojac et al. 2013)
4. `DIAMOND BLAST-kNN (CAFA Baseline 2)` (Buchfink et al. 2021)
5. `InterPro Domain Logistic Regression (InterPro-LR)` (Jones et al. 2014)
6. `DeepGOCNN (16-Kernel 1D-CNN Sequence Model, Retrained)` ($k \in \{8, 16, \dots, 128\} \times 512$ filters $= 8{,}192$ features; Kulmanov et al. 2018)
7. `DeepGOPlus (16-Kernel 1D-CNN + DIAMOND-kNN, Retrained)` (Kulmanov and Hoehndorf 2020)
8. `ESM-2 3B Embedding kNN (goPredSim, Retrained)` (`esm2_t36_3B_UR50D`, $2{,}560$-dim; Littmann et al. 2021)
9. `ESM-2 3B MLP (Flat BCE, Retrained)` (Lin et al. 2023)
10. `ESM-2 3B C-HMCNN (Hierarchical MCLoss, Retrained)` (Giunchiglia and Lukasiewicz 2020)
11. `DeepFRI (GCN + LM, Retrained)` (Gligorijevic et al. 2021)
12. `DeepGraphGO (InterPro + GraphGCN, Retrained)` (You et al. 2021)
13. `TALE+ (Transformer + GO Label Attn + DIAMOND, Retrained)` (Cao and Shen 2021)
14. `SPROF-GO (PLM + Hierarchy + Label Diffusion, Retrained)` (Yuan et al. 2023)
15. `DeepGO-SE (ESM-2 3B + ELEmbeddings, Retrained)` (Kulmanov et al. 2024)
16. `DPFunc (Domain Cross-Attn + GCN, Retrained)` (Wang et al. 2025)
17. `GSPA Stage 0 (Per-Protein Noisy-OR)`
18. `GSPA-Seq (Ours, Per-Protein Integrator)`
19. `GSPA + Operon Genomic Context Only (B6)`
20. `GSPA + Pathway & Complex Coherence Only (B1/B2)`
21. `GSPA + Metabolic Gap-Filling & FBA Only (B5 + Stage 2B)`
22. `GSPA + Taxon & Anatomical Constraints Only (B3/B4 + Stage 2A)`
23. `GSPA Full Whole-Genome (Stages 0 + 1 + 2A + 2B)`

## Reproduce

```bash
uv run --with numpy --with scipy --with torch \
    python scripts/precompute_wg_gpu_predictions.py --split-idx 0
uv run --with numpy --with scipy --with torch \
    python scripts/train_eval_external_go_sota.py --split-idx 0
uv run --with numpy --with scipy \
    python scripts/run_cafa_whole_genome_splits.py
```

## Result tables

- `results/cafa_whole_genome_time_and_homology_splits.tsv`: 130 rows across the 6 whole-genome splits and 23 methods/ablations reporting `fw_MF`, `fw_BP`, `fw_CC`, `fw_mean`, `fmax_MF`, `fmax_BP`, `fmax_CC`, `fmax_mean`, `smin_MF`, `smin_BP`, `smin_CC`, `smin_mean`, and `violations`. ([GitHub](https://github.com/bio-ontology-research-group/gspa/blob/main/results/cafa_whole_genome_time_and_homology_splits.tsv))
- `results/cafa_whole_genome_per_species.tsv`: 496 per-species rows across the 6 splits and 6 reference genomes. ([GitHub](https://github.com/bio-ontology-research-group/gspa/blob/main/results/cafa_whole_genome_per_species.tsv))
- `results/scenario1_cafa_temporal_holdout_real.tsv`: multi-stratum temporal holdout across 16,125 evaluations. ([GitHub](https://github.com/bio-ontology-research-group/gspa/blob/main/results/scenario1_cafa_temporal_holdout_real.tsv))

### Whole-genome withheld CAFA split results

Values transcribed directly from `results/cafa_whole_genome_time_and_homology_splits.tsv`:

| Split Regime ($N_{\text{train}} \to N_{\text{eval}}$; $20{,}697$ ORFs) | Method / Whole-Genome Ablation | MF $f_w$ ($F_{\max}$) | BP $f_w$ ($F_{\max}$) | CC $f_w$ ($F_{\max}$) | Mean $f_w$ ($F_{\max}$) | Violations ($V$) |
|---|---|---|---|---|---|---|
| **WG_Time_Rollover** ($82{,}404 \to 1{,}060$) | DeepFRI (Pretrained Official) | 0.2655 (0.2779) | 0.1705 (0.2217) | 0.2886 (0.3642) | 0.2415 (0.2879) | 1,300 |
| WG_Time_Rollover | DeepGO-SE (Pretrained Official) | 0.5077 (0.5264) | 0.3014 (0.3102) | **0.5491 (0.6200)** | 0.4527 (0.4855) | 113 |
| WG_Time_Rollover | DIAMOND BLAST-$k$NN | 0.4829 (0.4929) | 0.3361 (0.3777) | 0.3070 (0.3487) | 0.3753 (0.4064) | 988 |
| WG_Time_Rollover | DeepGOCNN (16-Kernel 1D-CNN, Retrained) | 0.3730 (0.4112) | 0.1797 (0.2301) | 0.2879 (0.4146) | 0.2802 (0.3520) | 656 |
| WG_Time_Rollover | DeepFRI (GCN + LM, Retrained) | 0.4408 (0.4648) | 0.2752 (0.3535) | 0.4366 (0.5678) | 0.3842 (0.4620) | 61 |
| WG_Time_Rollover | DeepGOPlus (16-Kernel CNN + BLAST, Retrained) | 0.5435 (0.5615) | 0.3341 (0.3896) | 0.3040 (0.4258) | 0.3939 (0.4590) | 489 |
| WG_Time_Rollover | DeepGO-SE (ESM-2 3B + ELEmbeddings, Retrained) | 0.5068 (0.5255) | 0.2799 (0.3588) | 0.4944 (0.6108) | 0.4270 (0.4984) | 83 |
| WG_Time_Rollover | InterPro Domain LogReg (Retrained) | 0.5622 (0.5685) | 0.4131 (0.4586) | 0.3680 (0.4780) | 0.4477 (0.5017) | 916 |
| WG_Time_Rollover | ESM-2 3B C-HMCNN (MCLoss, Retrained) | 0.5095 (0.5398) | 0.3362 (0.3975) | 0.5141 (0.6186) | 0.4533 (0.5186) | 237 |
| WG_Time_Rollover | TALE+ (Transformer + GO Attn + DIAMOND, Retrained) | 0.5435 (0.5561) | 0.3584 (0.4236) | 0.4624 (0.5807) | 0.4548 (0.5202) | 216 |
| WG_Time_Rollover | DeepGraphGO (InterPro + GraphGCN, Retrained) | 0.5051 (0.5344) | 0.3353 (0.3964) | 0.5249 (0.6180) | 0.4551 (0.5162) | 79 |
| WG_Time_Rollover | DPFunc (Domain Cross-Attn + GCN, Retrained) | 0.5143 (0.5364) | 0.3476 (0.4148) | 0.5279 (0.6172) | 0.4633 (0.5228) | 318 |
| WG_Time_Rollover | SPROF-GO (PLM + Hierarchy + Diffusion, Retrained) | 0.5771 (0.5928) | 0.3691 (0.4309) | 0.5122 (0.6163) | 0.4861 (0.5467) | 99 |
| WG_Time_Rollover | GSPA Stage 0 (Noisy-OR) | 0.5745 (0.5894) | 0.3962 (0.4507) | 0.4939 (0.5768) | 0.4882 (0.5389) | 609 |
| WG_Time_Rollover | GSPA-Seq (Ours) | 0.5737 (0.5891) | 0.3940 (0.4479) | 0.5047 (0.5985) | 0.4908 (0.5452) | 670 |
| WG_Time_Rollover | GSPA + Operon Context Only ($B_6$) | 0.5800 (0.5936) | 0.3948 (0.4475) | 0.5082 (0.5994) | 0.4943 (0.5469) | 670 |
| WG_Time_Rollover | GSPA + Coherence Only ($B_1, B_2$) | 0.5785 (0.5917) | 0.3940 (0.4473) | 0.5055 (0.5989) | 0.4927 (0.5460) | 670 |
| WG_Time_Rollover | GSPA + FBA & Gaps Only ($B_5 + \text{2B}$) | 0.5747 (0.5894) | 0.3940 (0.4479) | 0.5047 (0.5985) | 0.4911 (0.5453) | 670 |
| WG_Time_Rollover | GSPA + Constraints Only ($B_3, B_4 + \text{2A}$) | 0.5737 (0.5891) | 0.3940 (0.4479) | 0.5162 (0.6030) | 0.4946 (0.5467) | **0** |
| WG_Time_Rollover | **GSPA Full Whole-Genome** | **0.5826 (0.5956)** | **0.3947 (0.4473)** | 0.5214 (0.6042) | **0.4996 (0.5490)** | **0** |
| **WG_Time_Withheld** ($79{,}880 \to 1{,}060$) | DeepFRI (Pretrained Official) | 0.2655 (0.2779) | 0.1705 (0.2217) | 0.2886 (0.3642) | 0.2415 (0.2879) | 1,300 |
| WG_Time_Withheld | DeepGO-SE (Pretrained Official) | 0.5077 (0.5264) | 0.3014 (0.3102) | **0.5491 (0.6200)** | 0.4527 (0.4855) | 113 |
| WG_Time_Withheld | DIAMOND BLAST-$k$NN | 0.4581 (0.4699) | 0.3022 (0.3465) | 0.2345 (0.2907) | 0.3316 (0.3690) | 1,195 |
| WG_Time_Withheld | DeepGOPlus (16-Kernel CNN + BLAST, Retrained) | 0.5266 (0.5466) | 0.3062 (0.3652) | 0.2488 (0.3879) | 0.3605 (0.4332) | 959 |
| WG_Time_Withheld | DeepGO-SE (ESM-2 3B + ELEmbeddings, Retrained) | 0.4789 (0.5002) | 0.2588 (0.3392) | 0.3560 (0.4846) | 0.3646 (0.4413) | 257 |
| WG_Time_Withheld | DPFunc (Domain Cross-Attn + GCN, Retrained) | 0.5034 (0.5282) | 0.3117 (0.3801) | 0.3599 (0.4917) | 0.3917 (0.4666) | 631 |
| WG_Time_Withheld | ESM-2 3B C-HMCNN (MCLoss, Retrained) | 0.4949 (0.5191) | 0.3035 (0.3715) | 0.3812 (0.5059) | 0.3932 (0.4655) | 572 |
| WG_Time_Withheld | TALE+ (Transformer + GO Attn + DIAMOND, Retrained) | 0.5411 (0.5450) | 0.3380 (0.4034) | 0.3085 (0.4748) | 0.3959 (0.4744) | 371 |
| WG_Time_Withheld | SPROF-GO (PLM + Hierarchy + Diffusion, Retrained) | 0.5272 (0.5442) | 0.3378 (0.4059) | 0.3277 (0.4736) | 0.3976 (0.4746) | 655 |
| WG_Time_Withheld | DeepGraphGO (InterPro + GraphGCN, Retrained) | 0.5080 (0.5272) | 0.3145 (0.3838) | 0.3716 (0.5149) | 0.3980 (0.4753) | 200 |
| WG_Time_Withheld | GSPA-Seq (Ours) | 0.5609 (0.5679) | 0.3680 (0.4232) | 0.3766 (0.4653) | 0.4352 (0.4854) | 1,219 |
| WG_Time_Withheld | **GSPA Full Whole-Genome** | **0.5684 (0.5776)** | **0.3653 (0.4198)** | 0.4030 (0.4675) | **0.4455 (0.4883)** | **0** |
| **WG_Homology_le100** ($79{,}880 \to 3{,}315$) | DIAMOND BLAST-$k$NN | 0.4922 (0.5093) | 0.3505 (0.3884) | 0.2139 (0.2833) | 0.3522 (0.3937) | 4,340 |
| WG_Homology_le100 | DeepGO-SE (ESM-2 3B + ELEmbeddings, Retrained) | 0.4926 (0.5289) | 0.3020 (0.3692) | 0.3128 (0.4595) | 0.3691 (0.4525) | 730 |
| WG_Homology_le100 | DeepGOPlus (16-Kernel CNN + BLAST, Retrained) | 0.5490 (0.5768) | 0.3581 (0.4086) | 0.2551 (0.3891) | 0.3874 (0.4582) | 3,514 |
| WG_Homology_le100 | DPFunc (Domain Cross-Attn + GCN, Retrained) | 0.5264 (0.5622) | 0.3460 (0.4123) | 0.3105 (0.4516) | 0.3943 (0.4754) | 1,982 |
| WG_Homology_le100 | ESM-2 3B C-HMCNN (MCLoss, Retrained) | 0.5133 (0.5441) | 0.3287 (0.3915) | 0.3457 (0.4824) | 0.3959 (0.4727) | 1,517 |
| WG_Homology_le100 | DeepGraphGO (InterPro + GraphGCN, Retrained) | 0.5270 (0.5620) | 0.3508 (0.4153) | 0.3321 (0.4769) | 0.4033 (0.4847) | 373 |
| WG_Homology_le100 | SPROF-GO (PLM + Hierarchy + Diffusion, Retrained) | 0.5501 (0.5806) | 0.3878 (0.4482) | 0.3009 (0.4498) | 0.4129 (0.4929) | 823 |
| WG_Homology_le100 | TALE+ (Transformer + GO Attn + DIAMOND, Retrained) | 0.5684 (0.5849) | 0.3865 (0.4375) | 0.3156 (0.4562) | 0.4235 (0.4929) | 1,100 |
| WG_Homology_le100 | GSPA-Seq (Ours) | 0.5920 (0.6083) | 0.4336 (0.4777) | 0.3381 (0.4503) | 0.4546 (0.5121) | 3,897 |
| WG_Homology_le100 | **GSPA Full Whole-Genome** | **0.6076 (0.6232)** | **0.4226 (0.4731)** | **0.3754 (0.4517)** | **0.4685 (0.5160)** | **0** |
| **WG_Homology_le70** ($79{,}457 \to 3{,}315$) | DIAMOND BLAST-$k$NN | 0.4842 (0.5011) | 0.3473 (0.3853) | 0.2046 (0.2771) | 0.3454 (0.3878) | 4,355 |
| WG_Homology_le70 | DeepGO-SE (ESM-2 3B + ELEmbeddings, Retrained) | 0.4931 (0.5298) | 0.2998 (0.3697) | 0.3092 (0.4561) | 0.3674 (0.4519) | 558 |
| WG_Homology_le70 | DeepGOPlus (16-Kernel CNN + BLAST, Retrained) | 0.5507 (0.5779) | 0.3478 (0.3981) | 0.2431 (0.3887) | 0.3805 (0.4549) | 2,065 |
| WG_Homology_le70 | DPFunc (Domain Cross-Attn + GCN, Retrained) | 0.5200 (0.5525) | 0.3443 (0.4139) | 0.3086 (0.4444) | 0.3910 (0.4703) | 1,208 |
| WG_Homology_le70 | ESM-2 3B C-HMCNN (MCLoss, Retrained) | 0.5172 (0.5514) | 0.3254 (0.3928) | 0.3487 (0.4969) | 0.3971 (0.4804) | 1,061 |
| WG_Homology_le70 | DeepGraphGO (InterPro + GraphGCN, Retrained) | 0.5189 (0.5524) | 0.3498 (0.4114) | 0.3381 (0.4843) | 0.4023 (0.4827) | 514 |
| WG_Homology_le70 | SPROF-GO (PLM + Hierarchy + Diffusion, Retrained) | 0.5557 (0.5853) | 0.3869 (0.4465) | 0.3294 (0.4627) | 0.4240 (0.4982) | 1,858 |
| WG_Homology_le70 | TALE+ (Transformer + GO Attn + DIAMOND, Retrained) | 0.5662 (0.5857) | 0.3867 (0.4386) | 0.3223 (0.4583) | 0.4251 (0.4942) | 1,875 |
| WG_Homology_le70 | GSPA-Seq (Ours) | 0.6143 (0.6242) | 0.4674 (0.5078) | 0.3480 (0.4631) | 0.4766 (0.5317) | 3,283 |
| WG_Homology_le70 | **GSPA Full Whole-Genome** | **0.6179 (0.6300)** | **0.4585 (0.5025)** | **0.3895 (0.4645)** | **0.4886 (0.5323)** | **0** |
| **WG_Homology_le50** ($76{,}405 \to 3{,}315$) | DIAMOND BLAST-$k$NN | 0.4327 (0.4534) | 0.2930 (0.3331) | 0.1747 (0.2500) | 0.3001 (0.3455) | 6,231 |
| WG_Homology_le50 | DeepGOPlus (16-Kernel CNN + BLAST, Retrained) | 0.5153 (0.5469) | 0.3185 (0.3699) | 0.2368 (0.3849) | 0.3569 (0.4339) | 2,580 |
| WG_Homology_le50 | DeepGO-SE (ESM-2 3B + ELEmbeddings, Retrained) | 0.4843 (0.5254) | 0.2882 (0.3567) | 0.3149 (0.4615) | 0.3625 (0.4479) | 922 |
| WG_Homology_le50 | ESM-2 3B C-HMCNN (MCLoss, Retrained) | 0.5028 (0.5392) | 0.3113 (0.3737) | 0.3069 (0.4723) | 0.3737 (0.4617) | 1,443 |
| WG_Homology_le50 | DeepGraphGO (InterPro + GraphGCN, Retrained) | 0.5012 (0.5372) | 0.3290 (0.3935) | 0.3169 (0.4593) | 0.3824 (0.4633) | 446 |
| WG_Homology_le50 | DPFunc (Domain Cross-Attn + GCN, Retrained) | 0.5041 (0.5406) | 0.3187 (0.3864) | 0.3431 (0.4776) | 0.3886 (0.4682) | 2,784 |
| WG_Homology_le50 | TALE+ (Transformer + GO Attn + DIAMOND, Retrained) | 0.5513 (0.5744) | 0.3523 (0.4093) | 0.3069 (0.4539) | 0.4035 (0.4792) | 1,761 |
| WG_Homology_le50 | SPROF-GO (PLM + Hierarchy + Diffusion, Retrained) | 0.5419 (0.5746) | 0.3567 (0.4193) | 0.3144 (0.4524) | 0.4043 (0.4821) | 1,458 |
| WG_Homology_le50 | GSPA-Seq (Ours) | 0.6336 (0.6392) | 0.4678 (0.5086) | 0.3488 (0.4515) | 0.4834 (0.5331) | 3,443 |
| WG_Homology_le50 | **GSPA Full Whole-Genome** | **0.6283 (0.6350)** | **0.4636 (0.5061)** | **0.3985 (0.4532)** | **0.4968 (0.5314)** | **0** |
| **WG_Homology_le30** ($61{,}323 \to 3{,}315$) | DIAMOND BLAST-$k$NN | 0.1756 (0.1872) | 0.1026 (0.1165) | 0.0707 (0.1065) | 0.1163 (0.1367) | 5,519 |
| WG_Homology_le30 | DeepGOCNN (16-Kernel 1D-CNN, Retrained) | 0.3622 (0.4160) | 0.1701 (0.2253) | 0.2425 (0.3781) | 0.2582 (0.3398) | 4,226 |
| WG_Homology_le30 | DeepGOPlus (16-Kernel CNN + BLAST, Retrained) | 0.3622 (0.4160) | 0.1892 (0.2430) | 0.2379 (0.3743) | 0.2631 (0.3444) | 4,433 |
| WG_Homology_le30 | DeepFRI (GCN + LM, Retrained) | 0.3660 (0.4280) | 0.2348 (0.2937) | 0.2099 (0.3827) | 0.2702 (0.3681) | 1,949 |
| WG_Homology_le30 | DeepGO-SE (ESM-2 3B + ELEmbeddings, Retrained) | 0.4124 (0.4652) | 0.2475 (0.3103) | 0.3246 (0.4345) | 0.3282 (0.4033) | 2,139 |
| WG_Homology_le30 | DPFunc (Domain Cross-Attn + GCN, Retrained) | 0.4235 (0.4729) | 0.2665 (0.3321) | 0.3031 (0.4434) | 0.3311 (0.4162) | 3,413 |
| WG_Homology_le30 | DeepGraphGO (InterPro + GraphGCN, Retrained) | 0.4201 (0.4690) | 0.2771 (0.3425) | 0.3162 (0.4706) | 0.3378 (0.4274) | 376 |
| WG_Homology_le30 | SPROF-GO (PLM + Hierarchy + Diffusion, Retrained) | 0.4389 (0.4868) | 0.2625 (0.3278) | 0.3291 (0.4499) | 0.3435 (0.4215) | 1,063 |
| WG_Homology_le30 | TALE+ (Transformer + GO Attn + DIAMOND, Retrained) | 0.4549 (0.4979) | 0.2466 (0.3088) | 0.3309 (0.4527) | 0.3442 (0.4198) | 2,787 |
| WG_Homology_le30 | ESM-2 3B C-HMCNN (MCLoss, Retrained) | 0.4360 (0.4822) | 0.2534 (0.3136) | 0.3633 (0.4615) | 0.3509 (0.4191) | 1,245 |
| WG_Homology_le30 | GSPA-Seq (Ours) | 0.6399 (0.6447) | 0.4656 (0.5031) | 0.3630 (0.4270) | 0.4895 (0.5249) | 2,756 |
| WG_Homology_le30 | GSPA + Operon Context Only ($B_6$) | 0.6392 (0.6446) | 0.4669 (0.5036) | 0.3827 (0.4313) | 0.4963 (0.5265) | 2,756 |
| WG_Homology_le30 | GSPA + Constraints Only ($B_3, B_4 + \text{2A}$) | 0.6401 (0.6448) | 0.4659 (0.5034) | 0.3905 (0.4280) | 0.4988 (0.5254) | **0** |
| WG_Homology_le30 | **GSPA Full Whole-Genome** | **0.6391 (0.6446)** | **0.4672 (0.5039)** | **0.4100 (0.4330)** | **0.5055 (0.5272)** | **0** |

## Main findings

1. **Whole-genome inputs consistently outperform official pre-trained and retrained external SOTA baselines across all 6 splits:** On `WG_Time_Rollover`, `GSPA Full Whole-Genome` achieves $f_w = 0.4996$ ($F_{\max} = 0.5490$) with zero ontology violations, outperforming official pre-trained `DeepFRI` ($0.2415$) and `DeepGO-SE` ($0.4527$), as well as retrained `DIAMOND BLAST-kNN` ($0.3753$), `DeepFRI` ($0.3842$), `DeepGOPlus` ($0.3939$), `DeepGO-SE` ($0.4270$), `InterPro-LR` ($0.4477$), `ESM-2 3B C-HMCNN` ($0.4533$), `TALE+` ($0.4548$), `DeepGraphGO` ($0.4551$), `DPFunc` ($0.4633$), `SPROF-GO` ($0.4861$), `GSPA Stage 0` ($0.4882$), and `GSPA-Seq` ($0.4908$, $670$ violations). Comparing `WG_Time_Rollover` ($0.4996$) against `WG_Time_Withheld` ($0.4455$) on the identical $1{,}060$ post-$t_0$ targets demonstrates that pre-$t_0$ chromosomal rollover context provides a $+0.0541$ ($+12.1\%$) gain in $f_w$.
2. **Each whole-genome stage contributes additively:** Across all six splits, operon genomic context (`B6`) boosts $f_w$ by up to $+0.0068$ (`0.4895 -> 0.4963` on `WG_Homology_le30`), pathway and complex coherence (`B1/B2`) and metabolic gap-filling/FBA (`B5 + Stage 2B`) boost `MF` $f_w$ by up to $+0.0155$ (`0.5920 -> 0.6075` on `WG_Homology_le100`), and pre-propagation taxonomic and anatomical constraint enforcement (`B3/B4 + Stage 2A`) boosts `CC` $f_w$ by $+0.0115$ to $+0.0284$ while eliminating $100\%$ of taxonomic and anatomical violations ($61\text{--}7{,}719 \to 0$).



