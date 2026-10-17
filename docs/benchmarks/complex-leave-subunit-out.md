# Heteromeric complex leave-subunit-out

## Question

Can complex stoichiometry and co-directional operon architecture distinguish the true gene encoding a missing complex subunit from domain-matched decoy competitors when sequence similarity has diverged into the twilight zone?

Paper: Section 2.3, Table 2 (in latex `tab:complex_lso`), and Figure 4a to 4b (`paper/main.tex`).

## Datasets

### 175-case prokaryotic leave-subunit-out panel

175 trials across six prokaryotic genomes, evaluated on curated obligate heteromeric macromolecular complexes drawn from Complex Portal and CORUM:

- *Escherichia coli* K-12 MG1655 (30 trials, calibration genome)
- *Bacillus subtilis* 168 (30 trials, Gram-positive model)
- *Mycobacterium tuberculosis* H37Rv (30 trials, actinobacterium)
- *Pseudomonas aeruginosa* PAO1 (30 trials, gammaproteobacterium)
- *Synechocystis* sp. PCC 6803 (25 trials, photosynthetic cyanobacterium with dispersed subunits)
- *Escherichia coli* O157:H7 EDL933 (30 trials, enteric pathogen)

### Experimental design

In each trial:

1. One subunit of an obligate heteromeric complex is designated as held-out.
2. Its sequence-level confidence score is attenuated to simulate divergence into the low-identity twilight zone (\(\le 0.35\)).
3. Five domain-matched decoy proteins drawn from unrelated loci elsewhere in the host chromosome are added as competitors (yielding six candidate loci per trial).
4. Predictors must rank the true subunit above the decoys and assign high posterior confidence.

## Protocol

Four methods are evaluated across the 175 trials:

1. **Baseline Noisy-OR:** Unconstrained Stage 0 candidate ranking based on attenuated sequence and domain evidence alone.
2. **Stage 1 CoherencePrior:** Adds the genome-scale heteromeric complex coherence prior (\(B_2\)), which detects incomplete complexes and applies an odds boost to all candidate proteins matching the missing subunit's domain role.
3. **Stage 1 Coherence + Operon Context:** Combines the global stoichiometry prior with the co-directional operon context prior (\(B_6\)), boosting candidates physically adjacent to confirmed partner subunits in co-directional operon clusters.
4. **Stage 1 + Stage 2 + Stage 3 Agent:** Full GSPA pipeline incorporating Stage 2A obligate-complex singleton pruning and Stage 3 tool-mediated disambiguation for non-operonic, dispersed candidates.

### Metrics

- **Hit@1:** Fraction of trials where the true subunit is ranked Rank 1 among the six candidates.
- **Recall at 0.50:** Fraction of true subunits recovered with posterior probability \(\ge 0.50\).
- **Precision and \(F_1\) at 0.50:** Precision and harmonic mean of precision and recall at posterior threshold 0.50.
- **Complex coherence:** Genome-wide fraction of obligate complexes with all required subunits identified.

## Baselines

The study evaluates the sequence-only Stage 0 baseline, genome-wide CoherencePrior alone (global stoichiometry without genomic syntax), and Coherence plus Operon context.

## Reproduce

The benchmark script `benchmark/unified/pillar3/run_pillar3_complex_lso.py` coordinates the evaluation:

```bash
uv run --with numpy --with scipy \
    python benchmark/unified/pillar3/run_pillar3_complex_lso.py \
        --output-dir results/
```

This script reads the curated complex definitions, injects domain-matched decoys, runs each inference configuration, and outputs `results/pillar3_complex_lso.tsv`.

## Result tables

- `results/pillar3_complex_lso.tsv`: columns `genome`, `method`, `n_cases`, `hit_at_1`, `recall_at_05`, `precision_at_05`, `f1_at_05`, `complex_coherence`. ([GitHub](https://github.com/bio-ontology-research-group/gspa/blob/main/results/pillar3_complex_lso.tsv))

### Empirical results

Values transcribed directly from `results/pillar3_complex_lso.tsv`:

| Genome | Method | Trials | Hit@1 | Recall (\(\ge 0.50\)) | Precision (\(\ge 0.50\)) | \(F_1\) (\(\ge 0.50\)) | Complex Coherence |
|---|---|---|---|---|---|---|---|
| *E. coli* K-12 | Baseline Noisy-OR | 30 | 0.5667 | 0.1667 | 0.5000 | 0.2500 | 0.4706 |
| *E. coli* K-12 | Stage 1 CoherencePrior | 30 | 0.7000 | 0.5667 | 0.2833 | 0.3778 | 0.7938 |
| *E. coli* K-12 | Stage 1 Coherence + Operon | 30 | 0.9667 | 0.9667 | 0.6267 | 0.7089 | 0.7938 |
| *E. coli* K-12 | Stage 1 + 2 + 3 Agent | 30 | 0.9667 | 1.0000 | 0.6842 | 0.7619 | 0.8235 |
| *B. subtilis* 168 | Baseline Noisy-OR | 30 | 0.6000 | 0.2000 | 0.5455 | 0.2927 | 0.3378 |
| *B. subtilis* 168 | Stage 1 CoherencePrior | 30 | 0.7000 | 0.5333 | 0.2759 | 0.3636 | 0.6471 |
| *B. subtilis* 168 | Stage 1 Coherence + Operon | 30 | 1.0000 | 0.9667 | 0.6306 | 0.7090 | 0.6618 |
| *B. subtilis* 168 | Stage 1 + 2 + 3 Agent | 30 | 1.0000 | 0.9667 | 0.6905 | 0.7632 | 0.6912 |
| *M. tuberculosis* H37Rv | Baseline Noisy-OR | 30 | 0.4667 | 0.1333 | 0.4444 | 0.2051 | 0.3857 |
| *M. tuberculosis* H37Rv | Stage 1 CoherencePrior | 30 | 0.6000 | 0.4333 | 0.2167 | 0.2889 | 0.7797 |
| *M. tuberculosis* H37Rv | Stage 1 Coherence + Operon | 30 | 1.0000 | 0.9667 | 0.6589 | 0.7244 | 0.7797 |
| *M. tuberculosis* H37Rv | Stage 1 + 2 + 3 Agent | 30 | 1.0000 | 0.9667 | 0.7073 | 0.7733 | 0.7966 |
| *P. aeruginosa* PAO1 | Baseline Noisy-OR | 30 | 0.6000 | 0.1667 | 0.5000 | 0.2500 | 0.3776 |
| *P. aeruginosa* PAO1 | Stage 1 CoherencePrior | 30 | 0.7000 | 0.5333 | 0.2667 | 0.3556 | 0.6897 |
| *P. aeruginosa* PAO1 | Stage 1 Coherence + Operon | 30 | 0.9667 | 0.9667 | 0.6106 | 0.6862 | 0.7126 |
| *P. aeruginosa* PAO1 | Stage 1 + 2 + 3 Agent | 30 | 0.9667 | 0.9667 | 0.6744 | 0.7436 | 0.7356 |
| *Synechocystis* sp. | Baseline Noisy-OR | 25 | 0.5200 | 0.1600 | 0.4444 | 0.2353 | 0.3587 |
| *Synechocystis* sp. | Stage 1 CoherencePrior | 25 | 0.6400 | 0.4800 | 0.2500 | 0.3288 | 0.6125 |
| *Synechocystis* sp. | Stage 1 Coherence + Operon | 25 | 0.6400 | 0.6400 | 0.3393 | 0.4080 | 0.6125 |
| *Synechocystis* sp. | Stage 1 + 2 + 3 Agent | 25 | 0.6800 | 0.6800 | 0.3793 | 0.4480 | 0.6375 |
| *E. coli* O157:H7 | Baseline Noisy-OR | 30 | 0.6000 | 0.1667 | 0.5000 | 0.2500 | 0.4103 |
| *E. coli* O157:H7 | Stage 1 CoherencePrior | 30 | 0.7000 | 0.4667 | 0.2414 | 0.3182 | 0.7416 |
| *E. coli* O157:H7 | Stage 1 Coherence + Operon | 30 | 1.0000 | 0.9667 | 0.6622 | 0.7279 | 0.7416 |
| *E. coli* O157:H7 | Stage 1 + 2 + 3 Agent | 30 | 1.0000 | 0.9667 | 0.7073 | 0.7733 | 0.7640 |
| **All (Weighted Mean)** | **InterProScan / Pfam Domain Only** | **175** | **0.3886** | **0.1143** | **0.2857** | **0.1633** | **0.3120** |
| **All (Weighted Mean)** | **DIAMOND-BLASTp Only** | **175** | **0.4914** | **0.1486** | **0.4483** | **0.2232** | **0.3540** |
| **All (Weighted Mean)** | **Operon Neighborhood Vote (STRING-like)** | **175** | **0.7829** | **0.6514** | **0.4872** | **0.5574** | **0.5840** |
| **All (Weighted Mean)** | **Baseline Noisy-OR** | **175** | **0.5600** | **0.1657** | **0.4915** | **0.2479** | **0.3901** |
| **All (Weighted Mean)** | **Stage 1 CoherencePrior** | **175** | **0.6743** | **0.5029** | **0.2558** | **0.3391** | **0.7107** |
| **All (Weighted Mean)** | **Stage 1 Coherence + Operon** | **175** | **0.9371** | **0.9200** | **0.5951** | **0.6680** | **0.7170** |
| **All (Weighted Mean)** | **Stage 1 + 2 + 3 Agent** | **175** | **0.9429** | **0.9314** | **0.6485** | **0.7185** | **0.7445** |

## Main findings

1. **Outperforming sequence, domain, and neighborhood baselines:** Standalone InterProScan/Pfam (`0.3886` Hit@1), DIAMOND-BLASTp (`0.4914`), and Stage 0 Noisy-OR (`0.5600`) miss most twilight-zone subunits, while a STRING-style operon neighborhood vote achieves `0.7829` Hit@1 (`F1 = 0.5574`).
2. **Global stoichiometry + local operon synteny:** Applying `CoherencePrior` (\(B_2\)) raises Hit@1 from `0.5600` to `0.6743` and recall at 0.50 from `0.1657` to `0.5029`. Coupling complex stoichiometry with co-directional operon context (\(B_6\)) lifts Hit@1 to `0.9371` (\(p < 10^{-15}\), exact McNemar test; `0.9310` across the five out-of-sample species), recall at 0.50 to `0.9200`, and \(F_1\) to `0.6680`.
3. **Stage 2A singleton pruning and Stage 3 agent disambiguation:** Stage 2A obligate-complex singleton pruning and Stage 3 tool-mediated reasoning lift overall Hit@1 to `0.9429`, precision to `0.6485`, \(F_1\) to `0.7185`, and complex coherence to `0.7445` (raising *Synechocystis* sp. Hit@1 from `0.6400` to `0.6800`).

