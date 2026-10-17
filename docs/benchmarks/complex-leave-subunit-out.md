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
| *E. coli* K-12 | Baseline Noisy-OR | 30 | 0.5667 | 0.0000 | 0.0000 | 0.0000 | 0.4706 |
| *E. coli* K-12 | Stage 1 CoherencePrior | 30 | 0.5667 | 0.5667 | 0.1933 | 0.2756 | 0.7938 |
| *E. coli* K-12 | Stage 1 Coherence + Operon | 30 | 0.9667 | 1.0000 | 0.6267 | 0.7089 | 0.7938 |
| *E. coli* K-12 | Stage 1 + 2 + 3 Agent | 30 | 0.9667 | 1.0000 | 0.6267 | 0.7089 | 0.7938 |
| *B. subtilis* 168 | Baseline Noisy-OR | 30 | 0.6000 | 0.0000 | 0.0000 | 0.0000 | 0.3378 |
| *B. subtilis* 168 | Stage 1 CoherencePrior | 30 | 0.6000 | 0.5333 | 0.1972 | 0.2757 | 0.6471 |
| *B. subtilis* 168 | Stage 1 Coherence + Operon | 30 | 1.0000 | 0.9667 | 0.6306 | 0.7090 | 0.6618 |
| *B. subtilis* 168 | Stage 1 + 2 + 3 Agent | 30 | 1.0000 | 0.9667 | 0.6306 | 0.7090 | 0.6618 |
| *M. tuberculosis* H37Rv | Baseline Noisy-OR | 30 | 0.4667 | 0.0000 | 0.0000 | 0.0000 | 0.3857 |
| *M. tuberculosis* H37Rv | Stage 1 CoherencePrior | 30 | 0.4667 | 0.4333 | 0.1256 | 0.1911 | 0.7797 |
| *M. tuberculosis* H37Rv | Stage 1 Coherence + Operon | 30 | 1.0000 | 0.9667 | 0.6589 | 0.7244 | 0.7797 |
| *M. tuberculosis* H37Rv | Stage 1 + 2 + 3 Agent | 30 | 1.0000 | 0.9667 | 0.6589 | 0.7244 | 0.7797 |
| *P. aeruginosa* PAO1 | Baseline Noisy-OR | 30 | 0.6000 | 0.0000 | 0.0000 | 0.0000 | 0.3776 |
| *P. aeruginosa* PAO1 | Stage 1 CoherencePrior | 30 | 0.6000 | 0.5333 | 0.1772 | 0.2529 | 0.6897 |
| *P. aeruginosa* PAO1 | Stage 1 Coherence + Operon | 30 | 0.9667 | 0.9667 | 0.6106 | 0.6862 | 0.7126 |
| *P. aeruginosa* PAO1 | Stage 1 + 2 + 3 Agent | 30 | 0.9667 | 0.9667 | 0.6106 | 0.6862 | 0.7126 |
| *Synechocystis* sp. | Baseline Noisy-OR | 25 | 0.5200 | 0.0000 | 0.0000 | 0.0000 | 0.3587 |
| *Synechocystis* sp. | Stage 1 CoherencePrior | 25 | 0.5200 | 0.4800 | 0.1813 | 0.2507 | 0.6125 |
| *Synechocystis* sp. | Stage 1 Coherence + Operon | 25 | 0.6400 | 0.6400 | 0.3393 | 0.4080 | 0.6125 |
| *Synechocystis* sp. | Stage 1 + 2 + 3 Agent | 25 | 0.6800 | 0.6800 | 0.3793 | 0.4480 | 0.6125 |
| *E. coli* O157:H7 | Baseline Noisy-OR | 30 | 0.6000 | 0.0000 | 0.0000 | 0.0000 | 0.4103 |
| *E. coli* O157:H7 | Stage 1 CoherencePrior | 30 | 0.6000 | 0.4667 | 0.1639 | 0.2302 | 0.7416 |
| *E. coli* O157:H7 | Stage 1 Coherence + Operon | 30 | 1.0000 | 0.9667 | 0.6622 | 0.7279 | 0.7416 |
| *E. coli* O157:H7 | Stage 1 + 2 + 3 Agent | 30 | 1.0000 | 0.9667 | 0.6622 | 0.7279 | 0.7416 |
| **All (Weighted Mean)** | **Baseline Noisy-OR** | **175** | **0.5600** | **0.0000** | **0.0000** | **0.0000** | **0.3901** |
| **All (Weighted Mean)** | **Stage 1 CoherencePrior** | **175** | **0.5600** | **0.5029** | **0.1729** | **0.2459** | **0.7107** |
| **All (Weighted Mean)** | **Stage 1 Coherence + Operon** | **175** | **0.9371** | **0.9257** | **0.5951** | **0.6680** | **0.7170** |
| **All (Weighted Mean)** | **Stage 1 + 2 + 3 Agent** | **175** | **0.9429** | **0.9314** | **0.6009** | **0.6737** | **0.7170** |

## Main findings

1. **Global stoichiometry alone cannot rank candidates:** Applying the genome-wide CoherencePrior alone increases subunit recall at 0.50 from 0.0000 to 0.5029 and lifts complex completeness from 0.3901 to 0.7107, but leaves Hit@1 unchanged at the baseline level of 0.5600. Because all candidates matching the required domain role receive the same genome-wide boost, stoichiometry flags incomplete assemblies without discriminating among competing loci.
2. **Local operon synteny resolves locus ambiguity:** Coupling the complex coherence prior with co-directional operon context lifts Hit@1 from 0.5600 to 0.9371 (\(p < 10^{-15}\), exact McNemar test; 0.9310 across the five out-of-sample species), recall at 0.50 to 0.9257, and \(F_1\) to 0.6680. Stoichiometry determines *what* function is missing; chromosomal operon synteny determines *where* on the chromosome it is encoded.
3. **Dispersed complexes in cyanobacteria:** In *Synechocystis* sp. PCC 6803, photosystem and phycobilisome subunits are often encoded at dispersed genomic loci outside operons. Operon synteny increases Hit@1 from 0.5200 to 0.6400, and Stage 3 tool-mediated reasoning (inspecting Pfam architectures and pathway neighborhood evidence) lifts Hit@1 further to 0.6800 and \(F_1\) to 0.4480.
