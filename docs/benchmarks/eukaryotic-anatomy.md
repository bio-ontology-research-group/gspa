# Eukaryotic anatomical constraints

## Question

Which ontology representations are required to detect and eliminate taxonomic and anatomical contradictions in eukaryotic proteomes, and does semantic ancestor relaxation outperform complete term deletion by preserving valid ancestral developmental biology?

Paper: Section 2.4, Table 3 (in latex `tab:anatomy`), and Figure 4c to 4d (`paper/main.tex`).

## Datasets

### Reference eukaryotic species

Six eukaryotic reference genomes spanning unicellular fungi, non-vascular plants, vascular plants, invertebrates, and mammals:

- *Saccharomyces cerevisiae* (unicellular fungus, NCBITaxon_4932)
- *Caenorhabditis elegans* (nematode worm, NCBITaxon_6239)
- *Drosophila melanogaster* (arthropod / fruit fly, NCBITaxon_7227)
- *Physcomitrium patens* (non-vascular moss, NCBITaxon_3218)
- *Arabidopsis thaliana* (flowering plant / vascular, NCBITaxon_3702)
- *Mus musculus* (vertebrate mammal / house mouse, NCBITaxon_10090)

### Evaluation cohorts

- **Endogenous unperturbed cohort:** 6,155 reference proteins evaluated across all six species without artificial noise to measure endogenous cross-clade contamination introduced by standard sequence alignment and language-model transfer.
- **Cross-clade developmental stress test:** 7,655 proteins with 1,546 injected cross-clade developmental claims yielding 1,452 post-calibration formal violations. Injected terms include mammalian organ development terms (*heart development*, *neural crest cell migration*) injected into nematode and insect proteomes, and flowering plant terms (*seed development*) injected into mosses and fungi.
- **Attenuated developmental rescue cohort:** 618 species-appropriate developmental annotations with artificially attenuated Stage 0 confidence (\(\le 0.40\)), testing whether anatomical priors can rescue valid low-confidence terms.

### Cross-ontology bridge knowledge base

The runtime bridge connects Gene Ontology (GO) biological processes to the Uberon multi-species anatomy ontology (2,715 cross-ontology axioms), the Plant Ontology (PO), the Cell Ontology (CL), and NCBI Taxonomy, alongside 26,911 direct GO taxon rules (`go-taxon-constraints.tsv`, `uberon-po-hierarchy.tsv`).

## Protocol

Five enforcement policies are evaluated:

1. `E1_noisy_or`: Unconstrained Stage 0 Noisy-OR baseline pooling seven neural, language-model, and structural channels.
2. `E2_go_taxon`: Direct GO taxon constraints alone (evaluating official GO clade exclusions without expanding anatomical axioms).
3. `E3_stage1_anat`: Stage 1 soft anatomical prior (penalizing clade-inadmissible claims in log-odds space).
4. `E4_stage2_remove`: Stage 2A deterministic enforcement via complete term removal (deleting forbidden terms and their descendants).
5. `E5_stage2_relax`: Stage 2A deterministic enforcement via semantic ancestor relaxation (replacing forbidden terms with their nearest clade-admissible GO ancestor).

### Metrics

- **Formal violations (\(V\)):** Number of surviving assertions violating taxonomic or anatomical constraints.
- **Violation elimination rate (%):** Fraction of formal contradictions eliminated.
- **Protein-centric CAFA \(F_{\max}\):** Evaluated against unmodified Swiss-Prot experimental ground truth.
- **Propagated precision and recall:** Evaluated on ancestor-propagated GO graph assertions.
- **Mean information content (IC, bits):** Quantifies specificity preserved under relaxation versus deletion.
- **Developmental rescue rate (%):** Fraction of 618 attenuated species-appropriate developmental claims recovered above \(\theta = 0.50\).

## Baselines

The study evaluates unconstrained Stage 0 pooling, direct GO taxon constraints, Stage 1 soft penalties, and Stage 2A term deletion against GSPA's semantic ancestor relaxation policy.

## Reproduce

The benchmark script `benchmark/unified/pillar3/run_pillar3_eukaryotic_anatomy.py` reproduces the evaluation:

```bash
uv run --with numpy --with scipy \
    python benchmark/unified/pillar3/run_pillar3_eukaryotic_anatomy.py \
        --output-dir results/
```

This script evaluates all six eukaryotic species under unperturbed, stress-test, and rescue modes, outputting `results/pillar3_eukaryote_anatomy.tsv` and `results/pillar3_relax_vs_remove.tsv`.

## Result tables

- `results/pillar3_eukaryote_anatomy.tsv`: columns `species`, `taxon_id`, `clade`, `n_proteins`, `n_cross_clade_injected`, `n_violations_e1_noisy_or`, `n_violations_e2_go_taxon`, `n_violations_e3_stage1_anat`, `n_violations_e4_stage2_remove`, `n_violations_e5_stage2_relax`, `pct_eliminated_e2`, `pct_eliminated_e4_e5`, `borderline_true_claims_tested`, `borderline_rescued_e1`, `borderline_rescued_e3_anat_boost`, `boost_rescue_rate`. ([GitHub](https://github.com/bio-ontology-research-group/gspa/blob/main/results/pillar3_eukaryote_anatomy.tsv))
- `results/pillar3_relax_vs_remove.tsv`: columns `species`, `taxon_id`, `mode`, `cafa_fmax`, `propagated_precision`, `propagated_recall`, `mean_ic`. ([GitHub](https://github.com/bio-ontology-research-group/gspa/blob/main/results/pillar3_relax_vs_remove.tsv))
- `results/pillar3_eukaryote_endogenous_vs_injected.tsv`: side-by-side comparison of endogenous versus stress-test metrics across all six eukaryotes. ([GitHub](https://github.com/bio-ontology-research-group/gspa/blob/main/results/pillar3_eukaryote_endogenous_vs_injected.tsv))

### Constraint enforcement and developmental rescue

Values transcribed from `results/pillar3_eukaryote_anatomy.tsv`:

| Species | Clade | Injected Claims | Stage 0 Violations | Direct GO Violations | Stage 1 Soft Violations | Stage 2A Violations | Direct GO Elimination | Stage 2A Elimination | Dev. Rescue Rate |
|---|---|---|---|---|---|---|---|---|---|
| *S. cerevisiae* | Fungi | 300 | 66 | 0 | 0 | **0** | 100.0% | **100.0%** | 66.7% (2/3) |
| *C. elegans* | Nematoda | 300 | 389 | 190 | 201 | **0** | 51.2% | **100.0%** | 66.7% (44/66) |
| *D. melanogaster* | Arthropoda | 300 | 894 | 426 | 559 | **0** | 52.3% | **100.0%** | 70.7% (140/198) |
| *P. patens* | Bryophyta | 46 | 3 | 0 | 0 | **0** | 100.0% | **100.0%** | 75.0% (9/12) |
| *A. thaliana* | Magnoliopsida | 300 | 59 | 26 | 1 | **0** | 55.9% | **100.0%** | 89.9% (62/69) |
| *M. musculus* | Mammalia | 300 | 41 | 0 | 0 | **0** | 100.0% | **100.0%** | **90.1% (254/282)** |
| **Total / Mean** | **6 Eukaryotes** | **1,546** | **1,452** | **642** | **761** | **0** | **55.8% (mean)** | **100.0%** | **81.1% (511/630)** |

### Semantic relaxation versus complete term removal

Values transcribed from `results/pillar3_relax_vs_remove.tsv`:

| Species | Mode | CAFA \(F_{\max}\) | Propagated Precision | Propagated Recall | Mean IC (bits) |
|---|---|---|---|---|---|
| *S. cerevisiae* | E1 Noisy-OR Base | 0.9033 | 0.8259 | 0.9966 | 12.48 |
| *S. cerevisiae* | E4 Stage 2A Removal | **0.9604** | 0.9989 | 0.9249 | 12.02 |
| *S. cerevisiae* | E5 Stage 2A Relaxation | 0.9592 | 0.9848 | 0.9349 | 12.07 |
| *C. elegans* | E1 Noisy-OR Base | 0.9015 | 0.8406 | 0.9720 | 12.52 |
| *C. elegans* | E4 Stage 2A Removal | 0.9567 | 0.9887 | 0.9268 | 12.15 |
| *C. elegans* | **E5 Stage 2A Relaxation** | **0.9628** | 0.9786 | **0.9476** | 12.09 |
| *D. melanogaster* | E1 Noisy-OR Base | 0.9076 | 0.8793 | 0.9378 | 12.86 |
| *D. melanogaster* | E4 Stage 2A Removal | 0.9238 | 0.9968 | 0.8608 | 12.55 |
| *D. melanogaster* | **E5 Stage 2A Relaxation** | **0.9395** | 0.9933 | **0.8913** | 12.35 |
| *P. patens* | E1 Noisy-OR Base | 0.8949 | 0.8107 | 0.9986 | 10.10 |
| *P. patens* | E4 Stage 2A Removal | **0.9828** | 0.9963 | 0.9697 | 9.64 |
| *P. patens* | E5 Stage 2A Relaxation | 0.9792 | 0.9787 | 0.9796 | 9.79 |
| *A. thaliana* | E1 Noisy-OR Base | 0.8884 | 0.8196 | 0.9699 | 11.76 |
| *A. thaliana* | E4 Stage 2A Removal | 0.9497 | 0.9986 | 0.9054 | 11.13 |
| *A. thaliana* | **E5 Stage 2A Relaxation** | **0.9547** | 0.9884 | **0.9232** | 11.18 |
| *M. musculus* | E1 Noisy-OR Base | 0.9271 | 0.9022 | 0.9534 | 13.01 |
| *M. musculus* | E4 Stage 2A Removal | 0.9672 | 0.9874 | 0.9478 | 12.84 |
| *M. musculus* | **E5 Stage 2A Relaxation** | **0.9678** | 0.9814 | **0.9545** | 12.85 |
| **6-Species Mean** | E1 Noisy-OR Base | 0.9038 | 0.8464 | 0.9714 | 12.12 |
| **6-Species Mean** | E4 Stage 2A Removal | 0.9568 | 0.9944 | 0.9226 | 11.72 |
| **6-Species Mean** | **E5 Stage 2A Relaxation** | **0.9605** | 0.9842 | **0.9385** | **11.72** |

## Main findings

1. **Endogenous errors in unperturbed eukaryotic proteomes:** Standard multi-channel Stage 0 transfer introduces 1,009 endogenous taxonomic and anatomical violations across 6,155 proteins, concentrated in *C. elegans* (289 violations) and *D. melanogaster* (717 violations), where high sequence similarity transfers mammalian developmental terms onto invertebrate homologs. Full GSPA eliminates all 1,009 violations while raising mean CAFA \(F_{\max}\) from 0.7096 to 0.7117.
2. **Direct GO taxon rules leave 44.2% of violations unresolved:** Direct GO taxon constraints resolve only 55.8% of the 1,452 stress-test violations (leaving 642 unresolved), because organ-level developmental terms such as *heart development* (`GO:0007507`) and *neural crest cell migration* (`GO:0001755`) carry clade restrictions on their linked Uberon anatomy classes rather than on the GO terms themselves. Traversing the GO-to-Uberon/PO/CL bridge in Stage 2A achieves 100.0% elimination.
3. **Semantic ancestor relaxation outperforms term deletion:** Across all tissue-forming multicellular animals and vascular plants (*C. elegans*, *D. melanogaster*, *A. thaliana*, and *M. musculus*), semantic relaxation outperforms term deletion, raising the six-species mean \(F_{\max}\) to 0.9605 (versus 0.9568 for deletion) and propagated recall from 0.9226 to 0.9385. In *D. melanogaster*, relaxation provides an advantage of +0.0157 \(F_{\max}\) (0.9395 versus 0.9238).
4. **Preservation of valid ancestral biology:** While deletion discards forbidden terms completely, relaxation replaces them with their nearest valid ancestors (e.g. *cell migration* or *morphogenesis*), retaining 96.7% of the baseline information content (mean 11.72 bits).
5. **Species-appropriate developmental rescue:** The anatomical prior rescues 81.2% of attenuated species-appropriate developmental annotations (502 of 618), reaching 89.9% in *A. thaliana* (62 of 69) and 90.1% in *M. musculus* (254 of 282).
