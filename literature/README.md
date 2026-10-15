# Literature Collection for Genome-Scale Protein Annotation (GSPA)

This directory contains verified literature reviews, BibTeX reference files, and 68 open-access full-text PDFs supporting the GSPA manuscript.

## Synthesis Reports & BibTeX Files

1. **[01_genome_context_and_metabolic_gapfilling.md](01_genome_context_and_metabolic_gapfilling.md)** ([01_references.bib](01_references.bib))
   - Operon prediction, gene neighborhood/synteny conservation (Overbeek et al. 1999; Dandekar et al. 1998; Price et al. 2005)
   - Phylogenetic profiling (Pellegrini et al. 1999; Chen & Vitkup 2006)
   - Genomic language models: gLM (Hwang et al. 2024), gLM2 (Cornman et al. 2024), Evo (Nguyen et al. 2024), context-aware GNNs (Toibazar et al. 2024)
   - Metabolic network reconstruction & gap-filling: gapseq (Zimmermann et al. 2021), CarveMe (Machado et al. 2018), ModelSEED2 (Seaver et al. 2021), Pathway Tools Hole Filler (Green & Karp 2004), CanOE (Smith et al. 2012), and the parsimony gap-filling benchmark (Karp et al. 2018)

2. **[02_ontology_taxon_anatomical_and_cafa.md](02_ontology_taxon_anatomical_and_cafa.md)** ([02_references.bib](02_references.bib))
   - GO taxon constraints (`only_in_taxon`, `never_in_taxon`; Deegan et al. 2010)
   - Anatomical and developmental structure constraints in multicellular eukaryotes: GO–Uberon cross-product definitions (`results_in_development_of`, `occurs_in`; Mungall et al. 2011, 2012; Hill et al. 2010, 2013; Khodiyar et al. 2011)
   - GO Causal Activity Models (GO-CAMs; Thomas et al. 2019)
   - CAFA 1–5 evaluation protocols, Information Accretion ($IA$), IA-weighted $F_{\max}$, $S_{\min}$ (Clark & Radivojac 2013; Radivojac et al. 2013; Jiang et al. 2016; Zhou et al. 2019)
   - SOTA protein-level GO predictors: DeepGOZero (Kulmanov & Hoehndorf 2022), DeepGO-SE (Kulmanov et al. 2024), DeepFRI (Gligorijevic et al. 2021), Metagenomic-DeepFRI (Maranga et al. 2023 / Bezshapkin et al. 2026), NetGO 3.0 (Wang et al. 2023), SPROF-GO (Yuan et al. 2023), TALE (Cao & Shen 2021), TransFun (Boadu et al. 2023)

3. **[03_genome_annotators_complexes_and_pathways.md](03_genome_annotators_complexes_and_pathways.md)** ([03_references.bib](03_references.bib))
   - Genome annotation pipelines and their protein-by-protein independence assumption: NCBI PGAP (Tatusova et al. 2016; Li et al. 2021), Bakta (Schwengers et al. 2021), Prokka (Seemann 2014), eggNOG-mapper v2 (Cantalapiedra et al. 2021), InterProScan 5 (Jones et al. 2014; Paysan-Lafosse et al. 2023), DRAM (Shaffer et al. 2020), MicrobeAnnotator (Ruiz-Perez et al. 2021), KofamKOALA (Aramaki et al. 2020), Mantis (Queiros et al. 2021), RASTtk (Brettin et al. 2015)
   - Protein complex & pathway coherence: Complex Portal (Meldal et al. 2015, 2022), CORUM (Giurgiu et al. 2019), IntAct (Orchard et al. 2014; del-Toro et al. 2022), KEGG Modules (Kanehisa et al. 2016, 2023), MinPath (Ye & Doak 2009)

4. **[04_agentic_and_llm_annotation_workflows.md](04_agentic_and_llm_annotation_workflows.md)** ([04_references.bib](04_references.bib))
   - Tool-using LLM agents for gene set and function analysis: GeneAgent (Wang et al. 2025, Nature Methods), ProtAgents (Ghafarollahi & Buehler 2024, Digital Discovery), CRISPR-GPT (Qu et al. 2025, Nature Biomedical Engineering), The Virtual Lab (Swanson et al. 2025, Nature), ChemCrow (Bran et al. 2024, Nature Machine Intelligence)
   - LLM-based protein function / GO annotation and biocuration: OntoGPT / SPIRES (Caufield et al. 2024, Bioinformatics), GORetriever & Lit-KNN in NetGO 3.0 (Wang et al. 2023; Yan et al. 2024, Bioinformatics), ProtST (Xu et al. 2023, ICML), ProteinChat (Huo et al. 2024; Guo et al. 2023), Evolla (Zhou et al. 2025), ProLLaMA (Lv et al. 2025, IEEE Trans. AI), ProtGO (Wang et al. 2025, Bioinformatics), InstructProtein (Wang et al. 2024, ACL)
   - Agentic and LLM-assisted metabolic network curation & gap-filling: Human2 (Luo et al. 2025, PNAS), Biochemical Network Benchmarking (Tewari et al. 2026, eLife), MechAInistic (Loecker et al. 2026)
   - Scalable cascade routing & selective prediction: FrugalGPT (Chen et al. 2023, NeurIPS), BioMistral (Labrak et al. 2024, ACL), GSPA Three-Tier Adaptive Cascade Routing Architecture
   - Neuro-symbolic guardrails: ELK reasoner (Kazakov et al. 2014, JAR), DeepGOZero & DeepGO-SE (Kulmanov et al. 2022, 2024), Taxon Constraints (Deegan et al. 2010), Uberon cross-products (Mungall et al. 2012)

5. **[05_sota_8scenarios_and_4kingdoms.md](05_sota_8scenarios_and_4kingdoms.md)** ([05_references.bib](05_references.bib))
   - SOTA baselines, benchmark datasets, and computational efficiency across all 8 GSPA scenarios:
     * Scenario 1 (Suite 1): Prospective CAFA & enzyme temporal holdouts (NetGO 3.0, DeepGO-SE, SPROF-GO, AnnoPRO, TALE, ESM-2 kNN/MLP)
     * Scenario 2 (Suite 2): Multi-kingdom proteome-wide annotation across Bacteria, Archaea, Viruses, and Eukaryotes under Homology Masking Protocol H and Out-of-Clade Protocol D; viral Baltimore-class and essential core gene constraints
     * Scenarios 3, 4, 5 (Suites 3, 4, 5): Role-stratified leave-reaction-out (LRO) metabolic gap-filling across Core Biomass Stratum A, Alternative Catabolic Stratum B (B1 vs B2), and Secondary/Peripheral Stratum C; mathematical audit of reaction MILP defects vs GSPA StoichiometricLocusRefiner
     * Scenario 6 (Suite 6): Experimentally validated difficult bacterial enzymes (RB-TnSeq Price-149) across 6 SOTA EC predictors (CLEAN, DeepProZyme, EnzBERT, MAPred, GraphEC, TopEC)
     * Scenarios 7 & 8 (Suites 7 & 8): Large-scale phenotype and diagnostic enzyme reproduction in BacDive (4,788 strains, 171,984 assays; Track 1 blind genome-to-phenotype and Track 2 phenotype-assisted dark-matter attribution)
     * Dedicated sections on Obligate Heteromeric Complex Completion (CORUM, Complex Portal, STRING) and Multicellular Eukaryotic Anatomical/Taxonomic Constraint Enforcement (GO taxon constraints, Uberon/PO/CL cross-products, ELK reasoner)
     * Comprehensive computational efficiency benchmark (wall-clock runtime and peak RAM per genome across 25+ tools)

## Full-Text PDFs ([pdfs/](pdfs/))

63 verified open-access PDFs are stored in [`pdfs/`](pdfs/).
