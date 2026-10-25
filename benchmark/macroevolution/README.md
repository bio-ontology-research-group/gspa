# Macro-evolution of bacterial function from the GSPA annotation of AllTheBacteria

Uses the Stage D outputs of the AllTheBacteria run (`benchmark/allthebacteria/`, pipeline P11): for every
genome two per-protein posterior tables, `<gid>.noisy_or.tsv.gz` (Stage 0 only: the calibrated, pooled
homology channels) and `<gid>.gspa_full.tsv.gz` (Stages 0 to 2 with all six whole-genome priors and
Stage 2A relaxation). Every analysis below is run on both, and the difference between them is the
measure of what the whole-genome priors change at the scale of bacterial evolution.

Host: KAUST ibex, work area `/encrypted/e3008/gspa/macroevo` (`code/`, `ref/`, `lists/`, `presence/`,
`profiles/`, `tree/`, `results/`, `logs/`). `run.sh` only submits and copies.

## Steps

1. **Presence extraction** (`extract_presence.py`, `01_extract_loop.sbatch`; one process per Stage D
   batch). For every genome, configuration and GO term the maximum posterior with which any protein
   carries the term, after propagation to all `is_a`/`part_of` ancestors (`go-dag.tsv`), stored as a
   uint16 matrix (posterior x 1000) so that the later steps can apply any presence threshold. Per genome,
   configuration and threshold (0.50 and 0.30): the number of proteins with a term at or above it and
   per-protein counts for the union of the `goslim_prokaryote` and `goslim_metagenomics` slims. EC
   numbers at >= 0.50 go to a JSON-lines file. The Stage D posteriors never exceed 0.60, because the
   batch pipeline merges the Stage 0 channels into one evidence source per claim; at the paper's decision
   threshold of 0.50 a genome keeps about 480 terms and the essential functions sit at 0.35 to 0.48 under
   Stage 0 alone, so 0.50 is where the priors act and 0.30 is the sensitivity run (about 1,200 terms per
   genome).
2. **Genome table and species profiles** (`02_species_profiles.py`, `--threshold`). `genomes.tsv.gz`
   joins the batch manifests (species, phylum) and the Stage D per-genome table (ORFs, dark and twilight
   counts) and derives, at the threshold, the number of distinct terms and the fraction of the bacterial
   and archaeal essential-function profiles (gspa-core) present. Species are resolved to GTDB r220 names:
   the manifest name when r220 has it, otherwise the species of the GTDB genome with the same biosample
   (`species_source` column); AllTheBacteria genomes without a GTDB assignment keep the name `unknown`.
   `species_counts.npz` holds, per species and configuration, the number of genomes with each term.
3. **Tree** (`03_tree.py`). Species names map to their GTDB r220 representative (bac120 metadata) and to
   the tips of the GTDB-Tk r220 bac120 reference tree; the tree is pruned to the mapped species and
   rooted at the midpoint of its longest path. Reference species are those with at least one Swiss-Prot
   protein with experimental GO annotations (the Stage 0 DIAMOND database, `train_taxonomy.tsv`; strain
   taxids resolved to species with the NCBI taxdump in `ref/taxdump/`, then matched to GTDB genomes by
   `ncbi_taxid` or `ncbi_species_taxid`); the training set holds 19 bacterial taxa with 6,300 proteins,
   so the reference set is small, and `tips.tsv` keeps the protein count per tip. Every species gets its
   patristic distance to the nearest reference.
   Tip states: a term is present in a species when at least half of its genomes have it. Fitch parsimony
   over all variable terms gives the minimum number of changes per term; one most-parsimonious
   assignment (parent state kept when ambiguous; root takes the majority tip state) gives gains and
   losses per term and per branch.
4. **Tables** (`04_analyses.py`, `results/`): t1 annotation by phylum, t2 by distance to the nearest
   reference species, t3 parsimony summary with essential and taxon-forbidden terms and the paired
   comparison of the configurations, t4 the most labile, most convergent and most changed terms, t5
   power-law exponents of slim categories against genome size (one genome per species), t6 core and
   accessory functions within species with at least 10 genomes, t7 gain and loss rates per unit branch
   length by depth and the busiest branches, t8 functions reconstructed at the root and at phylum crown
   nodes, t9 forbidden terms, essential completeness and term counts split by whether the integrator's
   bundled taxonomy covered the genome's taxon.
5. **Threshold placement** (`05_rescue_floor.py`, `results/t10_rescue_floor.tsv`): for every genome-term
   pair, where full GSPA leaves the posterior relative to the 0.50 decision threshold, by the Stage 0
   value (below the rescue window, inside it, already above), for all terms and for the essential terms.

## Status

2026-10-25: extraction of the 242 batches with a Stage D partial (454,691 genomes) done in three looping
jobs (54509582-4; the 2,000-job QOS cap left no room for an array). Steps 2-4 queued as 54519736
(threshold 0.50, `profiles/ tree/ results/`) and 54519737 (0.30, `*_t30/`). The remaining 454 wave-2
batches are being integrated (P11 Stage D chain 54492370 and successors) and will be added when done.
