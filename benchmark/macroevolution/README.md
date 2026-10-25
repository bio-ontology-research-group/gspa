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
   batch). A genome has a GO term when at least one protein carries it at posterior >= 0.5; the set is
   propagated to all `is_a`/`part_of` ancestors (`go-dag.tsv`) and stored as a bit vector over the GO
   vocabulary. Per genome and configuration: number of annotated proteins, number of distinct terms,
   per-protein counts for the union of the `goslim_prokaryote` and `goslim_metagenomics` slims, and the
   fraction of the bacterial and archaeal essential-function profiles (gspa-core) present. EC numbers
   at >= 0.5 go to a JSON-lines file.
2. **Genome table and species profiles** (`02_species_profiles.py`). `genomes.tsv.gz` joins the batch
   manifests (GTDB species, phylum) and the Stage D per-genome table (ORFs, dark and twilight counts).
   `species_counts.npz` holds, per species and configuration, the number of genomes with each term.
3. **Tree** (`03_tree.py`). Species names map to their GTDB r220 representative (bac120 metadata) and to
   the tips of the GTDB-Tk r220 bac120 reference tree; the tree is pruned to the mapped species and
   rooted at the midpoint of its longest path. Reference species are those with at least 20 Swiss-Prot
   proteins with experimental GO annotations (the Stage 0 DIAMOND database, `train_taxonomy.tsv`,
   matched through NCBI taxids), and every species gets its patristic distance to the nearest reference.
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
   nodes.

## Status

2026-10-24: extraction running on the 242 batches with a Stage D partial (454,691 genomes; three looping
jobs 54494944-6 because the 2,000-job QOS cap left no room for an array). The remaining 454 wave-2
batches are being integrated (P11 Stage D chain 54492370 and successors) and will be added when done.
