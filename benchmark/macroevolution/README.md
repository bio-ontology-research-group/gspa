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

## Results (2026-10-25; 458,223 genomes, 11,812 species, 11,438 on the tree)

`results/` holds the run at the decision threshold 0.50 and `results_t30/` the sensitivity run at 0.30
(same code, `THR=0.3`). Numbers below are from those tables; the tables are the record.

**What the priors change at genome level (t1, t9, t10).** The set of annotated proteins is the same under
both configurations (annotated fraction +0.0005), as in the paper's Stage D summary, but the genome-level
function profile is not. At 0.50 full GSPA adds about 25 terms per genome and removes about 8, raises the
bacterial essential-function completeness from 0.27 to 0.48 in every phylum, and removes the eukaryotic
cellular-component terms that homology transfer puts into bacteria: Stage 0 alone places `chloroplast` in
11,436 of the 11,438 species, `nucleus` in 11,432 and `peroxisome` in 9,462 (t4a, t4d). Table t9 splits
the genomes by whether the integrator's bundled NCBI taxonomy covered their taxon: where it did (246,518
genomes), full GSPA removes every forbidden term (3.25 to 0.00 per genome; 99.9 % of genomes improved);
where it did not (211,152 genomes), the forbidden terms stay (2.99 to 2.96). The essential-function gain is
the same in both groups, because the essentiality prior needs no taxonomy. This is the cleanest
demonstration in the set: the same genomes, the same Stage 0 claims, and the prior acts exactly where its
input exists.

**Threshold placement of the essentiality prior (t10).** The rescue floor of the essentiality prior lifts a
borderline candidate "to the decision threshold", and t10 shows where that leaves the genome-level
posteriors: of the 7.0 M (genome, essential term) pairs whose Stage 0 value lies in the rescue window
[0.16, 0.50), full GSPA puts 2.5 M at 0.45-0.49, 1.4 M at 0.49-0.4999, 0.23 M at exactly 0.50 and 2.9 M
above, so 44 % cross the threshold and the rest stop just under it. For all terms in the window 2.6 % cross,
and 98.3 % of the pairs already above 0.50 stay there (the remainder are the suppressed taxon-forbidden
and co-annotation-conflict pairs). Whether an essential function counts as present at 0.50 under full GSPA
is therefore close to a coin flip, which inflates its apparent homoplasy (t3: essential losses 2,628 under
Stage 0, 5,361 under GSPA) and makes essential functions "accessory" within species (t6: 93 species with
an accessory essential function under Stage 0, 567 under GSPA). A floor strictly above the threshold
(0.55, say) would remove this artefact; it is an integrator finding, not a biological one, and the 0.30 run,
where the whole rescue window counts, shows the intended effect without it (essential completeness 0.56 to
0.74, essential terms at the root 16 to 22 of 32, t8).

**Distance to the experimentally annotated organisms (t2).** The Stage 0 training set contains experimental
GO annotations for 19 bacterial taxa (17 of them on the tree; E. coli alone 3,466 proteins). The annotated
fraction of a genome falls with patristic distance from the nearest of them (Spearman -0.60, both
configurations: the priors do not add proteins), and so does the essential-function completeness under
Stage 0 (-0.52). Under full GSPA the completeness no longer depends on that distance (-0.16): the
essentiality prior restores the core functions in the lineages farthest from any reference, and the
restoration is largest beyond one substitution per site from any reference (+0.23).

**Homoplasy (t3, t4).** Fitch parsimony over the 1,344 (Stage 0) and 1,393 (GSPA) terms that vary among
species gives a median consistency index near 0.01 for both; the paired comparison over the 1,262 common
terms is not significant at 0.50 (Wilcoxon p = 0.47) and is at 0.30 (p = 1e-8, GSPA higher). The terms
whose changes GSPA removes almost entirely (t4c) are those a prior makes universal: ribosome biogenesis,
rRNA processing and carboxylic-acid transport go from a scattered presence in 7,000-8,000 species to all
11,437; the terms whose changes GSPA increases (t4d) are the eukaryotic compartment terms, universal
under Stage 0 and removed only in the half of the genomes with a covered taxon, so their remaining
distribution follows taxonomy coverage, not phylogeny. Both directions are integrator effects and should be
read as such; the biological homoplasy sits in the terms unchanged between the configurations (ion
channels, N-acetyltransferases, UDP-glycosyltransferases, t4a).

**Scaling laws (t5).** One genome per species, 11,810 species with at least 500 ORFs: the number of
proteins in a functional category against the number of ORFs on log-log axes reproduces the known
power-law exponents of bacterial functional categories: transcription regulators and signal transduction
scale superlinearly (slopes 1.7 and 1.4), metabolism near 1.4, protein modification 1.7, and the ribosomal
and structural categories are flat (slope 0). The exponents are the same under both configurations for
every category that the priors do not touch; `nucleus` flips from a positive slope under Stage 0 to a
negative one under GSPA because it survives only in the genomes without taxon coverage.

**Branches and ancestral states (t7, t8).** Gains outnumber losses on every depth quartile and on terminal
as well as internal branches (ratio about 1.4), terminal branches carry the highest rate per unit length,
and the busiest internal branches are the crowns of Leptospirae, Cyanobacteriia and Myxococcota (gains)
and a few long terminal pairs with losses (t7b). The root reconstruction at 0.50 carries 341 (Stage 0) and
376 (GSPA) terms, with 4 and 13 of the 32 essential functions, and 16 and 22 at 0.30; the shortfall from 32
at the root is the threshold ceiling of the Stage D posteriors (0.59), which is the limit of the current
batch deployment, not of the method.

**Caveats to state with any of this.** The Stage D posteriors come from the batch pipeline with one merged
evidence source per claim, never exceed 0.59, and reach the 0.50 threshold for only about 480 terms per
genome; AllTheBacteria is dominated by clinical isolates (Pseudomonadota 263 k, Bacillota 92 k genomes;
Escherichia and Salmonella alone 251 k), so species-level profiles with the majority rule are the unit of
every phylogenetic analysis; 41,140 genomes have no GTDB species ("unknown") and are in the genome-level
tables only; the reference set is 17 species; and the 454 wave-2 batches still in Stage D are not yet
included.

## Status

2026-10-25: extraction of the 242 batches with a Stage D partial done in three looping jobs (54509582-4;
the 2,000-job QOS cap left no room for an array). Steps 2-5 run as 54525048 (0.50) and 54525050 (0.30);
tables fetched to `results/` and `results_t30/`. The remaining 454 wave-2 batches are being integrated
(P11 Stage D chain 54492370 and successors) and will be added when done (`run.sh list`, `run.sh extract`,
`run.sh profiles`).
