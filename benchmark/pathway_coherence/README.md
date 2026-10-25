# Pathway coherence: evidence-gated KEGG module completeness and module-based refinement

## Problem with the v2.1.0 metric

`PathwayDatabase.computePathwayCoherenceDetailed` treats a pathway as *triggered* when a single one of its
enzymes is annotated and then scores the fraction of its enzymes present. Shared enzymes of central
metabolism trigger dozens of pathways in every genome, KEGG maps are enzyme bags without step logic, and
an enzyme whose EC number has no GO mapping counts as missing. The number reports breadth of enzyme
coverage, not whether the organism can run a pathway. The coherence prior then boosted candidates for any
missing enzyme of any triggered pathway.

## What this adds

* **Module definitions with step logic** (`fetch_kegg_module_definitions.py` → `data/kegg_module_defs.json`):
  582 KEGG modules with their DEFINITION expressions parsed into required steps, alternative orthologs,
  complexes with optional subunits and KEGG's own gaps, plus KO → EC and the organisms KEGG lists as
  complete (`complete_organisms`). `data/kegg_organism_modules.json` holds gene-to-module links for the
  nine validation organisms. 2,549 steps; 69 % map to GO through ec2go, the rest are *unmappable* and are
  reported as such instead of as missing functions.
* **Assessment** (`gspa.ontology.ModuleDefinitions`, `gspa.ontology.ModuleAssessor`): per step a status
  (present, partial complex, sub-threshold evidence, missing, unmappable, forbidden for the taxon), module
  completeness over assessable required steps, and an activation rule weighted by enzyme specificity
  (1 / number of modules requiring the enzyme) so shared enzymes activate nothing on their own. Signature
  modules are assessed but never scored. `gspa evaluate --module-defs` adds `module_coherence`,
  `modules_activated`, `modules_complete` and `module_detail` to the quality JSON and keeps the old number
  under `pathway_coherence` and `pathway_map_coverage` as a breadth measure.
* **Refinement** (`CoherencePrior`, mode `module`): only activated modules at least half complete are
  refined, one candidate per open step and at most two steps per module; a candidate is a sub-threshold
  claim on a step enzyme scored by Stage 0 log-odds, corroboration, an unannotated-protein bonus, a penalty
  for a protein already confidently annotated elsewhere in the same aspect, an operon-partner bonus and the
  module's completeness. The chosen claim is pinned to 0.60. Steps without a candidate stay open.
* **Provenance**: `gspa integrate --module-defs` writes `<out>.modules.json` with every module before and
  after refinement, the proteins behind each present step, the annotations that changed with their
  evidence, and each unresolved gap with its best sub-threshold candidate and the reason it was left
  alone; rescues also go to `<out>.rescues.tsv`.

## Validation (`run_bench9.py`, unimatrix01 `/data/hohndor/gspa/pathway_coherence/`)

Nine reference genomes, three variants (`legacy`, `assess` = definitions loaded without module
refinement, `module`), judged by F-max against curated annotations, coherence-rescue precision against the
curated annotation of the rescued protein, and module-level agreement with the organisms KEGG lists as
complete. Results in `results/`.

## Round 1 findings (`results/` from job 9692, before the assessability rules)

Refinement was neutral on F-max (micro all 0.7734 to 0.7744, CAFA all 0.8057 to 0.8064, experimental
truth within 0.001) and made 331 coherence rescues across the nine genomes, 131 of the 324 with a curated
protein correct at term level (40 %). Module completeness judged against KEGG's complete-organism lists
reached recall 0.95 in E. coli but precision 0.50: 105 modules called complete that KEGG does not list.
Three causes, two of them fixable:

* 29 of the 105 had a single GO-mappable step and 26 had fewer than half their steps mappable, so one
  enzyme made them "complete". Fix: a module is *assessable* only with at least two mappable required
  steps covering at least half of its steps; the rest are reported as not assessable.
* Partial EC numbers (1.1.1.-) map through ec2go to enzyme-class terms (hydro-lyase, acid-thiol ligase,
  methyltransferase activity) that every genome carries, so coarse steps counted as present. Fix:
  partial ECs are treated as no mapping (337 of 4,107 orthologs).
* KEGG names organism-specific variants of one chemistry (isoprenoid biosynthesis for archaea, plants and
  non-plant eukaryotes; C1-unit interconversion, eukaryotes) whose orthologs differ below EC resolution,
  so a GO-based annotation satisfies all variants. Fix: the module name's domain hint against the genome's
  domain marks such modules as not expected; variants without a name hint remain indistinguishable and
  are a floor on the precision measured against KEGG.

Round 2 (`bench9_v2`) adds the `nocoh` baseline (all priors except coherence) so the annotations each
pathway refinement adds, and their precision against the curated truth, are measured directly.

## Rounds 2 to 4: safeguards and final validation

**Round 2** (`results_v2/`, after the assessability rules) showed which rescues are wrong
(`analyse_rescues.py`, `results_v2/rescue_features.tsv`): of 271 judged module rescues, those whose
candidate had an evidence-adjusted log-odds below zero were 17 % correct, those in modules that KEGG never
finds complete in another organism of the genome's domain 19 % (for example keratan and chondroitin sulfate
degradation, animal pathways that two shared glycosidases had activated in E. coli), and those on proteins
with a confident competing annotation 28 %. Two rules follow from what a rescue asserts, a posterior of 0.60,
rather than from fitting: the candidate's adjusted odds must be at least even, and the module must be
complete in at least five other organisms of the domain in KEGG (organism context; the test organism and,
for E. coli, both E. coli strains are left out of the count). On round 2 they keep 88 of the 124 correct
rescues and raise rescue precision from 46 % to 74 %. The organism context alone removes 56 of E. coli's 105
modules called complete without KEGG listing them, and none of its 103 correctly complete modules.

**Round 3** (`results_v3/`) added the `nocoh` baseline (all priors except coherence). The module rescues
were 69 % correct among the annotations they newly put above 0.50, but the coherence prior's has_part
process rescues and complex-partner boosts, unchanged so far, added 9 correct and 49 wrong annotations
(`results_v3/added_breakdown.tsv`). In module mode these lifts now need the same candidate evidence.

**Round 4** (`results_v4/`, final; nine genomes, curated truth from `truth_dual`):

| | F-max micro, all | F-max CAFA, all | F-max micro, exp | added correct | added wrong |
|---|---|---|---|---|---|
| no coherence prior (`nocoh`) | 0.7746 | 0.8066 | 0.2865 | | |
| v2.1.0 coherence prior (`legacy`) | 0.7734 | 0.8057 | 0.2864 | 191 | 674 |
| module mode (`module`) | 0.7748 | 0.8067 | 0.2866 | 42 | 29 |

Against the v2.1.0 prior, module mode removes 96 % of the wrong additions and has the higher micro F-max in
all nine genomes. Against no pathway refinement, it adds 42 correct functions for 29 wrong ones (precision
0.59, the posterior it assigns is 0.60); micro F-max on all annotations is equal or higher in eight genomes
and 0.0002 lower in M. tuberculosis. The gain over no refinement is small: refinement recovers a few dozen
functions per nine genomes, and the main improvement is the removal of the v2.1.0 prior's false additions.
Module completeness judged against KEGG's complete-module lists (assessable, expected modules only):
precision 0.81 and recall 0.78 at Stage 0, 0.78 and 0.90 after all priors in module mode
(`results_v4/runs.tsv`, `complete_*` columns). Per genome and per module rows are in `results_v4/`.

## Held-out module steps: recovery of correct functions (`held_out_steps.py`, `results_held_out/`)

Round 4 compares whole annotations, where pathway refinement touches few claims. To test recovery
directly, in each of the nine genomes one correct step per near-complete module was hidden: a step whose
annotated protein carries the step's function in the curated annotation (`truth_all`) at Stage 0 >= 0.50,
in an activated, expected module at least two-thirds complete. All claims of that protein on the step's GO
terms (with descendants and the EC numbers mapping to them) were set to a calibrated probability drawn from
U[0.25, 0.45], and the perturbed claims were integrated with refinement off (`nocoh`), with the v2.1.0
coherence prior and in module mode (leave-one-out organism context as before). 774 steps were hidden.

| | hidden steps recovered (>= 0.50) | other annotations added, correct | other annotations added, wrong |
|---|---|---|---|
| refinement off | 233 (30.1 %) | | |
| v2.1.0 coherence prior | 488 (63.0 %) | 215 | 688 |
| module mode | 384 (49.6 %) | 31 | 22 |

Module mode recovers 151 more hidden correct functions than no refinement (+65 %), more in every genome
(`results_held_out/held_out_runs.tsv`), while adding 22 wrong annotations elsewhere. The v2.1.0 prior
recovers more hidden steps, 104 more than module mode, but adds 688 wrong annotations to get them: about
6.6 wrong additions per extra correct recovery, against 0.15 for module mode. The 30 % recovered without
refinement comes from the other priors and from claims whose remaining support stays above threshold.
The targets were chosen from module-mode assessments of the unperturbed run, so they are steps the module
assessment can see; the comparison between variants uses the same targets.

## Final configuration (v2.2.2) and checks on the live pipeline

Round 4 showed that the coherence prior's has_part process rescues and complex-partner boosts add mostly wrong
annotations even when gated (29 % precise). From v2.2.2, module mode applies only the module-step rescues and the
complex-singleton penalty. Round 5 (`results_v5/`) and the held-out rerun with the same 774 targets
(`results_held_out_v2/`):

| | F-max micro, all | F-max CAFA, all | added correct | added wrong | hidden steps recovered | other wrong added |
|---|---|---|---|---|---|---|
| refinement off | 0.7746 | 0.8066 | | | 233 (30.1 %) | |
| v2.1.0 coherence prior | 0.7734 | 0.8057 | 191 | 674 | 488 (63.0 %) | 688 |
| module mode, v2.2.2 | 0.7748 | 0.8067 | 37 | 17 | 379 (49.0 %) | 10 |

On the unperturbed genomes refinement adds 37 correct and 17 wrong annotations (precision 0.69). With correct
steps hidden it recovers 146 more of them than no refinement, in every genome, for 10 wrong annotations elsewhere.

**Live pipeline (`live_check/`).** `validate_live.py` ran the deployed image on the E. coli K-12 MG1655 example
input against the curated annotation (RefSeq to UniProt). As deployed (GO predictions below 0.5 dropped before
integration), module refinement and no refinement give identical annotations: E. coli's module enzymes are
annotated from Swiss-Prot, which contains E. coli, and the 76 open steps of activated modules have no supporting
candidate. Passing the DIAMOND transfers between 30 % and 50 % identity through as claims produced one rescue and
lowered precision from 0.782 to 0.778 through the other priors, so the service keeps its 0.5 floor.

**Live examples on gspa.bio2vec.net.** `examples/ecoli` (E. coli K-12 MG1655): 148 activated modules, 105 complete,
no refinement change (see above). `examples/pdistasonis` (Parabacteroides distasonis, public AllTheBacteria assembly
SAMD00389893, 2.8 substitutions per site from the nearest experimentally annotated species): 13 rescues (12 module
steps, 1 essential function), 53 of 121 activated modules complete after refinement against 43 before, mean module
completeness 0.73 to 0.78; each change is listed with its evidence on the Pathways tab. This strain has no curated
annotation, so these changes are not validated individually; the precision of such rescues is the one measured above.
