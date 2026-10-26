# Cohort recurrence and reaction-role context

This page describes `gspa.integration.cohort`, the reusable core components for audit items 5, 6 and 7: joint inference of whether the genomic context of a candidate (gene, reaction) pair recurs across a cohort of genomes, and how that context is kept separate from the protein's role. The package emits **features**. It computes no posterior, applies no "how often" threshold and adds no fixed boost; those belong to the calibrated learned model and the decision threshold downstream. Everything below was checked on synthetic cohorts only, which shows that the features respond in the intended direction to controlled manipulations. It says nothing about empirical accuracy.

## Inputs (`CohortInput`)

| Input | Notes |
|---|---|
| `baselineSnapshotId` | Required. Names the pre-cohort sequence/evidence snapshot the confidences were taken from. Claims whose provenance starts with `cohort` or contains `feedback` are rejected, so the engine cannot consume its own output. |
| `genomes` (`GenomeEntry`) | `GenomeLayout`, optional `cladeId`, optional `strainClusterId`, optional `AssemblyInfo`, and `familySearchComplete` / `familySearchCoverage` (see below). Protein ids are FAA sequence ids, the same ones the layout uses. |
| `membership` (`FamilyAssignment`) | Per protein, soft weights (each in [0, 1], sum ≤ 1) over orthogroup and optional subfamily. Recurrence is resolved at subfamily level when one is given. |
| `claims` (`CandidateClaim`) | (gene, reaction) candidates with baseline confidence. Weak candidates are kept; no display threshold is applied. |
| `EcCall` + `EcCandidateExpander` | Direct EC → reaction expansion through `ReactionGraph.reactionsForEc`. No EC → GO step. |
| `roleEvidence` | Optional protein-intrinsic raw role scores in [0, 1] (sequence/domain derived) with their `source`. |
| `interactionEdges` | Optional independent, curated reaction-reaction edges. `directed = true` means A → B: A is upstream of B, and the orientation is kept from either side. |
| `reactionDirections` | Optional, per reaction. |
| `organismReactions` | Per genome: `setOrganismReactions(genomeId, reactions)`. A genome without an entry has an **unknown** plausible set; an empty collection is a **known empty** set (`organismSet(g)` returns `null` vs an empty set). |
| `treeNewick` | Optional tree whose leaves are exactly the genome ids, with finite non-negative branch lengths on every non-root node and a terminating `;`. Clades are the maximal subtrees whose largest leaf-to-leaf distance is at most `treeCutDistance`. |

`CohortInference` validates everything before computing and throws `IllegalArgumentException` for: non-finite, negative or > 1 membership weights; weight sums > 1; protein or genome ids that are not in the cohort or its layout; duplicate genomes, claims or assignments; claims with non-finite or out-of-range confidence; role scores outside [0, 1]; half-open interaction edges; contradictory family-search certification; a cluster assigned to several clades; and malformed trees (duplicate, missing or unknown leaves, negative or NaN branches, missing branch lengths, trailing text). Nothing falls back silently to independence. `CohortConfig.validate()` rejects `hopAlpha` outside (0, 1] and other degenerate settings.

The input has no label fields. Features are computed from sequence/evidence candidates, layouts, membership, phylogeny, assembly metadata and independent edges only.

## What is computed

**Reaction signature.** For reaction R the signature is the set of network neighbours within `maxHops` (default 1, via shared non-currency metabolites in `ReactionGraph`) plus curated interaction partners. Each link carries a hop distance, weight `alpha^d`, a flow class (upstream, downstream, bidirectional, curated, distal) that respects reaction direction, and, if an organism reaction set exists, whether the neighbour is plausible in the organism. This is not a named-pathway lookup.

**Windows and anchors.** A window is the genes whose midpoints lie within `halfWidthBp` of the candidate gene. A window is *context-positive* when it contains at least `minAnchors` distinct anchor reactions: confident baseline claims (≥ `anchorMinConfidence`) on other genes that belong to R's signature. The candidate gene itself, any gene of the same family (same-family leak) and R itself are never anchors.

**Peers and holdout.** Features for a gene in genome G use only the peer genomes: G and everything not independent of G (`HoldoutLevel`: genome, strain cluster by default, or clade) are removed. With no peers the output says so (`cohort_no_peers`, `cohort_insufficient_peers`), and nothing is imputed.

**Independent units and effective N.** A clade has total weight `cladeCap`, split equally over its strain clusters; a cluster's weight is split equally over its genomes. Adding a duplicate strain leaves every weighted quantity unchanged, and many clusters in one clade do not exceed one clade's weight. This is deliberately conservative. **Unknown metadata is never read as independence.** A genome whose clade is unknown (after inheriting the clade of its cluster, if that cluster has one) joins ONE unresolved clade block that carries a single clade's weight, however many such genomes there are; a genome with an unknown cluster joins one unresolved cluster per clade. So duplicated genomes with unknown metadata cannot raise `cohort_effective_n` above one block's weight. When no informative peer has a resolved clade, the phylogeny-adjusted outputs are withheld (`cohort_adjusted_available = 0`: no `_lo90`/`_hi90`, no permutation p) while descriptive unadjusted counts stay (`cohort_raw_genomes_*`, `cohort_support_fraction_raw`). Holdout is conservative in the same way: scoring an unknown-clade (or unknown-cluster) target also removes every genome of that unresolved block, so nothing leaks. `cohort_phylogeny_known_fraction` reports how much of the peer set had real metadata. `cohort_effective_n` is the summed weight of informative units and `cohort_coverage` the share of peer clades that were informative.

**Observation states per peer genome.**

| State | Meaning |
|---|---|
| `SUPPORTED` | a family member sits in a context-positive window |
| `NEGATIVE` | family members present, all fully observed windows lack the signature |
| `UNKNOWN_CONTEXT` | family present but the window is cut by a contig end, its extent is unknown, or a member has no coordinate (declared unplaced) |
| `TRUE_ABSENT` | no family member in an assembly marked `COMPLETE` **and** with `familySearchComplete = true` |
| `PATHWAY_WITHOUT_FAMILY` | as above, but the signature reactions are present elsewhere in the genome |
| `UNOBSERVED_ABSENT` | no family member, but the assembly is `DRAFT`/unknown or the family search is not certified over all proteins: unknown, not absence |

A complete assembly alone does not prove homolog absence, because the orthogroup assignment may cover only the input candidates. `familySearchComplete` is the per-genome certification that the search ran over every protein (a `familySearchCoverage` below 1 contradicts it and is rejected). Without certification all missing members are unknown, and copy number is unknown too, so such genomes are left out of the multi-copy and single-copy features.

**Proteins without a coordinate.** A real protein can lack a usable coordinate (absent from the layout, split across contigs, multi-locus). Declare it in `GenomeEntry.unplacedProteins` (metadata `unplacedProteins`) instead of inventing a locus: claims and family assignments on it are accepted, the family counts as **present** with `UNKNOWN_CONTEXT` (never `TRUE_ABSENT`/`UNOBSERVED_ABSENT`/`NEGATIVE`), and its own local features report `local_layout_missing = 1` with no window and no anchor. A genome with one placed copy lacking the signature plus an unresolved copy (unplaced, or in a window that is not fully observed) is `UNKNOWN_CONTEXT`, not `NEGATIVE`, because the unresolved copy may be the supported one; a supported placed copy still makes the genome `SUPPORTED`. Ids that are neither in the layout nor declared remain errors (probable id mismatch), as does declaring an id that has a coordinate. `familySearchComplete` certifies the family search over all proteins and is independent of location completeness, so a complete assembly may hold unplaced proteins without being marked a fragment.

Window state is `FULL`, `TRUNCATED` (known fragment end) or `UNKNOWN`. Missing contig length or closure gives `UNKNOWN`; completeness is never invented. Support is reported with and without the unknown windows (`cohort_support_fraction_unknown_as_negative` / `_as_support`) so a reader can see how much rests on unobservable windows.

**Support, with uncertainty.** `cohort_support_eff`, `cohort_negative_eff`, `cohort_effective_n`, the Jeffreys-smoothed `cohort_support_fraction` (0.5 prior; a singleton never reaches 1) and an *approximate* 90% Wilson interval on effective counts (`_lo90`, `_hi90`). The interval is not an exact credible interval.

**Matched background.** For each informative peer member, background windows are drawn from the *same genome* with the same window length, gene density (± `densityTolerance` genes), contig-end state and contig closure, excluding windows centred on or containing any family member. `cohort_bg_rate`, `cohort_bg_windows`, `cohort_enrichment_log` and `cohort_support_minus_bg` compare support with that background. Total genome presence is never used as the background.

**Joint, complementary patterns.** A signature reaction counts as jointly covered when independent units with at least `jointMinEffective` total weight carry it next to the family. `cohort_joint_coverage` versus `cohort_mean_single_coverage` gives `cohort_complementarity`: partial contexts in different lineages add up. `cohort_complementary_unannotated` is the weighted share of signature reactions with **no confident candidate in the baseline snapshot** in a peer genome, although other lineages carry them next to the family. This is an UNANNOTATED/unsupported label only: it is not evidence that the reaction is genuinely absent or lost, even in a complete assembly (the gene may be missed, divergent or unassembled), and no negative biological label is derived from it.

**Permutation null.** Within each peer genome, observed member windows are replaced by random matched background windows (same stratum), clade and cluster weights are preserved, and the seed is deterministic (`seed` combined with family, reaction and holdout key). Output: `cohort_perm_observed`, `cohort_perm_null_mean`, `cohort_perm_p` (add-one, `(k+1)/(B+1)`), `cohort_perm_n`. It is omitted when no matched windows exist.

**Paralogs.** Membership is soft and nothing transfers to all paralogs: a gene's evidence comes from its own family assignments, alternatives are returned per family with their own features, `cohort_multicopy_fraction` and the `cohort_single_copy_*` features show how much recurrence is ambiguous.

**Guards.** Empty signature (`signature_empty`), no peers, no informative peers, no matched background and no permutation pool produce explicit flags and omitted features, never `NaN`. Missing features are absent from the map rather than zero.

## Roles (item 7)

The package computes **no role probability, no null/unknown distribution and no role confidence threshold**; the learned stage owns the role posterior and null. `RoleEvidenceSummary` only passes raw protein-intrinsic evidence through: each `FeatureRecord` has `roleEvidence` (a list of items with role, source and the raw value, from supplied role scores and from direct EC calls), `roleEvidenceMissing` (true when there is no evidence at all, which is missing, not "non-catalytic"), and `roleAlternatives` (the roles with any evidence, in enum order, plus `UNKNOWN`, which always stays possible). A direct EC call is listed as raw support that `CATALYTIC` is possible, with its baseline confidence unchanged and its source; it is never converted into a probability of role. Features mirror this as `role_evidence_<role>_raw` (absent when there is no evidence) and `role_evidence_missing`. Context, operon association and cohort recurrence are not inputs and cannot create role evidence. Alternative families and reactions for the same gene stay in `alternatives`.

## Outputs (`CohortResult`)

`records`: one `FeatureRecord` per input (gene, reaction) claim, with `features`, `states`, raw role evidence, alternatives and provenance (engine version, baseline snapshot id, baseline provenance, holdout level and the genomes held out, peer count, family resolution, seed). `aggregates`: one per (family, reaction) over **all** genomes, labelled `holdout = none`; they describe the cohort and must not be used to score a member.

## Wiring into the typed workflow and learned features

*Done in `gspa.integration.contextbridge` and the `cohort_context` stage of `IntegrationWorkflow`; see
[context-workflow.md](context-workflow.md) for the options, metadata schema and the two-pass batch. The list below is what
that wiring had to supply and is kept as the contract.* Required by the caller:

1. Build one `CohortInput` per cohort **from the pre-cohort snapshot** (Stage 0/evidence claims and EC calls before any context stage runs), and set `baselineSnapshotId` to that snapshot's identifier. Never pass confidences from a previous cohort/feedback pass.
2. Metadata the workflow must supply: a layout per genome keyed by FAA seqid; orthogroup (and subfamily where available) membership with weights; clade ids or a tree, and strain-cluster ids (ANI-based clusters are suitable); per-genome `AssemblyInfo` with `status`, and per-contig `length` and `closed`. Also required: per-genome `familySearchComplete = true` only when the orthogroup/subfamily search covered every protein (not only the candidates); without it no absence is called true. Without clade/cluster metadata all such genomes share one unresolved block (conservative, never independent) and the adjusted CI/permutation p are withheld; record `cohort_phylogeny_known_fraction` and `cohort_adjusted_available` in the learned feature table.
3. Reaction knowledge: a `ReactionGraph` with currency set and EC bindings, optional reaction directions, curated interaction edges (with `directed` orientation) and per-genome plausible reaction sets (`setOrganismReactions`; leave a genome out when its set is unknown).
4. Call `new CohortInference(input, graph, config).run()` once per batch, after all genomes' baselines exist and before the learned scorer. Join `FeatureRecord.features` on `(genomeId, proteinId, reactionId)`; proteins use FAA seqids. Pass the maps to the learned model unchanged, treating an absent key as missing, and let calibration choose any "often" threshold. Use `holdout = STRAIN_CLUSTER` or `CLADE` consistently with the train/test split used for the learned model, and never feed `aggregates` as per-gene features.
5. Persist `CohortResult.provenance` and each record's `provenance` in the run audit.
6. Role posterior and null are the learned stage's: feed it `roleEvidence`/`role_evidence_*` and `role_evidence_missing` as raw inputs.
7. A single-genome run produces `cohort_no_peers = 1`; the workflow must route that case to a cohort-free feature set rather than interpret it as negative evidence.

## Cost and limits

No candidate is capped or dropped, so weak candidates survive even for dense EC-derived inputs (a million candidates yield a million records). Cost structure: per genome, window membership and the confident-reaction set of each window are built once; per (genome, reaction), one boolean per window is cached; matched-background lookups use density buckets plus a per-(genome, family) set of contaminated windows, so a lookup scales with bucket size and not with the number of genes; peer features are cached per (family, reaction, holdout key). The remaining growth is roughly (distinct family-reaction pairs) x (peer genomes), plus O(genes) memory per distinct (genome, reaction) pair; the permutation null multiplies the per-peer cost by `permutations`. Everything runs in one process in memory. For very large panels, run the engine per batch of reactions or genomes and concatenate records (peers must still include the whole cohort for the batch's reactions), or set `permutations = 0`. The synthetic scale test only checks that nothing is dropped; no throughput figure is claimed.
