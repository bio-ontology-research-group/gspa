# Cohort recurrence and reaction-role context

This page describes `gspa.integration.cohort`, the reusable core components for audit items 5, 6 and 7: joint inference of whether the genomic context of a candidate (gene, reaction) pair recurs across a cohort of genomes, and how that context is kept separate from the protein's role. The package emits **features**. It computes no posterior, applies no "how often" threshold and adds no fixed boost; those belong to the calibrated learned model and the decision threshold downstream. Everything below was checked on synthetic cohorts only, which shows that the features respond in the intended direction to controlled manipulations. It says nothing about empirical accuracy.

## Inputs (`CohortInput`)

| Input | Notes |
|---|---|
| `baselineSnapshotId` | Required. Names the pre-cohort sequence/evidence snapshot the confidences were taken from. Claims whose provenance starts with `cohort` or contains `feedback` are rejected, so the engine cannot consume its own output. |
| `genomes` (`GenomeEntry`) | `GenomeLayout`, optional `cladeId`, optional `strainClusterId`, optional `AssemblyInfo`. Protein ids are FAA sequence ids, the same ones the layout uses. |
| `membership` (`FamilyAssignment`) | Per protein, soft weights (sum ≤ 1) over orthogroup and optional subfamily. Recurrence is resolved at subfamily level when one is given. |
| `claims` (`CandidateClaim`) | (gene, reaction) candidates with baseline confidence. Weak candidates are kept; no display threshold is applied. |
| `EcCall` + `EcCandidateExpander` | Direct EC → reaction expansion through `ReactionGraph.reactionsForEc`. No EC → GO step. |
| `roleEvidence` | Optional protein-intrinsic role scores (sequence/domain derived). |
| `interactionEdges` | Optional independent, curated reaction-reaction edges. |
| `reactionDirections`, `organismReactions` | Optional. `organismReactions == null` means unknown, not empty. |
| `treeNewick` | Optional tree; leaves are genome ids. Clades are the maximal subtrees whose largest leaf-to-leaf distance is at most `treeCutDistance`. |

The input has no label fields. Features are computed from sequence/evidence candidates, layouts, membership, phylogeny, assembly metadata and independent edges only.

## What is computed

**Reaction signature.** For reaction R the signature is the set of network neighbours within `maxHops` (default 1, via shared non-currency metabolites in `ReactionGraph`) plus curated interaction partners. Each link carries a hop distance, weight `alpha^d`, a flow class (upstream, downstream, bidirectional, curated, distal) that respects reaction direction, and, if an organism reaction set exists, whether the neighbour is plausible in the organism. This is not a named-pathway lookup.

**Windows and anchors.** A window is the genes whose midpoints lie within `halfWidthBp` of the candidate gene. A window is *context-positive* when it contains at least `minAnchors` distinct anchor reactions: confident baseline claims (≥ `anchorMinConfidence`) on other genes that belong to R's signature. The candidate gene itself, any gene of the same family (same-family leak) and R itself are never anchors.

**Peers and holdout.** Features for a gene in genome G use only the peer genomes: G and everything not independent of G (`HoldoutLevel`: genome, strain cluster by default, or clade) are removed. With no peers the output says so (`cohort_no_peers`, `cohort_insufficient_peers`), and nothing is imputed.

**Independent units and effective N.** A clade has total weight `cladeCap`, split equally over its strain clusters; a cluster's weight is split equally over its genomes. Adding a duplicate strain leaves every weighted quantity unchanged, and many clusters in one clade do not exceed one clade's weight. This is deliberately conservative. A genome without clade metadata is its own clade and without cluster metadata its own cluster; `cohort_phylogeny_known_fraction` reports how much of the peer set had real metadata. `cohort_effective_n` is the summed weight of informative units and `cohort_coverage` the share of peer clades that were informative.

**Observation states per peer genome.**

| State | Meaning |
|---|---|
| `SUPPORTED` | a family member sits in a context-positive window |
| `NEGATIVE` | family members present, all fully observed windows lack the signature |
| `UNKNOWN_CONTEXT` | family present but the window is cut by a contig end or its extent is unknown |
| `TRUE_ABSENT` | family absent from an assembly marked `COMPLETE` |
| `PATHWAY_WITHOUT_FAMILY` | as above, but the signature reactions are present elsewhere in the genome |
| `UNOBSERVED_ABSENT` | family absent from a `DRAFT` or metadata-less assembly: not evidence of absence |

Window state is `FULL`, `TRUNCATED` (known fragment end) or `UNKNOWN`. Missing contig length or closure gives `UNKNOWN`; completeness is never invented. Support is reported with and without the unknown windows (`cohort_support_fraction_unknown_as_negative` / `_as_support`) so a reader can see how much rests on unobservable windows.

**Support, with uncertainty.** `cohort_support_eff`, `cohort_negative_eff`, `cohort_effective_n`, the Jeffreys-smoothed `cohort_support_fraction` (0.5 prior; a singleton never reaches 1) and an *approximate* 90% Wilson interval on effective counts (`_lo90`, `_hi90`). The interval is not an exact credible interval.

**Matched background.** For each informative peer member, background windows are drawn from the *same genome* with the same window length, gene density (± `densityTolerance` genes), contig-end state and contig closure, excluding windows centred on or containing any family member. `cohort_bg_rate`, `cohort_bg_windows`, `cohort_enrichment_log` and `cohort_support_minus_bg` compare support with that background. Total genome presence is never used as the background.

**Joint, complementary patterns.** A signature reaction counts as jointly covered when independent units with at least `jointMinEffective` total weight carry it next to the family. `cohort_joint_coverage` versus `cohort_mean_single_coverage` gives `cohort_complementarity`: partial contexts in different lineages add up. `cohort_complementary_missing` is the weighted share of signature reactions that complete peer genomes lack entirely but that other lineages carry next to the family.

**Permutation null.** Within each peer genome, observed member windows are replaced by random matched background windows (same stratum), clade and cluster weights are preserved, and the seed is deterministic (`seed` combined with family, reaction and holdout key). Output: `cohort_perm_observed`, `cohort_perm_null_mean`, `cohort_perm_p` (add-one, `(k+1)/(B+1)`), `cohort_perm_n`. It is omitted when no matched windows exist.

**Paralogs.** Membership is soft and nothing transfers to all paralogs: a gene's evidence comes from its own family assignments, alternatives are returned per family with their own features, `cohort_multicopy_fraction` and the `cohort_single_copy_*` features show how much recurrence is ambiguous.

**Guards.** Empty signature (`signature_empty`), no peers, no informative peers, no matched background and no permutation pool produce explicit flags and omitted features, never `NaN`. Missing features are absent from the map rather than zero.

## Roles (item 7)

`RoleResolver` turns protein-intrinsic evidence (supplied role scores, direct EC calls) into a distribution over catalytic, accessory, transporter, regulator, complex subunit and unknown, with an unknown pseudo-mass so no evidence means unknown with probability 1. Context, operon association and cohort recurrence are not arguments to it and cannot make a protein catalytic. Each `FeatureRecord` carries the distribution, the resolved role (argmax if ≥ 0.5, else unknown), `role_noncatalytic_mass` (the explicit null), and `alternatives`: other family assignments with their features and other reactions claimed for the same gene.

## Outputs (`CohortResult`)

`records`: one `FeatureRecord` per input (gene, reaction) claim, with `features`, `states`, role distribution, alternatives and provenance (engine version, baseline snapshot id, baseline provenance, holdout level and the genomes held out, peer count, family resolution, seed). `aggregates`: one per (family, reaction) over **all** genomes, labelled `holdout = none`; they describe the cohort and must not be used to score a member.

## Wiring into the typed workflow and learned features

Not done in this package. Required by the caller:

1. Build one `CohortInput` per cohort **from the pre-cohort snapshot** (Stage 0/evidence claims and EC calls before any context stage runs), and set `baselineSnapshotId` to that snapshot's identifier. Never pass confidences from a previous cohort/feedback pass.
2. Metadata the workflow must supply: a layout per genome keyed by FAA seqid; orthogroup (and subfamily where available) membership with weights; clade ids or a tree, and strain-cluster ids (ANI-based clusters are suitable); per-genome `AssemblyInfo` with `status`, and per-contig `length` and `closed`. Without them the engine emits explicit unknown states and unit weights fall back to one genome = one clade, which overstates independence; record `cohort_phylogeny_known_fraction` in the learned feature table.
3. Reaction knowledge: a `ReactionGraph` with currency set and EC bindings, optional reaction directions, curated interaction edges and an organism reaction set.
4. Call `new CohortInference(input, graph, config).run()` once per batch, after all genomes' baselines exist and before the learned scorer. Join `FeatureRecord.features` on `(genomeId, proteinId, reactionId)`; proteins use FAA seqids. Pass the maps to the learned model unchanged, treating an absent key as missing, and let calibration choose any "often" threshold. Use `holdout = STRAIN_CLUSTER` or `CLADE` consistently with the train/test split used for the learned model, and never feed `aggregates` as per-gene features.
5. Persist `CohortResult.provenance` and each record's `provenance` in the run audit.
6. A single-genome run produces `cohort_no_peers = 1`; the workflow must route that case to a cohort-free feature set rather than interpret it as negative evidence.
