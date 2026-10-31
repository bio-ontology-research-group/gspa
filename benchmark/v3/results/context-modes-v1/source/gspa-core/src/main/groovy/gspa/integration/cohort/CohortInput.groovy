package gspa.integration.cohort

import gspa.model.GenomeLayout
import groovy.transform.CompileStatic

/**
 * Everything cohort inference reads. There are deliberately no label fields:
 * the input carries sequence/evidence-derived candidates, layouts, family
 * membership, phylogeny, assembly metadata and independent interaction
 * edges only. All protein identifiers are the FAA sequence ids that the
 * layout uses (not GFF ID attributes).
 *
 * <p>Baseline confidences must come from the pre-cohort sequence/evidence
 * snapshot named by {@link #baselineSnapshotId}; claims whose provenance
 * marks them as cohort/feedback output are rejected by {@link #validate}.</p>
 */
@CompileStatic
class CohortInput {

    static class ContigInfo {
        /** Contig length in bp; null = unknown. */
        Integer length
        /** True = closed replicon / complete contig, false = known fragment, null = unknown. */
        Boolean closed
    }

    static class AssemblyInfo {
        AssemblyStatus status = AssemblyStatus.UNKNOWN
        /** Optional completeness estimate (informational; never promotes UNKNOWN to COMPLETE). */
        Double completenessFraction
        Map<String, ContigInfo> contigs = new LinkedHashMap<>()
    }

    static class GenomeEntry {
        String genomeId
        GenomeLayout layout
        /** Clade identifier (null = unknown). */
        String cladeId
        /** Strain-cluster identifier (null = unknown). */
        String strainClusterId
        /** Null = no assembly metadata (status UNKNOWN, window completeness UNKNOWN). */
        AssemblyInfo assembly
        /**
         * Certification that the orthogroup/subfamily search was run over ALL
         * proteins of this genome (not only the input candidates), so that a
         * family without a member here is truly absent from the annotated
         * proteome. null/false = not certified: missing members are unknown.
         */
        Boolean familySearchComplete
        /** Optional fraction of proteins searched; must be 1.0 when familySearchComplete is true. */
        Double familySearchCoverage
        /**
         * Proteins explicitly declared as known but without a usable coordinate in {@link #layout} (absent from the
         * layout, split across contigs, multi-locus). They are real proteins of the genome: a family member or EC
         * claim on one is accepted, counts as family PRESENT, and its local context is UNKNOWN (no window, no
         * anchor). No coordinate is invented for them. An id that is neither in the layout nor declared here is
         * still rejected as a probable typo; an id in both is a contradiction and rejected too. Whether the family
         * search covered these proteins is {@link #familySearchComplete}, not this set.
         */
        Set<String> unplacedProteins = new LinkedHashSet<>()
    }

    /** Soft membership of a protein in an orthogroup (optionally a subfamily). */
    static class FamilyAssignment {
        String proteinId
        String orthogroupId
        String subfamilyId
        double weight = 1.0d

        /** Recurrence is resolved at subfamily level when one is given. */
        String familyKey() {
            subfamilyId ? orthogroupId + '#' + subfamilyId : orthogroupId
        }
    }

    /** A (gene, reaction) candidate from the pre-cohort snapshot. */
    static class CandidateClaim {
        String genomeId
        String proteinId
        String reactionId
        double baselineConfidence
        /** EC number the candidate came from (direct EC -> reaction support, no GO needed). */
        String ecNumber
        boolean ecDirect
        /** Number of reactions the EC maps to (specificity of the direct support). */
        int ecReactionCount = 1
        /** Source of the baseline confidence, e.g. "diamond+interpro@snapshot". */
        String provenance
    }

    /** A direct EC call for a protein, expanded to reactions by {@link EcCandidateExpander}. */
    static class EcCall {
        String genomeId
        String proteinId
        String ecNumber
        double confidence
        String provenance
    }

    /** Protein-intrinsic role evidence (sequence/domain derived), raw 0..1 per role, with its source. */
    static class RoleEvidence {
        String genomeId
        String proteinId
        String source
        Map<Role, Double> scores = new LinkedHashMap<>()
    }

    /** Independent curated interaction between two reactions' catalysts. */
    static class InteractionEdge {
        String reactionA
        String reactionB
        String source
        /** If true the edge points A -> B: A is upstream of B. */
        boolean directed
    }

    /** One named-pathway reaction group, whose reaction mappings are alternatives. */
    static class PathwayGroup {
        String pathwayId
        String groupId
        String source
        Set<String> reactionIds = new LinkedHashSet<>()
    }

    /** Identifier of the pre-cohort snapshot the baseline confidences were drawn from. */
    String baselineSnapshotId

    final Map<String, GenomeEntry> genomes = new LinkedHashMap<>()
    /** genome -> protein -> assignments. */
    final Map<String, Map<String, List<FamilyAssignment>>> membership = new LinkedHashMap<>()
    final List<CandidateClaim> claims = new ArrayList<>()
    final Map<String, RoleEvidence> roleEvidence = new LinkedHashMap<>()
    final List<InteractionEdge> interactionEdges = new ArrayList<>()
    final Map<String, ReactionDirection> reactionDirections = new LinkedHashMap<>()
    /** Null = unavailable; [] = supplied, known-empty pathway membership. */
    List<PathwayGroup> pathwayGroups = null
    /**
     * Per-genome plausible reaction sets. A genome absent from the map has an
     * UNKNOWN set; a genome mapped to an empty set has a known empty set.
     */
    final Map<String, Set<String>> organismReactions = new LinkedHashMap<>()
    /** Optional supplied tree in Newick form; leaves must be exactly the genome ids. */
    String treeNewick

    CohortInput addGenome(GenomeEntry g) {
        if (!g.genomeId) throw new IllegalArgumentException('genomeId is required')
        if (genomes.containsKey(g.genomeId)) throw new IllegalArgumentException("Duplicate genome ${g.genomeId}")
        genomes.put(g.genomeId, g)
        this
    }

    CohortInput assign(String genomeId, FamilyAssignment a) {
        membership.computeIfAbsent(genomeId, { new LinkedHashMap<String, List<FamilyAssignment>>() })
                .computeIfAbsent(a.proteinId, { new ArrayList<FamilyAssignment>() }).add(a)
        this
    }

    CohortInput addClaim(CandidateClaim c) {
        claims.add(c)
        this
    }

    CohortInput addRoleEvidence(RoleEvidence r) {
        roleEvidence.put(r.genomeId + '|' + r.proteinId, r)
        this
    }

    CohortInput addEdge(InteractionEdge e) {
        interactionEdges.add(e)
        this
    }

    /** Sets a genome's plausible reaction set; pass null to mark it unknown. */
    CohortInput setOrganismReactions(String genomeId, Collection<String> reactions) {
        if (reactions == null) organismReactions.remove(genomeId)
        else organismReactions.put(genomeId, new LinkedHashSet<String>(reactions))
        this
    }

    /** Plausible reaction set of a genome; null when unknown (distinct from empty). */
    Set<String> organismSet(String genomeId) {
        organismReactions.get(genomeId)
    }

    List<FamilyAssignment> assignments(String genomeId, String proteinId) {
        Map<String, List<FamilyAssignment>> m = membership.get(genomeId)
        List<FamilyAssignment> l = m == null ? null : m.get(proteinId)
        l == null ? Collections.<FamilyAssignment>emptyList() : l
    }

    private static boolean finite(double d) { !Double.isNaN(d) && !Double.isInfinite(d) }

    private boolean knownProtein(String genomeId, String proteinId) {
        GenomeEntry e = genomes.get(genomeId)
        e?.layout == null || e.layout.get(proteinId) != null || e.unplacedProteins.contains(proteinId)
    }

    /** Rejects inconsistent or malformed input before any computation. */
    void validate() {
        if (!baselineSnapshotId) {
            throw new IllegalArgumentException('baselineSnapshotId is required: baseline confidences must come from a named pre-cohort snapshot')
        }
        genomes.each { String id, GenomeEntry g ->
            if (g.genomeId != id) throw new IllegalArgumentException("Genome key ${id} != genomeId ${g.genomeId}")
            if (g.layout != null) {
                if (g.layout.genomeId && g.layout.genomeId != id) {
                    throw new IllegalArgumentException("Layout of ${id} carries genomeId ${g.layout.genomeId}")
                }
                g.layout.byProtein.each { String p, gspa.model.ProteinLocus l ->
                    if (l.genomeId && l.genomeId != id) throw new IllegalArgumentException("Locus ${p} of ${id} carries genomeId ${l.genomeId}")
                }
            }
            g.unplacedProteins.each { String p ->
                if (!p?.trim()) throw new IllegalArgumentException("Genome ${id}: blank id in unplacedProteins")
                if (g.layout?.get(p) != null) {
                    throw new IllegalArgumentException("Genome ${id}: protein ${p} is declared unplaced but has a coordinate in the layout")
                }
            }
            if (g.familySearchComplete == Boolean.TRUE && g.familySearchCoverage != null && g.familySearchCoverage < 1.0d) {
                throw new IllegalArgumentException("Genome ${id}: familySearchComplete=true contradicts familySearchCoverage=${g.familySearchCoverage}")
            }
            if (g.familySearchCoverage != null && !(g.familySearchCoverage >= 0d && g.familySearchCoverage <= 1d)) {
                throw new IllegalArgumentException("Genome ${id}: familySearchCoverage out of [0,1]")
            }
        }
        Set<String> seen = new HashSet<>()
        for (CandidateClaim c : claims) {
            if (!c.reactionId || !c.proteinId) throw new IllegalArgumentException('Claim without proteinId/reactionId')
            if (!finite(c.baselineConfidence) || c.baselineConfidence < 0d || c.baselineConfidence > 1d) {
                throw new IllegalArgumentException("baselineConfidence out of [0,1] for ${c.genomeId}/${c.proteinId}/${c.reactionId}")
            }
            String prov = (c.provenance ?: '').toLowerCase()
            if (prov.startsWith('cohort') || prov.contains('feedback')) {
                throw new IllegalArgumentException("Claim ${c.genomeId}/${c.proteinId}/${c.reactionId} has provenance '${c.provenance}': cohort/feedback output cannot be a baseline")
            }
            if (!genomes.containsKey(c.genomeId)) {
                throw new IllegalArgumentException("Claim refers to unknown genome ${c.genomeId}")
            }
            if (!knownProtein(c.genomeId, c.proteinId)) {
                throw new IllegalArgumentException("Claim protein ${c.proteinId} is not in the layout of ${c.genomeId} and not declared unplaced")
            }
            if (!seen.add(c.genomeId + '\t' + c.proteinId + '\t' + c.reactionId)) {
                throw new IllegalArgumentException("Duplicate claim ${c.genomeId}/${c.proteinId}/${c.reactionId}")
            }
        }
        membership.each { String g, Map<String, List<FamilyAssignment>> m ->
            if (!genomes.containsKey(g)) throw new IllegalArgumentException("Membership refers to unknown genome ${g}")
            m.each { String p, List<FamilyAssignment> l ->
                if (!knownProtein(g, p)) throw new IllegalArgumentException("Membership protein ${p} is not in the layout of ${g} and not declared unplaced")
                double total = 0d
                Set<String> keys = new HashSet<>()
                for (FamilyAssignment a : l) {
                    if (!a.orthogroupId) throw new IllegalArgumentException("Assignment of ${g}/${p} has no orthogroupId")
                    if (!finite(a.weight) || a.weight < 0d || a.weight > 1d) {
                        throw new IllegalArgumentException("Membership weight of ${g}/${p} must be finite and in [0,1]")
                    }
                    if (!keys.add(a.familyKey())) throw new IllegalArgumentException("Duplicate assignment ${a.familyKey()} for ${g}/${p}")
                    total += a.weight
                }
                if (total > 1.0d + 1e-9d) {
                    throw new IllegalArgumentException("Membership weights of ${g}/${p} sum to ${total} > 1")
                }
            }
        }
        roleEvidence.each { String k, RoleEvidence r ->
            if (!genomes.containsKey(r.genomeId)) throw new IllegalArgumentException("Role evidence refers to unknown genome ${r.genomeId}")
            if (!knownProtein(r.genomeId, r.proteinId)) {
                throw new IllegalArgumentException("Role evidence protein ${r.proteinId} is not in the layout of ${r.genomeId} and not declared unplaced")
            }
            r.scores.each { Role role, Double v ->
                if (v == null || !finite(v.doubleValue()) || v < 0d || v > 1d) {
                    throw new IllegalArgumentException("Role evidence ${role} of ${r.genomeId}/${r.proteinId} must be finite and in [0,1]")
                }
            }
        }
        for (InteractionEdge e : interactionEdges) {
            if (!e.reactionA || !e.reactionB) throw new IllegalArgumentException('Interaction edge needs both reactions')
        }
        Set<List<String>> pathwayKeys = new HashSet<>()
        for (PathwayGroup group : pathwayGroups ?: Collections.<PathwayGroup>emptyList()) {
            if (!group.pathwayId?.trim() || !group.groupId?.trim() || !group.source?.trim() ||
                !group.reactionIds || group.reactionIds.any { String r -> !r?.trim() }) {
                throw new IllegalArgumentException('Pathway group requires pathwayId, groupId, source and reaction alternatives')
            }
            if (!pathwayKeys.add([group.pathwayId, group.groupId])) {
                throw new IllegalArgumentException("Duplicate pathway group ${group.pathwayId}/${group.groupId}")
            }
        }
        organismReactions.keySet().each { String g ->
            if (!genomes.containsKey(g)) throw new IllegalArgumentException("Organism reaction set refers to unknown genome ${g}")
        }
    }
}
