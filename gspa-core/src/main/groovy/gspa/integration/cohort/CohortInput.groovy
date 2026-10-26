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
        /** Strain-cluster identifier (null = unknown, treated as a singleton cluster). */
        String strainClusterId
        /** Null = no assembly metadata (status UNKNOWN, window completeness UNKNOWN). */
        AssemblyInfo assembly
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

    /** Protein-intrinsic role evidence (sequence/domain derived), 0..1 per role. */
    static class RoleEvidence {
        String genomeId
        String proteinId
        Map<Role, Double> scores = new LinkedHashMap<>()
    }

    /** Independent curated interaction between two reactions' catalysts. */
    static class InteractionEdge {
        String reactionA
        String reactionB
        String source
        /** If true the edge points A -> B (A upstream of B). */
        boolean directed
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
    /** Reactions plausibly present in the organism (null = unknown, not "empty"). */
    Set<String> organismReactions
    /** Optional supplied tree in Newick form; leaves are genome ids. */
    String treeNewick

    CohortInput addGenome(GenomeEntry g) {
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

    List<FamilyAssignment> assignments(String genomeId, String proteinId) {
        Map<String, List<FamilyAssignment>> m = membership.get(genomeId)
        List<FamilyAssignment> l = m == null ? null : m.get(proteinId)
        l == null ? Collections.<FamilyAssignment>emptyList() : l
    }

    void validate() {
        if (!baselineSnapshotId) {
            throw new IllegalArgumentException('baselineSnapshotId is required: baseline confidences must come from a named pre-cohort snapshot')
        }
        for (CandidateClaim c : claims) {
            if (c.baselineConfidence < 0d || c.baselineConfidence > 1d || Double.isNaN(c.baselineConfidence)) {
                throw new IllegalArgumentException("baselineConfidence out of [0,1] for ${c.genomeId}/${c.proteinId}/${c.reactionId}")
            }
            String prov = (c.provenance ?: '').toLowerCase()
            if (prov.startsWith('cohort') || prov.contains('feedback')) {
                throw new IllegalArgumentException("Claim ${c.genomeId}/${c.proteinId}/${c.reactionId} has provenance '${c.provenance}': cohort/feedback output cannot be a baseline")
            }
            if (!genomes.containsKey(c.genomeId)) {
                throw new IllegalArgumentException("Claim refers to unknown genome ${c.genomeId}")
            }
        }
        membership.each { String g, Map<String, List<FamilyAssignment>> m ->
            if (!genomes.containsKey(g)) throw new IllegalArgumentException("Membership refers to unknown genome ${g}")
            m.each { String p, List<FamilyAssignment> l ->
                double total = 0d
                l.each { total += it.weight }
                if (total > 1.0d + 1e-6d) {
                    throw new IllegalArgumentException("Membership weights of ${g}/${p} sum to ${total} > 1")
                }
            }
        }
    }
}
