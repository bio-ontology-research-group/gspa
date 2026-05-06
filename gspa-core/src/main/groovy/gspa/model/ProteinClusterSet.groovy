package gspa.model

/**
 * A collection of protein clusters produced by a {@code ProteinClusterer}.
 *
 * Single-genome (Phase 10 Part 1) and cross-genome (Phase 11) use the same
 * class; the only difference is whether member {@link ProteinRef}s share
 * one genomeId. Lookups are indexed for constant-time protein → cluster
 * resolution in the propagation path.
 */
class ProteinClusterSet {

    List<ProteinCluster> clusters = []
    double identityThreshold
    double coverageThreshold

    // Lazily-built indexes keyed by ProteinRef.
    private Map<ProteinRef, ProteinCluster> byMember
    private Map<ProteinRef, ProteinCluster> byRep

    private void ensureIndexes() {
        if (byMember != null) return
        byMember = [:]
        byRep = [:]
        for (ProteinCluster c : clusters) {
            byRep[c.representative] = c
            for (ProteinRef m : c.members) byMember[m] = c
        }
    }

    /** Cluster containing the given protein, or null if unclustered. */
    ProteinCluster clusterFor(ProteinRef protein) {
        ensureIndexes()
        byMember[protein]
    }

    /** Representative for the given protein; returns the protein itself if unclustered. */
    ProteinRef representativeOf(ProteinRef protein) {
        ProteinCluster c = clusterFor(protein)
        c != null ? c.representative : protein
    }

    /** True if the protein IS the representative of its cluster. */
    boolean isRepresentative(ProteinRef protein) {
        ensureIndexes()
        byRep.containsKey(protein)
    }

    /** All representatives in this set. */
    List<ProteinRef> representatives() {
        clusters.collect { it.representative }
    }

    int clusterCount() { clusters.size() }

    /** Total number of proteins covered (sum over clusters). */
    int proteinCount() { clusters.sum { it.size() } ?: 0 }
}
