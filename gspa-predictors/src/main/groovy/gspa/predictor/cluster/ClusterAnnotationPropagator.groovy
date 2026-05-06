package gspa.predictor.cluster

import gspa.model.Annotation
import gspa.model.Genome
import gspa.model.Protein
import gspa.model.ProteinCluster
import gspa.model.ProteinClusterSet
import gspa.model.ProteinRef
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Copy annotations from cluster representatives to their members.
 *
 * <p>Called by {@link gspa.predictor.AnnotationPipeline} after each
 * per-protein predictor finishes running on the representative FASTA.
 * For each cluster, every annotation on the representative is copied to
 * each non-representative member, stamped with provenance metadata so
 * the propagation is auditable.</p>
 *
 * <p>Agnostic to the number of input genomes — works identically in
 * Part 1 (one genome) and Phase 11 (many genomes).</p>
 */
class ClusterAnnotationPropagator {

    private static final Logger log = LoggerFactory.getLogger(ClusterAnnotationPropagator)

    /**
     * Only propagate annotations whose {@code source} is in this set, if set.
     * Null = propagate all.
     */
    Set<String> sourcesToPropagate = null

    /**
     * Walk each cluster, copy representative annotations to members.
     *
     * @param clusters       cluster set with representative → members
     * @param genomes        genomes whose proteins are being annotated
     *                       (indexed by genome id for O(1) lookup)
     * @param sourcePredictor predictor whose annotations are being propagated;
     *                       used to filter (copy only annotations with that
     *                       source). Pass null to propagate every annotation
     *                       on the representative.
     * @return number of annotations copied
     */
    int propagate(ProteinClusterSet clusters, List<Genome> genomes, String sourcePredictor = null) {
        if (clusters == null || clusters.clusters.isEmpty() || genomes == null || genomes.isEmpty()) return 0
        Map<String, Map<String, Protein>> proteinIndex = indexGenomes(genomes)

        int copied = 0
        for (ProteinCluster c : clusters.clusters) {
            Protein rep = resolve(proteinIndex, c.representative)
            if (rep == null || rep.annotations == null) continue

            List<Annotation> source = new ArrayList<>(rep.annotations.annotations)
            if (sourcePredictor != null) {
                source = source.findAll { it.source == sourcePredictor }
            } else if (sourcesToPropagate != null) {
                source = source.findAll { it.source in sourcesToPropagate }
            }
            if (source.isEmpty()) continue

            for (ProteinRef memberRef : c.members) {
                if (memberRef == c.representative) continue
                Protein member = resolve(proteinIndex, memberRef)
                if (member == null) continue
                for (Annotation a : source) {
                    member.annotations.add(copyAnnotation(a, c.clusterId, c.representative.proteinId))
                    copied++
                }
            }
        }
        log.debug("Propagated ${copied} annotations across ${clusters.clusterCount()} clusters (source=${sourcePredictor ?: 'all'})")
        copied
    }

    private static Annotation copyAnnotation(Annotation a, String clusterId, String repId) {
        Map<String, Object> md = new LinkedHashMap<>(a.metadata ?: [:])
        md.propagatedFrom = repId
        md.clusterId = clusterId
        new Annotation(
            type: a.type,
            value: a.value,
            score: a.score,
            source: a.source,
            evidence: a.evidence,
            metadata: md,
            goAspect: a.goAspect,
            evidenceType: a.evidenceType,
        )
    }

    /** Build {genomeId → {proteinId → Protein}} lookup from the input genomes. */
    private static Map<String, Map<String, Protein>> indexGenomes(List<Genome> genomes) {
        Map<String, Map<String, Protein>> out = new LinkedHashMap<>()
        for (Genome g : genomes) {
            Map<String, Protein> byId = new LinkedHashMap<>()
            for (Protein p : g.proteins) byId[p.id] = p
            out[g.id] = byId
        }
        out
    }

    private static Protein resolve(Map<String, Map<String, Protein>> idx, ProteinRef ref) {
        if (ref == null) return null
        Map<String, Protein> m = idx[ref.genomeId]
        m == null ? null : m[ref.proteinId]
    }
}
