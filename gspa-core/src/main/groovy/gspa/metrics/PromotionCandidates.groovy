package gspa.metrics

import gspa.model.Annotation
import gspa.model.Protein
import gspa.model.Genome
import gspa.ontology.GoOntology

/**
 * Shared "best candidate to promote a missing term onto" search, used by the
 * completeness and coherence enforcers. A truly-missing term F has no protein
 * annotated with F or a descendant, so the available evidence is a near-ancestor
 * of F: the candidate is the protein whose best annotation among F's ancestors
 * (within {@code maxHops} hops up the GO DAG) has the highest score, decayed by
 * hop distance. Returns {@code [protein, score, basis]} or {@code null}.
 */
class PromotionCandidates {

    static Map best(Genome genome, String term, GoOntology go,
                    int maxHops, double hopDecay, double minScore,
                    Set<String> excludeProteins = [] as Set) {
        Map<String, Integer> anc = ancestorsWithinHops(term, go, maxHops)
        if (anc.isEmpty()) return null
        Protein bestP = null
        double bestScore = -1.0d
        String bestBasis = null
        genome.proteins.each { Protein p ->
            if (excludeProteins.contains(p.id)) return
            p.annotations.goAnnotations().each { Annotation a ->
                Integer h = anc[a.value]
                if (h != null) {
                    double s = a.score * Math.pow(hopDecay, h)
                    if (s > bestScore) {
                        bestScore = s
                        bestP = p
                        bestBasis = "${a.value}@${String.format(Locale.ROOT, '%.3f', a.score)} (${h} hop${h > 1 ? 's' : ''} up) on ${p.id}"
                    }
                }
            }
        }
        (bestP != null && bestScore >= minScore) ?
            [protein: bestP, score: bestScore, basis: bestBasis] : null
    }

    /** Ancestors of {@code f} within {@code maxHops} hops, mapped to hop distance. */
    static Map<String, Integer> ancestorsWithinHops(String f, GoOntology go, int maxHops) {
        Map<String, Integer> hops = [:]
        Set<String> frontier = [f] as Set
        for (int h = 1; h <= maxHops; h++) {
            Set<String> next = new HashSet<>()
            frontier.each { String t ->
                Set<String> parents = new HashSet<>()
                parents.addAll(go.getDirectParents(t))
                parents.addAll(go.getPartOfParents(t))
                parents.each { p -> if (!hops.containsKey(p)) { hops[p] = h; next.add(p) } }
            }
            frontier = next
            if (frontier.isEmpty()) break
        }
        hops
    }
}
