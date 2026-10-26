package gspa.integration.cohort

import gspa.ontology.ReactionGraph
import groovy.transform.CompileStatic

/**
 * The reaction signature of a candidate reaction R: the reactions whose
 * catalysts are expected near (or interacting with) R's catalyst. Built from
 * the metabolic network (not named pathways), curated independent
 * interaction edges with their orientation, and reaction directions.
 *
 * <p>The signature is structural and genome-independent, so it is cached per
 * reaction. Genome-specific knowledge (plausible reaction sets) is applied by
 * the caller through {@link CohortInput#organismSet}.</p>
 */
@CompileStatic
class ReactionContext {

    enum Flow { UPSTREAM, DOWNSTREAM, BOTH, DISTAL, CURATED }

    static class Link {
        String reactionId
        /** Real hop count (breadth-first levels), independent of the weight. */
        int distance
        /** alpha^distance. */
        double weight
        Flow flow
        boolean curated
    }

    private final ReactionGraph graph
    private final CohortInput input
    private final CohortConfig config
    private final Map<String, Map<String, Link>> cache = new HashMap<>()
    private Map<String, List<CohortInput.InteractionEdge>> edgesByReaction

    ReactionContext(ReactionGraph graph, CohortInput input, CohortConfig config) {
        this.graph = graph
        this.input = input
        this.config = config
    }

    private Map<String, List<CohortInput.InteractionEdge>> edgeIndex() {
        if (edgesByReaction == null) {
            edgesByReaction = new HashMap<>()
            for (CohortInput.InteractionEdge e : input.interactionEdges) {
                edgesByReaction.computeIfAbsent(e.reactionA, { new ArrayList<CohortInput.InteractionEdge>() }).add(e)
                if (e.reactionB != e.reactionA) {
                    edgesByReaction.computeIfAbsent(e.reactionB, { new ArrayList<CohortInput.InteractionEdge>() }).add(e)
                }
            }
        }
        edgesByReaction
    }

    /** Signature links of {@code reactionId}; empty when it has no neighbours or edges. */
    Map<String, Link> links(String reactionId) {
        Map<String, Link> hit = cache.get(reactionId)
        if (hit != null) return hit
        Map<String, Link> out = new LinkedHashMap<>()
        if (graph != null && graph.reactions.containsKey(reactionId)) {
            Set<String> seen = new HashSet<>()
            seen.add(reactionId)
            List<String> frontier = [reactionId]
            double weight = 1.0d
            for (int d = 1; d <= config.maxHops && !frontier.isEmpty(); d++) {
                weight *= config.hopAlpha
                List<String> next = new ArrayList<>()
                for (String f : frontier) {
                    for (String n : graph.neighbors(f)) {
                        if (!seen.add(n)) continue
                        next.add(n)
                        out.put(n, new Link(reactionId: n, distance: d, weight: weight,
                                flow: d == 1 ? flowOf(reactionId, n) : Flow.DISTAL, curated: false))
                    }
                }
                frontier = next
            }
        }
        List<CohortInput.InteractionEdge> edges = edgeIndex().get(reactionId)
        if (edges != null) {
            for (CohortInput.InteractionEdge e : edges) {
                String other = e.reactionA == reactionId ? e.reactionB : e.reactionA
                if (other == reactionId) continue
                // Orientation relative to the focal reaction: A -> B means A is upstream of B.
                Flow orient = !e.directed ? Flow.CURATED : (e.reactionB == reactionId ? Flow.UPSTREAM : Flow.DOWNSTREAM)
                Link existing = out.get(other)
                if (existing == null) {
                    out.put(other, new Link(reactionId: other, distance: 1, weight: config.hopAlpha, flow: orient, curated: true))
                } else {
                    existing.curated = true
                    existing.flow = merge(existing.flow, orient)
                    if (existing.distance > 1) { existing.distance = 1; existing.weight = config.hopAlpha }
                }
            }
        }
        cache.put(reactionId, out)
        out
    }

    private static Flow merge(Flow network, Flow curated) {
        if (curated == Flow.CURATED || network == curated) return network
        if (network == Flow.DISTAL || network == Flow.CURATED) return curated
        Flow.BOTH
    }

    /** Effective (substrates, products) after applying the reaction's direction. */
    private List<Set<String>> effective(String rxn) {
        ReactionGraph.ReactionSpec spec = graph.reactions.get(rxn)
        Set<String> subs = new LinkedHashSet<>(spec.substrates ?: (Set<String>) [])
        Set<String> prods = new LinkedHashSet<>(spec.products ?: (Set<String>) [])
        subs.removeAll(graph.currencyMetabolites)
        prods.removeAll(graph.currencyMetabolites)
        ReactionDirection d = directionOf(rxn)
        if (d == ReactionDirection.REVERSE) return [prods, subs]
        if (d == ReactionDirection.FORWARD) return [subs, prods]
        // REVERSIBLE / UNKNOWN: either direction may operate.
        Set<String> both = new LinkedHashSet<>(subs)
        both.addAll(prods)
        [both, both]
    }

    ReactionDirection directionOf(String rxn) {
        ReactionDirection d = input.reactionDirections.get(rxn)
        d == null ? ReactionDirection.UNKNOWN : d
    }

    private Flow flowOf(String focal, String neighbour) {
        List<Set<String>> f = effective(focal)
        List<Set<String>> n = effective(neighbour)
        // neighbour upstream: its products feed the focal's substrates.
        boolean up = !n[1].intersect(f[0]).isEmpty()
        boolean down = !n[0].intersect(f[1]).isEmpty()
        up && down ? Flow.BOTH : (up ? Flow.UPSTREAM : (down ? Flow.DOWNSTREAM : Flow.DISTAL))
    }
}
