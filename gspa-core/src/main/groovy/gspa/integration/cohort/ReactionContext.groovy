package gspa.integration.cohort

import gspa.ontology.ReactionGraph
import groovy.transform.CompileStatic

/**
 * The reaction signature of a candidate reaction R: the reactions whose
 * catalysts are expected near (or interacting with) R's catalyst. Built from
 * the metabolic network (not named pathways), curated independent
 * interaction edges, reaction directions, and an optional organism-specific
 * plausible reaction set.
 */
@CompileStatic
class ReactionContext {

    enum Flow { UPSTREAM, DOWNSTREAM, BOTH, DISTAL, CURATED }

    static class Link {
        String reactionId
        int distance
        double weight
        Flow flow
        boolean curated
        /** true / false when an organism reaction set exists, null when unknown. */
        Boolean inOrganismSet
    }

    private final ReactionGraph graph
    private final CohortInput input
    private final CohortConfig config
    private final Map<String, Map<String, Link>> cache = new HashMap<>()

    ReactionContext(ReactionGraph graph, CohortInput input, CohortConfig config) {
        this.graph = graph
        this.input = input
        this.config = config
    }

    /** Signature links of {@code reactionId}; empty when it has no neighbours or edges. */
    Map<String, Link> links(String reactionId) {
        Map<String, Link> hit = cache.get(reactionId)
        if (hit != null) return hit
        Map<String, Link> out = new LinkedHashMap<>()
        if (graph != null && graph.reactions.containsKey(reactionId)) {
            Map<String, Double> bfs = graph.bfs(reactionId, config.maxHops, config.hopAlpha)
            Set<String> adjacent = graph.neighbors(reactionId)
            bfs.each { String r, Double w ->
                int d = 1
                double ww = w.doubleValue()
                // distance recovered from weight = alpha^d
                if (config.hopAlpha > 0d && config.hopAlpha < 1d) {
                    d = (int) Math.round(Math.log(ww) / Math.log(config.hopAlpha))
                }
                out.put(r, new Link(reactionId: r, distance: d, weight: ww,
                        flow: d == 1 ? flowOf(reactionId, r) : Flow.DISTAL, curated: false))
            }
        }
        for (CohortInput.InteractionEdge e : input.interactionEdges) {
            String other = e.reactionA == reactionId ? e.reactionB : (e.reactionB == reactionId ? e.reactionA : null)
            if (other == null || other == reactionId) continue
            Flow f = Flow.CURATED
            Link existing = out.get(other)
            if (existing == null) {
                out.put(other, new Link(reactionId: other, distance: 1, weight: config.hopAlpha, flow: f, curated: true))
            } else {
                existing.curated = true
            }
        }
        Set<String> org = input.organismReactions
        out.values().each { Link l -> l.inOrganismSet = org == null ? (Boolean) null : (Boolean) org.contains(l.reactionId) }
        cache.put(reactionId, out)
        out
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
