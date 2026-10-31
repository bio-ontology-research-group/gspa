package gspa.integration.cohort

import gspa.ontology.ReactionGraph
import groovy.transform.CompileStatic

/**
 * The reaction signature of a candidate reaction R: the reactions whose
 * catalysts are expected near (or interacting with) R's catalyst. Built from
 * a universal or organism-restricted metabolic network, or named pathway groups.
 * Independent interaction edges and orientation apply to network modes only.
 *
 * <p>Organism mode restricts traversal, including intermediate nodes, and caches
 * by genome and reaction. Other modes share a structural signature across genomes.</p>
 */
@CompileStatic
class ReactionContext {

    enum Flow { UPSTREAM, DOWNSTREAM, BOTH, DISTAL, CURATED, PATHWAY }

    static class Link {
        String reactionId
        /** Correlated named-pathway mappings; null means this reaction alone. */
        Set<String> alternativeReactionIds
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
    private final Map<String, Map<String, Map<String, Link>>> organismCache = new HashMap<>()
    private Map<String, List<CohortInput.InteractionEdge>> edgesByReaction
    private Map<String, List<CohortInput.PathwayGroup>> groupsByReaction
    private Map<String, List<CohortInput.PathwayGroup>> groupsByPathway
    private final Map<Map<String, Link>, Map<String, Link>> aliasIndexes = new IdentityHashMap<>()

    ReactionContext(ReactionGraph graph, CohortInput input, CohortConfig config) {
        this.graph = graph
        this.input = input
        this.config = config
        if (config.contextMode == ContextMode.PATHWAY && input.pathwayGroups == null) {
            throw new IllegalArgumentException('PATHWAY context requires supplied pathway groups (empty is distinct from unavailable)')
        }
        if (config.contextMode == ContextMode.PATHWAY && !input.interactionEdges.isEmpty()) {
            throw new IllegalArgumentException('PATHWAY context is pathway-only; curated interactions require UNIVERSAL or ORGANISM context')
        }
        if (config.contextMode == ContextMode.ORGANISM) {
            for (String g : input.genomes.keySet()) {
                if (input.organismSet(g) == null) throw new IllegalArgumentException("ORGANISM context requires a known plausible reaction set for ${g}")
                for (String r : input.organismSet(g)) {
                    if (graph == null || !graph.reactions.containsKey(r)) throw new IllegalArgumentException("ORGANISM context for ${g} names unknown graph reaction ${r}")
                }
            }
        }
        for (CohortInput.PathwayGroup group : input.pathwayGroups ?: Collections.<CohortInput.PathwayGroup>emptyList()) {
            for (String r : group.reactionIds) {
                if (graph == null || !graph.reactions.containsKey(r)) {
                    throw new IllegalArgumentException("Pathway group ${group.pathwayId}/${group.groupId} names unknown graph reaction ${r}")
                }
            }
        }
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
        if (config.contextMode == ContextMode.ORGANISM) {
            throw new IllegalArgumentException('ORGANISM reaction context requires a genome id')
        }
        links(reactionId, null)
    }

    /** Each organism uses its own supplied network; the focal reaction remains a hypothesis. */
    Map<String, Link> links(String reactionId, String genomeId) {
        Set<String> allowed = null
        Map<String, Map<String, Link>> selectedCache = cache
        if (config.contextMode == ContextMode.ORGANISM) {
            if (genomeId == null || input.organismSet(genomeId) == null) {
                throw new IllegalArgumentException("ORGANISM context requires a known plausible reaction set for ${genomeId}")
            }
            allowed = input.organismSet(genomeId)
            selectedCache = organismCache.computeIfAbsent(genomeId, { new HashMap<String, Map<String, Link>>() })
        }
        Map<String, Link> hit = selectedCache.get(reactionId)
        if (hit != null) return hit
        Map<String, Link> out = new LinkedHashMap<>()
        if (config.contextMode == ContextMode.PATHWAY) {
            out.putAll(pathwayLinks(reactionId))
        } else if (graph != null && graph.reactions.containsKey(reactionId)) {
            Set<String> seen = new HashSet<>()
            seen.add(reactionId)
            List<String> frontier = [reactionId]
            double weight = 1.0d
            for (int d = 1; d <= config.maxHops && !frontier.isEmpty(); d++) {
                weight *= config.hopAlpha
                List<String> next = new ArrayList<>()
                for (String f : frontier) {
                    for (String n : graph.neighbors(f)) {
                        // Do not traverse an excluded intermediate and then filter only endpoints.
                        if (allowed != null && !allowed.contains(n)) continue
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
        Map<String, Link> frozen = Collections.unmodifiableMap(out)
        selectedCache.put(reactionId, frozen)
        frozen
    }

    /** Co-membership is not a biochemical direction or independent interaction. */
    private Map<String, Link> pathwayLinks(String focal) {
        if (groupsByReaction == null) {
            groupsByReaction = new HashMap<>()
            groupsByPathway = new HashMap<>()
            for (CohortInput.PathwayGroup group : input.pathwayGroups) {
                groupsByPathway.computeIfAbsent(group.pathwayId, { new ArrayList<CohortInput.PathwayGroup>() }).add(group)
                for (String r : group.reactionIds) {
                    groupsByReaction.computeIfAbsent(r, { new ArrayList<CohortInput.PathwayGroup>() }).add(group)
                }
            }
        }
        List<CohortInput.PathwayGroup> own = groupsByReaction.get(focal)
        if (own == null) return Collections.emptyMap()
        Set<String> alternatives = new HashSet<>()
        Map<String, String> parents = new HashMap<>()
        for (CohortInput.PathwayGroup group : own) alternatives.addAll(group.reactionIds)
        for (CohortInput.PathwayGroup group : own) {
            for (CohortInput.PathwayGroup other : groupsByPathway.get(group.pathwayId)) {
                List<String> members = new ArrayList<>()
                members.addAll(other.reactionIds)
                for (String r : members) parents.putIfAbsent(r, r)
                for (int i = 1; i < members.size(); i++) {
                    String a = rootOf(parents, members[0]), b = rootOf(parents, members[i])
                    if (a != b) parents.put(a < b ? b : a, a < b ? a : b)
                }
            }
        }
        Map<String, Set<String>> components = new TreeMap<>()
        for (String r : parents.keySet()) components.computeIfAbsent(rootOf(parents, r), { new TreeSet<String>() }).add(r)
        Map<String, Link> result = new LinkedHashMap<>()
        components.each { String key, Set<String> members ->
            if (Collections.disjoint(members, alternatives)) {
                result.put(key, new Link(reactionId: key, alternativeReactionIds: Collections.unmodifiableSet(members),
                        distance: 1, weight: config.hopAlpha, flow: Flow.PATHWAY, curated: false))
            }
        }
        result
    }

    private static String rootOf(Map<String, String> parents, String node) {
        String root = node
        while (parents.get(root) != root) root = parents.get(root)
        while (node != root) { String next = parents.get(node); parents.put(node, root); node = next }
        root
    }

    static Set<String> alternativesOf(Link link) {
        link.alternativeReactionIds ?: Collections.singleton(link.reactionId)
    }

    /** Resolve every retained alternative to its shared context-evidence group. */
    Link findLink(Map<String, Link> links, String reaction) {
        if (config.contextMode != ContextMode.PATHWAY) return links.get(reaction)
        Map<String, Link> aliases = aliasIndexes.get(links)
        if (aliases == null) {
            aliases = new HashMap<>()
            for (Link link : links.values()) for (String r : alternativesOf(link)) aliases.put(r, link)
            aliasIndexes.put(links, aliases)
        }
        aliases.get(reaction)
    }

    private static Flow merge(Flow network, Flow curated) {
        if (curated == Flow.CURATED || network == curated) return network
        if (network == Flow.DISTAL || network == Flow.CURATED || network == Flow.PATHWAY) return curated
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
