package gspa.ontology

import groovy.transform.Canonical
import groovy.transform.CompileStatic

/**
 * Panel-level reaction graph: reactions connected via shared
 * non-currency metabolites. Built by {@link ReactionGraphLoader}
 * from gapsmith's seed_reactions.tsv + diffusion_mets.tsv +
 * degree-threshold currency filter.
 *
 * <p>Two lookup modes:</p>
 * <ul>
 *   <li>{@link #neighbors} — 1-hop neighbors of a reaction (cached).</li>
 *   <li>{@link #bfs} — reach up to distance {@code k} with
 *       distance-weighted values {@code alpha^d}.</li>
 * </ul>
 *
 * <p>The stoichiometry-based graph is directed implicitly
 * ({@code substrate → product}) but here we collapse to undirected
 * adjacency for the M1 suggester; directional variants are an M3
 * feature refinement.</p>
 */
@CompileStatic
class ReactionGraph {

    @Canonical
    static class ReactionSpec {
        String rxnId
        String ecNumber           // may be null / empty
        Set<String> substrates    // compound IDs (cpd*)
        Set<String> products
    }

    /** rxnId -> spec */
    final Map<String, ReactionSpec> reactions = new LinkedHashMap<>()

    /** metaboliteId -> rxnIds that touch it (after currency pruning) */
    final Map<String, Set<String>> metaboliteToReactions = new LinkedHashMap<>()

    /** currency metabolites (excluded from adjacency) */
    final Set<String> currencyMetabolites = new LinkedHashSet<>()

    /** EC number → reactions with that EC (populated by EC-aliases loader). */
    final Map<String, Set<String>> ecToReactions = new LinkedHashMap<>()

    /** rxnId → EC numbers (from ec_numbers column or EC-aliases file). */
    final Map<String, Set<String>> reactionToEcs = new LinkedHashMap<>()

    /** 1-hop neighbor cache */
    private final Map<String, Set<String>> neighborCache = new HashMap<>()

    int size() { reactions.size() }

    void addReaction(ReactionSpec spec) {
        reactions[spec.rxnId] = spec
    }

    void markCurrency(String metaboliteId) {
        currencyMetabolites << metaboliteId
    }

    /** Associate an EC number with a reaction (both-ways). */
    void bindEc(String rxnId, String ec) {
        if (!ec || !rxnId) return
        ecToReactions.computeIfAbsent(ec, { new LinkedHashSet<String>() }) << rxnId
        reactionToEcs.computeIfAbsent(rxnId, { new LinkedHashSet<String>() }) << ec
        // Keep ReactionSpec.ecNumber as the first-assigned value (back-compat).
        ReactionSpec spec = reactions[rxnId]
        if (spec != null && !spec.ecNumber) spec.ecNumber = ec
    }

    Set<String> reactionsForEc(String ec) {
        ecToReactions.getOrDefault(ec, Collections.<String>emptySet())
    }

    Set<String> ecsForReaction(String rxnId) {
        reactionToEcs.getOrDefault(rxnId, Collections.<String>emptySet())
    }

    /**
     * Call after all reactions + currency metabolites are loaded.
     * Builds metabolite → reaction indices excluding currency.
     */
    void build() {
        metaboliteToReactions.clear()
        for (ReactionSpec spec : reactions.values()) {
            for (String m : (Set<String>)(spec.substrates + spec.products)) {
                if (currencyMetabolites.contains(m)) continue
                metaboliteToReactions.computeIfAbsent(m, { new LinkedHashSet<String>() }) << spec.rxnId
            }
        }
    }

    /** 1-hop neighbors via any shared non-currency metabolite. */
    Set<String> neighbors(String rxnId) {
        Set<String> cached = neighborCache[rxnId]
        if (cached != null) return cached
        Set<String> out = new LinkedHashSet<>()
        ReactionSpec spec = reactions[rxnId]
        if (spec == null) {
            neighborCache[rxnId] = out
            return out
        }
        for (String m : (Set<String>)(spec.substrates + spec.products)) {
            if (currencyMetabolites.contains(m)) continue
            Set<String> others = metaboliteToReactions[m]
            if (others == null) continue
            for (String r : others) {
                if (r != rxnId) out << r
            }
        }
        neighborCache[rxnId] = out
        out
    }

    /**
     * BFS from {@code rxnId}: returns neighbor reaction ID → weight
     * {@code alpha^distance} for distances {@code 1..maxK}.
     * The source reaction itself is NOT included.
     */
    Map<String, Double> bfs(String rxnId, int maxK, double alpha) {
        Map<String, Double> out = new LinkedHashMap<>()
        if (maxK <= 0) return out
        Set<String> seen = new HashSet<>()
        seen << rxnId
        List<String> frontier = [rxnId]
        double weight = 1.0d
        for (int d = 1; d <= maxK; d++) {
            weight *= alpha
            List<String> next = []
            for (String f : frontier) {
                for (String n : neighbors(f)) {
                    if (seen.add(n)) {
                        out[n] = weight
                        next << n
                    }
                }
            }
            if (next.isEmpty()) break
            frontier = next
        }
        out
    }
}
