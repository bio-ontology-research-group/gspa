package gspa.integration.prior

import gspa.integration.IntegrationState
import gspa.integration.MetabolicGap
import gspa.integration.cohort.CohortInput
import gspa.integration.cohort.ContextMode
import gspa.model.GenomeLayout
import gspa.model.ProteinLocus
import gspa.ontology.ReactionGraph

/**
 * Computes reaction-local context support {@code c(P', gap) in [0, 1]} for
 * candidate proteins of a metabolic gap under UNIVERSAL, PATHWAY, or ORGANISM
 * reaction context modes.
 */
class GapContextSupport {

    ContextMode mode = ContextMode.UNIVERSAL
    int windowBp = 10000
    double anchorThreshold = 0.3d
    double alpha = 0.5d
    int maxHops = 2
    List<CohortInput.PathwayGroup> pathwayGroups = []
    Set<String> organismReactions = [] as Set

    Map<String, Object> lastReport = new LinkedHashMap<>()

    GapContextSupport() {
        resetReport()
    }

    GapContextSupport(ContextMode mode,
                      int windowBp = 10000,
                      double anchorThreshold = 0.3d,
                      double alpha = 0.5d,
                      int maxHops = 2,
                      List<CohortInput.PathwayGroup> pathwayGroups = null,
                      Set<String> organismReactions = null) {
        this.mode = mode
        this.windowBp = windowBp
        this.anchorThreshold = anchorThreshold
        this.alpha = alpha
        this.maxHops = maxHops
        this.pathwayGroups = pathwayGroups ?: []
        this.organismReactions = organismReactions ?: ([] as Set<String>)
        resetReport()
    }

    void resetReport() {
        lastReport = new LinkedHashMap<String, Object>([
            mode                             : mode?.name()?.toLowerCase(Locale.ROOT),
            window_bp                        : windowBp,
            anchor_threshold                 : anchorThreshold,
            alpha                            : alpha,
            max_hops                         : maxHops,
            gaps_evaluated                   : 0,
            gaps_with_context                : 0,
            context_reactions_per_gap        : new ArrayList<Integer>(),
            anchored_context_reactions_per_gap: new ArrayList<Integer>(),
            candidates_with_c_gt_0           : 0,
            gaps                             : new ArrayList<Map<String, Object>>(),
        ])
    }

    static Set<String> readOrganismReactions(File file) {
        if (file == null || !file.isFile()) {
            throw new IllegalArgumentException("file does not exist: ${file}")
        }
        Set<String> rxns = new LinkedHashSet<>()
        file.eachLine('UTF-8') { String raw ->
            String line = raw.trim()
            if (!line || line.startsWith('#')) return
            int hash = line.indexOf('#')
            if (hash >= 0) line = line.substring(0, hash).trim()
            if (line) rxns << line
        }
        rxns
    }

    Map<String, Double> contextSupport(MetabolicGap gap,
                                       Collection<String> candidateProteinIds,
                                       IntegrationState state) {
        if (state?.genomeLayout == null) {
            throw new IllegalStateException("Gap prior context mode ${mode} requires a genome layout (--genome-layout)")
        }
        if (lastReport == null || lastReport.isEmpty()) {
            resetReport()
        }
        ReactionGraph graph = state.reactionGraph
        if (graph == null || gap == null) {
            recordGap(gap, 0, 0, candidateProteinIds?.size() ?: 0, 0)
            return Collections.emptyMap()
        }

        // 1. Focal graph reactions
        Set<String> focalReactions = new LinkedHashSet<>()
        if (gap.ecNumber) {
            for (String rawEc : gap.ecNumber.split('[/,;\\s]+')) {
                String ec = rawEc.trim().replaceFirst(/^EC:/, '')
                if (!ec) continue
                focalReactions.addAll(graph.reactionsForEc(ec))
                if (rawEc.trim() != ec) {
                    focalReactions.addAll(graph.reactionsForEc(rawEc.trim()))
                }
            }
        }
        if (gap.reactionId && graph.reactions.containsKey(gap.reactionId)) {
            focalReactions << gap.reactionId
        }
        if (focalReactions.isEmpty()) {
            recordGap(gap, 0, 0, candidateProteinIds?.size() ?: 0, 0)
            return Collections.emptyMap()
        }

        // 2. Context reactions with weights, per mode, union over focal reactions r0, excluding focal reactions
        Map<String, Double> contextWeights = new LinkedHashMap<>()
        for (String r0 : focalReactions) {
            Map<String, Double> reached = contextReactionsForFocal(r0, graph)
            for (Map.Entry<String, Double> e : reached.entrySet()) {
                if (focalReactions.contains(e.key)) continue
                Double cur = contextWeights.get(e.key)
                if (cur == null || e.value > cur) {
                    contextWeights.put(e.key, e.value)
                }
            }
        }
        if (contextWeights.isEmpty()) {
            recordGap(gap, 0, 0, candidateProteinIds?.size() ?: 0, 0)
            return Collections.emptyMap()
        }

        // 3. Gap's own GO terms (gap.goTerm and GO terms of gap's ECs)
        Set<String> gapOwnGoTerms = new LinkedHashSet<>()
        if (gap.goTerm) {
            for (String rawGo : gap.goTerm.split('[,;\\s]+')) {
                String go = rawGo.trim()
                if (go) gapOwnGoTerms << go
            }
        }
        if (gap.ecNumber) {
            for (String rawEc : gap.ecNumber.split('[/,;\\s]+')) {
                String ec = rawEc.trim().replaceFirst(/^EC:/, '')
                if (ec) gapOwnGoTerms.addAll(lookupGoTerms(ec, state))
            }
        }

        // 4. Functions and anchors of each context reaction r
        GenomeLayout layout = state.genomeLayout
        Map<String, Map<String, Double>> posteriorsByGo = indexPosteriorsByGo(state)
        Map<String, Double> anchoredWeights = new LinkedHashMap<>()
        Map<String, List<ProteinLocus>> anchorsByReaction = new LinkedHashMap<>()

        for (Map.Entry<String, Double> entry : contextWeights.entrySet()) {
            String r = entry.key
            Set<String> ecs = new LinkedHashSet<>(graph.ecsForReaction(r) ?: Collections.<String>emptySet())
            if (ecs.isEmpty()) {
                String fallbackEc = graph.reactions.get(r)?.ecNumber
                if (fallbackEc) {
                    for (String tok : fallbackEc.split('[/,;\\s]+')) {
                        String trimmed = tok.trim()
                        if (trimmed) ecs << trimmed
                    }
                }
            }
            Set<String> funcs = new LinkedHashSet<>()
            for (String ec : ecs) {
                String cleanEc = ec.trim().replaceFirst(/^EC:/, '')
                if (cleanEc) funcs.addAll(lookupGoTerms(cleanEc, state))
            }
            funcs.removeAll(gapOwnGoTerms)
            if (funcs.isEmpty()) continue

            Set<String> anchorIds = new LinkedHashSet<>()
            for (String term : funcs) {
                Map<String, Double> byProt = posteriorsByGo.get(term)
                if (byProt == null) continue
                for (Map.Entry<String, Double> pe : byProt.entrySet()) {
                    if (pe.value != null && sigmoid(pe.value) >= anchorThreshold && layout.get(pe.key) != null) {
                        anchorIds << pe.key
                    }
                }
            }
            if (!anchorIds.isEmpty()) {
                anchoredWeights.put(r, entry.value)
                anchorsByReaction.put(r, anchorIds.collect { layout.get(it) })
            }
        }

        if (anchoredWeights.isEmpty()) {
            recordGap(gap, contextWeights.size(), 0, candidateProteinIds?.size() ?: 0, 0)
            return Collections.emptyMap()
        }

        // 5. Compute c(P') for each candidate protein
        double totalWeight = 0.0d
        for (Double w : anchoredWeights.values()) {
            totalWeight += w
        }
        Map<String, Double> support = new LinkedHashMap<>()
        int positiveCandidates = 0
        for (String candId : (candidateProteinIds ?: Collections.<String>emptyList())) {
            ProteinLocus loc = layout.get(candId)
            if (loc == null || totalWeight <= 0.0d) {
                support.put(candId, 0.0d)
                continue
            }
            int candMid = loc.midpoint()
            double num = 0.0d
            for (Map.Entry<String, Double> re : anchoredWeights.entrySet()) {
                List<ProteinLocus> anchors = anchorsByReaction.get(re.key)
                boolean present = false
                for (ProteinLocus a : anchors) {
                    if (a.proteinId != candId && a.contig == loc.contig && Math.abs(a.midpoint() - candMid) <= windowBp) {
                        present = true
                        break
                    }
                }
                if (present) {
                    num += re.value
                }
            }
            double cVal = num / totalWeight
            support.put(candId, cVal)
            if (cVal > 0.0d) positiveCandidates++
        }

        recordGap(gap, contextWeights.size(), anchoredWeights.size(), support.size(), positiveCandidates)
        support
    }

    private Map<String, Double> contextReactionsForFocal(String r0, ReactionGraph graph) {
        if (mode == ContextMode.UNIVERSAL) {
            return graph.bfs(r0, maxHops, alpha)
        }
        if (mode == ContextMode.ORGANISM) {
            Map<String, Double> out = new LinkedHashMap<>()
            if (maxHops <= 0 || !graph.reactions.containsKey(r0)) return out
            Set<String> allowed = organismReactions ?: Collections.<String>emptySet()
            Set<String> seen = new HashSet<>()
            seen.add(r0)
            List<String> frontier = [r0]
            double weight = 1.0d
            for (int d = 1; d <= maxHops && !frontier.isEmpty(); d++) {
                weight *= alpha
                List<String> next = new ArrayList<>()
                for (String f : frontier) {
                    for (String n : graph.neighbors(f)) {
                        if (!allowed.contains(n)) continue
                        if (!seen.add(n)) continue
                        out.put(n, weight)
                        next.add(n)
                    }
                }
                frontier = next
            }
            return out
        }
        if (mode == ContextMode.PATHWAY) {
            Map<String, Double> out = new LinkedHashMap<>()
            List<CohortInput.PathwayGroup> groups = pathwayGroups ?: Collections.<CohortInput.PathwayGroup>emptyList()
            List<CohortInput.PathwayGroup> ownGroups = groups.findAll { it.reactionIds != null && it.reactionIds.contains(r0) }
            if (ownGroups.isEmpty()) return out
            Set<String> ownReactions = new HashSet<>()
            Set<String> ownPathways = new HashSet<>()
            for (CohortInput.PathwayGroup g : ownGroups) {
                ownPathways.add(g.pathwayId)
                ownReactions.addAll(g.reactionIds)
            }
            for (CohortInput.PathwayGroup g : groups) {
                if (!ownPathways.contains(g.pathwayId) || g.reactionIds == null) continue
                for (String r : g.reactionIds) {
                    if (!ownReactions.contains(r)) {
                        out.put(r, alpha)
                    }
                }
            }
            return out
        }
        Collections.emptyMap()
    }

    private static List<String> lookupGoTerms(String cleanEc, IntegrationState state) {
        if (!cleanEc) return Collections.emptyList()
        if (state.ecToGoAlternatives != null) {
            List<String> terms = state.ecToGoAlternatives.get(cleanEc)
            if (terms == null) terms = state.ecToGoAlternatives.get('EC:' + cleanEc)
            return terms != null ? terms.findAll { it } : Collections.emptyList()
        }
        String legacy = state.pathwayDatabase?.ec2go?.get(cleanEc) ?: state.pathwayDatabase?.ec2go?.get('EC:' + cleanEc)
        legacy ? [legacy] : Collections.emptyList()
    }

    private static Map<String, Map<String, Double>> indexPosteriorsByGo(IntegrationState state) {
        Map<String, Map<String, Double>> idx = new LinkedHashMap<>()
        if (state?.posteriorLogOdds == null) return idx
        for (Map.Entry<String, Double> e : state.posteriorLogOdds.entrySet()) {
            String[] parts = IntegrationState.splitFunctionKey(e.key)
            if (parts == null || parts[1] != 'GO' || e.value == null) continue
            idx.computeIfAbsent(parts[2], { new LinkedHashMap<String, Double>() })[parts[0]] = e.value
        }
        idx
    }

    private void recordGap(MetabolicGap gap, int nContext, int nAnchored, int nCandidates, int nPositive) {
        lastReport.mode = mode?.name()?.toLowerCase(Locale.ROOT)
        lastReport.window_bp = windowBp
        lastReport.anchor_threshold = anchorThreshold
        lastReport.alpha = alpha
        lastReport.max_hops = maxHops
        lastReport.gaps_evaluated = ((lastReport.gaps_evaluated as Integer) ?: 0) + 1
        if (nAnchored > 0) {
            lastReport.gaps_with_context = ((lastReport.gaps_with_context as Integer) ?: 0) + 1
        }
        (lastReport.context_reactions_per_gap as List<Integer>) << nContext
        (lastReport.anchored_context_reactions_per_gap as List<Integer>) << nAnchored
        lastReport.candidates_with_c_gt_0 = ((lastReport.candidates_with_c_gt_0 as Integer) ?: 0) + nPositive
        (lastReport.gaps as List<Map<String, Object>>) << [
            pathway_id                : gap?.pathwayId,
            reaction_id               : gap?.reactionId,
            ec_number                 : gap?.ecNumber,
            go_term                   : gap?.goTerm,
            context_reactions         : nContext,
            anchored_context_reactions: nAnchored,
            candidates_evaluated      : nCandidates,
            candidates_with_c_gt_0    : nPositive,
        ]
    }

    private static double sigmoid(double x) {
        if (x >= 500.0d) return 1.0d
        if (x <= -500.0d) return 0.0d
        1.0d / (1.0d + Math.exp(-x))
    }
}
