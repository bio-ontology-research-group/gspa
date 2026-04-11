package gspa.integration.prior

import gspa.integration.IntegrationState
import gspa.integration.MetabolicGap
import gspa.integration.Prior
import gspa.ontology.PathwayGraph
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Boosts candidate claims in operons whose high-confidence members share
 * a pathway membership — the "same-operon → same-pathway" heuristic,
 * refined with a quantitative consensus strength.
 *
 * <p>Plan §A.4 of Phase 7. Phase 7 implements boost-only mode: this prior
 * only boosts <i>existing</i> candidate claims, it never creates new
 * claims for unannotated proteins. Phase 8's {@code DarkMatterSuggester}
 * is the one that turns operon + gap context into brand-new claims for
 * previously-uncovered functions.</p>
 *
 * <p>For each operon O, compute the pathway consensus by looking at GO
 * terms of its high-confidence members (posterior > MAP). For each
 * pathway P that at least two members are annotated with, measure the
 * consensus strength as the fraction of operon members that agree. Then,
 * for each weak candidate claim {@code (protein_in_O, F)} where F is in
 * the pathway's required set, add
 * {@code +alphaCtx * consensusStrength}.</p>
 *
 * <p>Special interaction: if the operon contains a gapseq-flagged gap
 * candidate <i>and</i> the candidate's function would fill that gap, the
 * boost is multiplied by {@code alphaGapCtx}.</p>
 */
class GenomicContextPrior implements Prior {

    private static final Logger log = LoggerFactory.getLogger(GenomicContextPrior)

    /** Per-match boost magnitude (log-odds). Learned by the benchmark. */
    double alphaCtx = 1.0d

    /** Extra multiplier when the candidate function would fill a gapseq gap. */
    double alphaGapCtx = 1.5d

    /** Cached per iteration: (proteinId, goTerm) → boost magnitude. */
    private Map<String, Double> boosts = new LinkedHashMap<>()

    @Override
    String name() { 'genomic_context' }

    @Override
    Set<String> inputs() { ['operons', 'pathwayDatabase', 'metabolicGaps', 'posteriors'] as Set }

    @Override
    void beginIteration(IntegrationState state) {
        boosts = new LinkedHashMap<>()
        if (state.operons == null || state.operons.isEmpty()) return
        if (state.pathwayDatabase == null) return

        // Pathway GO term sets, precomputed once.
        Map<String, Set<String>> pathwayTerms = new LinkedHashMap<>()
        for (Map.Entry<String, PathwayGraph> e : state.pathwayDatabase.pathways.entrySet()) {
            Set<String> req = e.value.getRequiredGoTerms() ?: new LinkedHashSet<String>()
            if (!req.isEmpty()) pathwayTerms[e.key] = req
        }
        if (pathwayTerms.isEmpty()) return

        // Gap-filling set: GO terms that would close a metabolic gap.
        Set<String> gapFillerTerms = new LinkedHashSet<>()
        for (MetabolicGap gap : state.metabolicGaps ?: []) {
            if (gap.goTerm) gapFillerTerms << gap.goTerm
        }

        // For each operon: determine pathway consensus and emit boosts.
        for (List<String> operon : state.operons) {
            if (operon == null || operon.size() < 2) continue

            // Count how many operon members have a posterior > MAP for each pathway.
            Map<String, Integer> pathwayCounts = new LinkedHashMap<>()
            Set<String> strongOperonGoTerms = new LinkedHashSet<>()
            for (String pid : operon) {
                for (String fkey : state.functionKeysForProtein(pid)) {
                    if (state.getProbability(fkey) <= IntegrationState.MAP_THRESHOLD) continue
                    String[] parts = IntegrationState.splitFunctionKey(fkey)
                    if (parts == null || parts[1] != 'GO') continue
                    strongOperonGoTerms << parts[2]
                }
            }

            // For each strong GO term in operon, find pathways containing it.
            for (Map.Entry<String, Set<String>> pe : pathwayTerms.entrySet()) {
                int members = 0
                for (String goId : strongOperonGoTerms) {
                    if (pe.value.contains(goId)) members++
                }
                if (members >= 1) pathwayCounts[pe.key] = members
            }
            if (pathwayCounts.isEmpty()) continue

            // Consensus strength = top pathway's member fraction.
            String topPathway = null
            int topCount = 0
            for (Map.Entry<String, Integer> c : pathwayCounts.entrySet()) {
                if (c.value > topCount) { topCount = c.value; topPathway = c.key }
            }
            double consensusStrength = (double) topCount / operon.size()
            if (consensusStrength <= 0.0) continue

            Set<String> topPathwayTerms = pathwayTerms[topPathway]

            // Emit boosts for every operon member's candidate GO claims
            // whose functionId is in the top pathway's required set.
            for (String pid : operon) {
                for (String fkey : state.functionKeysForProtein(pid)) {
                    String[] parts = IntegrationState.splitFunctionKey(fkey)
                    if (parts == null || parts[1] != 'GO') continue
                    String goId = parts[2]
                    if (!topPathwayTerms.contains(goId)) continue

                    double boost = alphaCtx * consensusStrength
                    // Gap interaction bonus: candidate would close a gapseq gap.
                    if (gapFillerTerms.contains(goId)) boost *= alphaGapCtx

                    String boostKey = fkey
                    boosts.merge(boostKey, boost, { a, b -> Math.max(a as double, b as double) })
                }
            }
        }

        log.debug("GenomicContextPrior: ${boosts.size()} claims boosted across ${state.operons.size()} operons")
    }

    @Override
    double logOddsBoost(String proteinId, String functionKey, IntegrationState state) {
        return boosts.getOrDefault(functionKey, 0.0d)
    }
}
