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

    /** Relative boost when another protein in the same operon carries the same heteromeric complex ('n'). */
    double complexOperonPartnerWeight = 0.95d

    @Override
    void beginIteration(IntegrationState state) {
        boosts = new LinkedHashMap<>()
        if (state.operons == null || state.operons.isEmpty()) return

        // Pathway GO term sets, precomputed once.
        Map<String, Set<String>> pathwayTerms = new LinkedHashMap<>()
        if (state.pathwayDatabase != null) {
            for (Map.Entry<String, PathwayGraph> e : state.pathwayDatabase.pathways.entrySet()) {
                Set<String> req = e.value.getRequiredGoTerms() ?: new LinkedHashSet<String>()
                if (!req.isEmpty()) pathwayTerms[e.key] = req
            }
        }

        // Gap-filling set: GO terms that would close a metabolic gap.
        Set<String> gapFillerTerms = new LinkedHashSet<>()
        for (MetabolicGap gap : state.metabolicGaps ?: []) {
            if (gap.goTerm) gapFillerTerms << gap.goTerm
        }

        // For each operon: determine pathway consensus and co-operonic complex partners.
        for (List<String> operon : state.operons) {
            if (operon == null || operon.size() < 2) continue

            // Track strong GO terms per protein in this operon.
            Map<String, Set<String>> strongByProt = new LinkedHashMap<>()
            for (String pid : operon) {
                Set<String> pTerms = new LinkedHashSet<>()
                for (String fkey : state.functionKeysForProtein(pid)) {
                    if (state.getProbability(fkey) <= IntegrationState.MAP_THRESHOLD) continue
                    String[] parts = IntegrationState.splitFunctionKey(fkey)
                    if (parts == null || parts[1] != 'GO') continue
                    pTerms << parts[2]
                }
                if (!pTerms.isEmpty()) strongByProt[pid] = pTerms
            }
            if (strongByProt.isEmpty()) continue

            // 1. Pathway consensus across all active pathways in this operon
            if (!pathwayTerms.isEmpty()) {
                for (Map.Entry<String, Set<String>> pe : pathwayTerms.entrySet()) {
                    Set<String> reqTerms = pe.value
                    int activeProts = 0
                    for (Set<String> pTerms : strongByProt.values()) {
                        if (pTerms.any { reqTerms.contains(it) }) activeProts++
                    }
                    if (activeProts < 1) continue
                    double consensusStrength = (double) activeProts / operon.size()
                    for (String pid : operon) {
                        for (String fkey : state.functionKeysForProtein(pid)) {
                            String[] parts = IntegrationState.splitFunctionKey(fkey)
                            if (parts == null || parts[1] != 'GO') continue
                            String goId = parts[2]
                            if (!reqTerms.contains(goId)) continue
                            double boost = alphaCtx * consensusStrength
                            if (gapFillerTerms.contains(goId)) boost *= alphaGapCtx
                            boosts.merge(fkey, boost, { a, b -> Math.max(a as double, b as double) })
                        }
                    }
                }
            }

            // 2. Co-operonic obligate heteromeric complex subunit boost
            if (state.complexClassification != null && !state.complexClassification.isEmpty()) {
                for (String pid : operon) {
                    for (String fkey : state.functionKeysForProtein(pid)) {
                        String[] parts = IntegrationState.splitFunctionKey(fkey)
                        if (parts == null || parts[1] != 'GO') continue
                        String goId = parts[2]
                        if (state.complexClassification[goId] != 'n') continue
                        int partnerCount = 0
                        for (Map.Entry<String, Set<String>> e : strongByProt.entrySet()) {
                            if (e.key != pid && e.value.contains(goId)) partnerCount++
                        }
                        if (partnerCount >= 1) {
                            double cBoost = alphaCtx * complexOperonPartnerWeight
                            boosts.merge(fkey, cBoost, { a, b -> (a as double) + (b as double) })
                        }
                    }
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
