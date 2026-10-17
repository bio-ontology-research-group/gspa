package gspa.integration.prior

import gspa.integration.IntegrationState
import gspa.integration.Prior
import gspa.metrics.Coherence
import gspa.metrics.ProcessCoherenceResult
import gspa.ontology.PathwayGraph
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Boosts candidate claims that would close a triggered-but-unsatisfied
 * process dependency, or complete a partially-covered pathway.
 *
 * <p>Plan §A.4 of Phase 7. Uses
 * {@link Coherence#evaluateProcessCoherence(Set)} to find the unsatisfied
 * {@code (C, F)} has_part pairs, then for every candidate claim
 * {@code (protein, F)} returns {@code +alphaCoh * (1 - fracAnnotated)} —
 * where {@code fracAnnotated} reflects how much of the triggering process
 * is already covered (closer to completion → bigger boost).</p>
 *
 * <p>Pathway coherence: for each triggered pathway with partial coverage,
 * boost candidates whose GO term is in the pathway's required set,
 * scaled by {@code (1 - pathwayCoverage)}.</p>
 *
 * <p>No-ops cleanly when either {@code goReasoner} or
 * {@code pathwayDatabase} is missing from state.</p>
 */
class CoherencePrior implements Prior {

    private static final Logger log = LoggerFactory.getLogger(CoherencePrior)

    /** Per-match boost magnitude (log-odds). Learned by the benchmark. */
    double alphaCoh = 1.0d

    /** Relative weight for boosting missing partner subunits of a 1-subunit obligate heteromeric complex. */
    double complexPartnerWeight = 0.85d

    /** Relative weight for penalizing an isolated singleton obligate heteromeric complex with no candidate partner. */
    double complexSingletonPenaltyWeight = 0.75d

    /** Minimum sub-threshold probability for a protein to qualify as a candidate partner subunit. */
    double complexCandidateMinProb = 0.08d

    /**
     * Cached per iteration: unsatisfied F terms that would complete a
     * process dependency, each paired with the triggering C term and the
     * fraction already satisfied (so the boost can be weighted by
     * (1 - f)).
     */
    private Map<String, Double> processMissingTerms = new LinkedHashMap<>()

    /** Cached per iteration: pathway GO terms needed (term → (1 - coverage)). */
    private Map<String, Double> pathwayMissingTerms = new LinkedHashMap<>()

    /** Cached per iteration: functionKey -> complex coherence weight (+boost or -penalty). */
    private Map<String, Double> complexKeyWeights = new LinkedHashMap<>()

    /** Cached per iteration: functionKey -> targeted coherence log-odds adjustment when claimCounts is populated. */
    private Map<String, Double> targetedKeyBoosts = Collections.emptyMap()

    @Override
    String name() { 'coherence' }

    @Override
    Set<String> inputs() { ['goReasoner', 'pathwayDatabase', 'posteriors'] as Set }

    @Override
    void beginIteration(IntegrationState state) {
        processMissingTerms = new LinkedHashMap<>()
        pathwayMissingTerms = new LinkedHashMap<>()
        complexKeyWeights = new LinkedHashMap<>()
        targetedKeyBoosts = Collections.emptyMap()

        Set<String> annotated = (state.goOntology != null || state.goReasoner != null)
            ? state.currentlyAnnotatedGoTermsPropagated()
            : state.currentlyAnnotatedGoTerms()
        Set<String> baselineAnnotated = (state.goOntology != null || state.goReasoner != null)
            ? state.baselineAnnotatedGoTermsPropagated()
            : state.baselineAnnotatedGoTerms()

        Map<String, Double> missingRootWeights = new LinkedHashMap<>()
        Set<String> allPathwayDomain = new LinkedHashSet<>()

        // --- Process coherence: uses GoReasoner if available, or GoOntology has_part axioms ---
        if (state.goReasoner != null || state.goOntology != null) try {
            def coherence = new Coherence(state.goOntology, state.goReasoner)
            ProcessCoherenceResult result = coherence.evaluateProcessCoherence(baselineAnnotated)
            int triggered = result.triggered
            int satisfied = result.satisfied
            double fracAnnotated = triggered == 0 ? 1.0 : (satisfied / (double) triggered)
            double weight = Math.max(0.25d, fracAnnotated)

            for (Map.Entry<String, String> unsat : result.unsatisfied) {
                String fTerm = unsat.value
                if (fTerm != null) {
                    missingRootWeights.merge(fTerm, weight, { a, b -> Math.max(a as double, b as double) })
                    processMissingTerms.merge(fTerm, weight, { a, b -> Math.max(a as double, b as double) })
                    allPathwayDomain << fTerm
                    if (state.goOntology != null) {
                        for (String desc : state.goOntology.getDescendants(fTerm)) {
                            processMissingTerms.merge(desc, weight, { a, b -> Math.max(a as double, b as double) })
                        }
                    }
                }
            }
        } catch (Exception ex) {
            log.warn("CoherencePrior: process-coherence check failed: ${ex.message}")
        }

        // --- Pathway coherence: for each pathway, find missing required terms ---
        if (state.pathwayDatabase != null) {
            try {
                for (PathwayGraph pathway : state.pathwayDatabase.pathways.values()) {
                    Set<String> required = pathway.getRequiredGoTerms()
                    if (required == null || required.isEmpty()) continue
                    allPathwayDomain.addAll(required)
                    boolean triggeredPw = required.any { annotated.contains(it) || baselineAnnotated.contains(it) }
                    if (!triggeredPw) continue
                    int requiredCount = required.size()
                    int annotatedCount = required.count { baselineAnnotated.contains(it) }
                    double coverage = (double) annotatedCount / requiredCount
                    double weight = Math.max(0.35d, 0.35d + 0.65d * coverage)
                    for (String term : required) {
                        if (baselineAnnotated.contains(term)) continue
                        if (coverage >= 0.15d) {
                            missingRootWeights.merge(term, weight, { a, b -> Math.max(a as double, b as double) })
                        }
                        pathwayMissingTerms.merge(term, weight, { a, b -> Math.max(a as double, b as double) })
                        if (state.goOntology != null) {
                            for (String desc : state.goOntology.getDescendants(term)) {
                                pathwayMissingTerms.merge(desc, weight, { a, b -> Math.max(a as double, b as double) })
                            }
                        }
                    }
                }
            } catch (Exception ex) {
                log.warn("CoherencePrior: pathway check failed: ${ex.message}")
            }
        }

        // --- Complex coherence: obligate heteromeric complexes ('n') require >= 2 distinct proteins ---
        if (state.complexClassification != null && !state.complexClassification.isEmpty()) {
            try {
                Map<String, List<String>> strongByComplex = [:].withDefault { [] }
                Map<String, List<String>> candByComplex = [:].withDefault { [] }

                for (Map.Entry<String, Double> entry : state.posteriorLogOdds.entrySet()) {
                    String[] parts = IntegrationState.splitFunctionKey(entry.key)
                    if (parts == null || parts[1] != 'GO') continue
                    String goId = parts[2]
                    if (state.complexClassification[goId] != 'n') continue
                    double likeProb = state.getLikelihoodProbability(entry.key)
                    double postProb = state.getProbability(entry.key)
                    if (likeProb > IntegrationState.MAP_THRESHOLD) {
                        strongByComplex[goId] << entry.key
                    } else if (likeProb >= complexCandidateMinProb || postProb >= complexCandidateMinProb) {
                        candByComplex[goId] << entry.key
                    }
                }

                for (Map.Entry<String, List<String>> e : strongByComplex.entrySet()) {
                    String goId = e.key
                    List<String> strongKeys = e.value
                    List<String> candKeys = candByComplex[goId] ?: []
                    if (strongKeys.size() == 1) {
                        if (!candKeys.isEmpty()) {
                            for (String ck : candKeys) {
                                complexKeyWeights[ck] = complexPartnerWeight
                            }
                        } else {
                            complexKeyWeights[strongKeys[0]] = -complexSingletonPenaltyWeight
                        }
                    } else if (strongKeys.size() >= 2 && !candKeys.isEmpty()) {
                        for (String ck : candKeys) {
                            complexKeyWeights[ck] = 0.65d * complexPartnerWeight
                        }
                    }
                }
            } catch (Exception ex) {
                log.warn("CoherencePrior: complex check failed: ${ex.message}")
            }
        }

        if (state.claimCounts != null && !state.claimCounts.isEmpty()) {
            Map<String, Double> boosts = state.computeProteinSpecificityPrior(null, allPathwayDomain, 0.38d * alphaCoh, 0.42d * alphaCoh)
            if (!missingRootWeights.isEmpty()) {
                Map<String, List<Map.Entry<String, Double>>> candsByGo = [:].withDefault { [] }
                for (String fkey : state.posteriorLogOdds.keySet()) {
                    String[] parts = IntegrationState.splitFunctionKey(fkey)
                    if (parts == null || parts[1] != 'GO') continue
                    double likeProb = state.getLikelihoodProbability(fkey)
                    if (likeProb >= 0.18d && likeProb <= IntegrationState.MAP_THRESHOLD) {
                        candsByGo[parts[2]] << new AbstractMap.SimpleImmutableEntry<>(fkey, likeProb)
                    }
                }
                for (Map.Entry<String, Double> me : missingRootWeights.entrySet()) {
                    String rootTerm = me.key
                    double w = me.value
                    Set<String> tSet = new LinkedHashSet<>()
                    tSet << rootTerm
                    if (state.goOntology != null) {
                        tSet.addAll(state.goOntology.getDescendants(rootTerm))
                    }
                    List<Map.Entry<String, Double>> cands = []
                    for (String gt : tSet) {
                        List<Map.Entry<String, Double>> list = candsByGo.get(gt)
                        if (list != null) cands.addAll(list)
                    }
                    if (!cands.isEmpty()) {
                        state.sortCandidatesForRescue(cands)
                        String fk = cands[0].key
                        double pVal = cands[0].value as double
                        if (state.candidateRescueTier(fk, pVal) >= 1) {
                            double bVal = Math.max(1.80d * alphaCoh * w, 0.14d - state.getLikelihoodLogOdds(fk))
                            boosts[fk] = Math.max(boosts.getOrDefault(fk, 0.0d), bVal)
                        }
                    }
                }
            }
            for (Map.Entry<String, Double> ce : complexKeyWeights.entrySet()) {
                boosts[ce.key] = boosts.getOrDefault(ce.key, 0.0d) + alphaCoh * ce.value
            }
            targetedKeyBoosts = boosts
        }

        log.debug("CoherencePrior: ${processMissingTerms.size()} process-missing, ${pathwayMissingTerms.size()} pathway-missing, ${complexKeyWeights.size()} complex-adjusted keys")
    }

    @Override
    double logOddsBoost(String proteinId, String functionKey, IntegrationState state) {
        if (state.claimCounts != null && !state.claimCounts.isEmpty()) {
            return targetedKeyBoosts.getOrDefault(functionKey, 0.0d)
        }
        String[] parts = IntegrationState.splitFunctionKey(functionKey)
        if (parts == null || parts[1] != 'GO') return 0.0d
        String goId = parts[2]

        double processWeight = processMissingTerms.getOrDefault(goId, 0.0d)
        double pathwayWeight = pathwayMissingTerms.getOrDefault(goId, 0.0d)
        double complexWeight = complexKeyWeights.getOrDefault(functionKey, 0.0d)
        double posWeight = Math.max(processWeight, pathwayWeight)
        double totalWeight = posWeight + complexWeight
        if (totalWeight == 0.0d) return 0.0d
        alphaCoh * totalWeight
    }
}
