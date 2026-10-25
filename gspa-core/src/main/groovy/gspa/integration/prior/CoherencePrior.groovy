package gspa.integration.prior

import gspa.integration.ClaimKey
import gspa.ontology.ModuleAssessor
import gspa.ontology.ModuleDefinitions

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

    // ---- module-based pathway refinement (used when the state carries KEGG module definitions) ----
    /** 'module' (evidence-gated KEGG modules, pinned rescues) or 'legacy' (v2.1.0 triggered-pathway boosts). 'auto' = module when definitions are loaded. */
    String pathwayMode = 'auto'
    /** Posterior a rescued step enzyme is pinned to. */
    double rescueTarget = 0.60d
    /** A module is refined only when activated and at least this complete. */
    double refineMinCompleteness = 0.5d
    /** Steps rescued per module at most (highest-scoring first). */
    int maxRescuePerModule = 2
    /** Lowest candidate Stage 0 likelihood considered. */
    double minCandidate = 0.16d
    /** Candidates whose evidence-adjusted log-odds is below this are left alone (0 = adjusted odds at least even). */
    double minScore = 0.0d
    double wSupport = 0.5d
    double wFree = 0.5d
    double wCompeting = 1.0d
    double wOperon = 1.0d
    double wModule = 1.0d
    /** Rescues already logged (one record per claim per run). */
    private final Set<String> logged = new HashSet<>()
    private static final Set<String> ASPECT_ROOTS = ['GO:0008150', 'GO:0003674', 'GO:0005575'] as Set

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

        boolean moduleMode = state.moduleDefinitions != null && pathwayMode != 'legacy'
        safeguardStrong = null
        Map<String, Double> moduleBoosts = [:]
        if (moduleMode) {
            try {
                moduleBoosts = moduleRescues(state)
            } catch (Exception ex) {
                log.warn("CoherencePrior: module refinement failed: ${ex.message}")
            }
        }
        // --- Pathway coherence (legacy): for each pathway, find missing required terms ---
        if (state.pathwayDatabase != null && !moduleMode) {
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
                        // Module mode leaves has_part process rescues out: validated at 16-29 % precision (benchmark/pathway_coherence, rounds 3-4).
                        if (state.candidateRescueTier(fk, pVal) >= 1 && !moduleMode) {
                            double bVal = Math.max(1.80d * alphaCoh * w, 0.14d - state.getLikelihoodLogOdds(fk))
                            boosts[fk] = Math.max(boosts.getOrDefault(fk, 0.0d), bVal)
                        }
                    }
                }
            }
            for (Map.Entry<String, Double> ce : complexKeyWeights.entrySet()) {
                // in module mode a positive complex-partner boost needs the same candidate evidence as a module rescue
                // Module mode keeps the singleton penalty but not the partner boosts (same validation).
                if (moduleMode && ce.value > 0) continue
                boosts[ce.key] = boosts.getOrDefault(ce.key, 0.0d) + alphaCoh * ce.value
            }
            for (Map.Entry<String, Double> mb : moduleBoosts.entrySet()) {
                boosts[mb.key] = Math.max(boosts.getOrDefault(mb.key, 0.0d), mb.value)
            }
            targetedKeyBoosts = boosts
        }

        log.debug("CoherencePrior: ${processMissingTerms.size()} process-missing, ${pathwayMissingTerms.size()} pathway-missing, ${complexKeyWeights.size()} complex-adjusted keys")
    }

    /** Proteins' confident (Stage 0 >= 0.50) GO terms, built lazily per iteration for the safeguard. */
    private Map<String, Set<String>> safeguardStrong = null

    /**
     * Candidate safeguard shared by every coherence lift in module mode: the candidate's evidence-adjusted
     * log-odds (Stage 0 log-odds + wSupport ln n_sup + wFree for an otherwise unannotated protein - wCompeting
     * when the protein already carries a confident annotation of the same aspect outside the candidate's
     * lineage) must reach {@link #minScore}.
     */
    private boolean passesSafeguard(IntegrationState state, String fk) {
        String[] parts = IntegrationState.splitFunctionKey(fk)
        if (parts == null || parts[1] != 'GO') return false
        if (safeguardStrong == null) {
            safeguardStrong = [:].withDefault { [] as Set }
            Map<String, Double> src = (state.likelihoodLogOdds != null && !state.likelihoodLogOdds.isEmpty()) ? state.likelihoodLogOdds : state.posteriorLogOdds
            for (Map.Entry<String, Double> e : src.entrySet()) {
                if (e.value < 0.0d) continue
                String[] p = IntegrationState.splitFunctionKey(e.key)
                if (p != null && p[1] == 'GO') safeguardStrong[p[0]] << p[2]
            }
        }
        String pid = parts[0], g = parts[2]
        double sc = state.getLikelihoodLogOdds(fk) + wSupport * Math.log(Math.max(1, state.getClaimCount(fk)))
        boolean free = state.getMaxLikelihoodProbForProtein(pid) < IntegrationState.MAP_THRESHOLD
        if (free) sc += wFree
        else if (state.goOntology != null) {
            String asp = state.goOntology.getAspect(g)
            Set<String> anc = state.goOntology.getAncestors(g) ?: ([] as Set)
            Set<String> desc = state.goOntology.getDescendants(g) ?: ([] as Set)
            boolean competing = safeguardStrong[pid].any { String o -> o != g && !anc.contains(o) && !desc.contains(o) && !ASPECT_ROOTS.contains(o) && state.goOntology.getAspect(o) == asp }
            if (competing) sc -= wCompeting
        }
        sc >= minScore
    }

    /**
     * Module-based refinement. Modules are assessed on the Stage 0 likelihoods (so the targets do not
     * move between iterations); a module activated by its evidence and at least {@link #refineMinCompleteness}
     * complete may have up to {@link #maxRescuePerModule} steps rescued, one candidate per step: the
     * sub-threshold claim on a step enzyme with the best score, where the score is the Stage 0 log-odds,
     * corroboration, a bonus for an otherwise unannotated protein, a penalty for a protein already confidently
     * annotated to another function of the same aspect, a bonus when an operon partner carries a present
     * step of the same module, and the module's completeness. The chosen claim is pinned to
     * {@link #rescueTarget}; every rescue is logged with its evidence. Steps without any candidate stay
     * open: the module is never forced to completion.
     */
    private Map<String, Double> moduleRescues(IntegrationState state) {
        Map<String, Double> boosts = [:]
        ModuleDefinitions defs = state.moduleDefinitions
        Set<String> baseline = (state.goOntology != null) ? state.baselineAnnotatedGoTermsPropagated() : state.baselineAnnotatedGoTerms()
        // evidence: best sub-threshold likelihood per GO term; candidates per GO term; proteins per present term
        Map<String, Double> evidence = [:]
        Map<String, List<String>> candsByGo = [:].withDefault { [] }
        Map<String, Set<String>> presentProtsByGo = [:].withDefault { [] as Set }
        Map<String, Set<String>> strongByProt = [:].withDefault { [] as Set }
        Map<String, Double> src = (state.likelihoodLogOdds != null && !state.likelihoodLogOdds.isEmpty()) ? state.likelihoodLogOdds : state.posteriorLogOdds
        for (Map.Entry<String, Double> e : src.entrySet()) {
            String[] parts = IntegrationState.splitFunctionKey(e.key)
            if (parts == null || parts[1] != 'GO') continue
            double p = 1.0d / (1.0d + Math.exp(-e.value))
            if (p >= IntegrationState.MAP_THRESHOLD) {
                presentProtsByGo[parts[2]] << parts[0]
                strongByProt[parts[0]] << parts[2]
            } else if (p >= minCandidate) {
                candsByGo[parts[2]] << e.key
                if (p > evidence.getOrDefault(parts[2], 0.0d)) evidence[parts[2]] = p
            }
        }
        Set<String> forbidden = [] as Set
        try {
            if (state.satConsistencyChecker?.organismTaxon) forbidden = state.satConsistencyChecker.getForbiddenConstrainedTerms() ?: ([] as Set)
        } catch (Exception ignored) {}
        def assessor = new ModuleAssessor(defs)
        assessor.domain = gspa.integration.ModuleReport.domainOf(state)
        def assessment = assessor.assess(baseline, evidence, forbidden)
        Map<String, List<String>> mates = [:]
        for (List<String> op : (state.operons ?: [])) {
            for (String pid : op) mates[pid] = op.findAll { it != pid }
        }
        double target = Math.log(rescueTarget / (1.0d - rescueTarget))
        int nRescued = 0
        for (ModuleAssessor.ModuleResult m : assessment.modules.values()) {
            if (!m.activated || m.signature || Double.isNaN(m.completeness) || m.completeness < refineMinCompleteness || m.complete) continue
            // proteins carrying a present step of this module (operon context)
            Set<String> moduleProts = [] as Set
            for (ModuleAssessor.StepResult st : m.steps) for (String g : st.presentGo) moduleProts.addAll(presentProtsByGo[g])
            List<Map> options = []
            for (ModuleAssessor.StepResult st : m.steps) {
                if (!st.required || st.status == 'present' || st.status == 'unmappable' || st.status == 'forbidden') continue
                Set<String> terms = new LinkedHashSet<>(st.goTerms)
                if (state.goOntology != null) for (String g : st.goTerms) terms.addAll(state.goOntology.getDescendants(g) ?: [])
                String bestKey = null; double bestScore = Double.NEGATIVE_INFINITY; Map bestEv = null
                for (String g : terms) {
                    for (String fk : candsByGo[g]) {
                        String[] parts = IntegrationState.splitFunctionKey(fk)
                        String pid = parts[0]
                        double l0 = state.getLikelihoodLogOdds(fk)
                        double sc = l0 + wSupport * Math.log(Math.max(1, state.getClaimCount(fk)))
                        boolean free = state.getMaxLikelihoodProbForProtein(pid) < IntegrationState.MAP_THRESHOLD
                        if (free) sc += wFree
                        boolean competing = false
                        if (!free && state.goOntology != null) {
                            String asp = state.goOntology.getAspect(g)
                            Set<String> anc = state.goOntology.getAncestors(g) ?: ([] as Set)
                            competing = strongByProt[pid].any { String o -> !terms.contains(o) && !anc.contains(o) && !ASPECT_ROOTS.contains(o) && state.goOntology.getAspect(o) == asp }
                            if (competing) sc -= wCompeting
                        }
                        String partner = (mates[pid] ?: []).find { moduleProts.contains(it) }
                        if (partner != null) sc += wOperon
                        sc += wModule * m.completeness
                        if (sc > bestScore) {
                            bestScore = sc; bestKey = fk
                            bestEv = [module: m.id, module_name: m.name, step: st.index, step_label: st.label, module_completeness: m.completeness,
                                      activation_score: m.activationScore, n_sup: state.getClaimCount(fk), free: free, competing: competing, operon_partner: partner]
                        }
                    }
                }
                if (bestKey != null && bestScore >= minScore) options << [key: bestKey, score: bestScore, ev: bestEv, like: state.getLikelihoodProbability(bestKey)]
            }
            options.sort { a, b -> (b.score as double) <=> (a.score as double) }
            for (Map o : options.take(maxRescuePerModule)) {
                String fk = o.key as String
                double l0 = state.getLikelihoodLogOdds(fk)
                boosts[fk] = Math.max(boosts.getOrDefault(fk, 0.0d), Math.max(1.80d * alphaCoh, target - l0))
                Double prev = state.getPinnedFloor(ClaimKey.parse(fk))
                if (prev == null || prev < target) state.setPinnedFloor(ClaimKey.parse(fk), target)
                String id = "${m.id}|${fk}".toString()
                if (logged.add(id)) {
                    state.rescueLog << [prior: name(), mode: 'module', category: "${m.id}:step${o.ev.step}".toString(), functionKey: fk,
                                        likelihood: o.like, score: o.score, target: rescueTarget, evidence: o.ev]
                }
                nRescued++
            }
        }
        if (nRescued) log.debug("CoherencePrior: ${nRescued} module step rescues")
        boosts
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
