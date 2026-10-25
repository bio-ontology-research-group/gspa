package gspa.integration.prior

import gspa.integration.ClaimKey
import gspa.integration.IntegrationState
import gspa.integration.Prior
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Completes the genome's essential-function profile.
 *
 * <p>The essential function list comes from {@link gspa.config.EssentialFunctions};
 * a category is <i>uncovered</i> when no protein carries it, or a subclass of it,
 * at the Stage 0 likelihood above the decision threshold (propagated through the
 * GO ontology when one is wired). For every uncovered category the prior selects
 * the borderline candidate claim (likelihood in {@code [minCandidate, 0.50]}) that
 * is most likely to be the missing gene and rescues it. Two rescue modes exist:</p>
 *
 * <ul>
 *   <li>{@code targeted} (default). Candidates are scored on the evidence that
 *   distinguishes a missing essential gene from a spurious alignment: the Stage 0
 *   log-odds; corroboration across claims ({@code wSupport * ln n_sup}); whether
 *   the protein is otherwise unannotated, so the rescue gives it a first function
 *   rather than a second one ({@code wFree}); whether the protein already carries
 *   a confident annotation of the same aspect outside the category, which the
 *   rescue would contradict ({@code wCompeting}); and the fraction of its operon
 *   partners annotated to the essential domain, since essential genes cluster
 *   ({@code wContext}). The top {@code maxRescuePerCategory} candidates are pinned
 *   to the posterior {@code rescueTarget}: the refiner's floor mechanism keeps them
 *   there regardless of damping, iteration count and the genome-scale cap on prior
 *   lifts, so a rescued function is present at the threshold by construction
 *   instead of landing on either side of it.</li>
 *   <li>{@code legacy}. The v2.1.0 behaviour: tier ranking
 *   ({@link IntegrationState#candidateRescueTier}) and a boost of
 *   {@code max(1.30 alphaEss, 0.14 - L0)} on the top candidate, which the refiner
 *   caps at log-odds 0.12 and damps, so the rescued posterior ends between about
 *   0.47 and 0.53.</li>
 * </ul>
 *
 * <p>In both modes every claim on a term in the essential domain receives the
 * specificity-weighted boost of {@code 0.30 alphaEss} from
 * {@link IntegrationState#computeProteinSpecificityPrior}. Every rescue is recorded
 * in {@link IntegrationState#rescueLog} (category, claim, likelihood, score, target).</p>
 *
 * <p>Without claim counts (unit tests, legacy callers) the prior falls back to a
 * flat {@code +alphaEss} on every claim under an uncovered category. If the state
 * has no {@code essentialFunctions}, the prior no-ops.</p>
 */
class EssentialityPrior implements Prior {

    private static final Logger log = LoggerFactory.getLogger(EssentialityPrior)
    private static final Set<String> ASPECT_ROOTS = ['GO:0008150', 'GO:0003674', 'GO:0005575'] as Set

    /** Per-descendant boost magnitude (log-odds). Learned by the benchmark. */
    double alphaEss = 1.5d

    /** Maximum total boost per claim, as a multiple of alphaEss. */
    double maxBoostMultiplier = 2.0d

    /** 'targeted' (scored candidate, pinned to rescueTarget) or 'legacy' (v2.1.0 tier ranking and boost). */
    String rescueMode = 'targeted'

    /** Posterior probability a rescued claim is pinned to (targeted mode). */
    double rescueTarget = 0.60d

    /** Lowest Stage 0 likelihood probability a candidate may have. */
    double minCandidate = 0.16d

    /** Candidates rescued per uncovered category (targeted mode). */
    int maxRescuePerCategory = 1

    /** Score weights (targeted mode), on the log-odds scale. */
    double wSupport = 0.5d
    double wFree = 0.5d
    double wCompeting = 1.0d
    double wContext = 1.0d

    /** Cached per-iteration: union of GO descendants of every uncovered essential function (fallback when claimCounts is empty). */
    private Set<String> boostableDescendants = Collections.emptySet()

    /** Cached per-iteration: functionKey -> log-odds adjustment when claimCounts is populated. */
    private Map<String, Double> keyBoosts = Collections.emptyMap()

    /** Cached per-iteration: number of boost sources, for per-claim capping. */
    private int uncoveredCount = 0

    /** Rescues already written to the state's log (one entry per claim per run). */
    private final Set<String> logged = new HashSet<>()

    @Override
    String name() { 'essentiality' }

    @Override
    Set<String> inputs() { ['essentialFunctions', 'goOntology', 'goReasoner', 'posteriors'] as Set }

    private Set<String> closure(IntegrationState state, String term) {
        Set<String> s = new LinkedHashSet<>()
        s << term
        if (state.goReasoner != null) {
            try { s.addAll(state.goReasoner.getSubClasses(term, false)) }
            catch (Exception ex) { log.warn("EssentialityPrior: getSubClasses failed for ${term}: ${ex.message}") }
        } else if (state.goOntology != null) {
            try { s.addAll(state.goOntology.getDescendants(term)) }
            catch (Exception ex) { log.warn("EssentialityPrior: getDescendants failed for ${term}: ${ex.message}") }
        }
        s
    }

    private static Set<String> ancestorsOf(IntegrationState state, String term) {
        if (state.goOntology == null) return Collections.emptySet()
        try { return state.goOntology.getAncestors(term) ?: Collections.emptySet() }
        catch (Exception ignored) { return Collections.emptySet() }
    }

    @Override
    void beginIteration(IntegrationState state) {
        boostableDescendants = Collections.emptySet()
        keyBoosts = Collections.emptyMap()
        uncoveredCount = 0

        if (state.essentialFunctions == null) return

        Set<String> annotated = (state.goOntology != null || state.goReasoner != null)
            ? state.baselineAnnotatedGoTermsPropagated()
            : state.baselineAnnotatedGoTerms()
        Set<String> essentials = state.essentialFunctions.getGoTerms()
        Set<String> uncovered = new LinkedHashSet<>()
        for (String e : essentials) {
            if (!annotated.contains(e)) uncovered << e
        }
        uncoveredCount = uncovered.size()

        Map<String, Set<String>> descByUncovered = new LinkedHashMap<>()
        Set<String> descendants = new LinkedHashSet<>()
        for (String u : uncovered) {
            Set<String> uSet = closure(state, u)
            descByUncovered[u] = uSet
            descendants.addAll(uSet)
        }
        if (state.goReasoner == null && state.goOntology == null && !uncovered.isEmpty()) {
            log.info("EssentialityPrior: no GoOntology/GoReasoner, boosting ${uncovered.size()} exact essential terms (no descendant expansion)")
        }
        boostableDescendants = descendants

        if (state.claimCounts != null && !state.claimCounts.isEmpty()) {
            Set<String> allEssDomain = new LinkedHashSet<>(essentials)
            for (String e : essentials) allEssDomain.addAll(closure(state, e))
            Map<String, Double> boosts = state.computeProteinSpecificityPrior(null, allEssDomain, 0.30d * alphaEss, 0.30d * alphaEss)
            if (rescueMode == 'legacy') {
                legacyRescue(state, descByUncovered, boosts)
            } else {
                targetedRescue(state, descByUncovered, allEssDomain, boosts)
            }
            keyBoosts = boosts
        }
        log.debug("EssentialityPrior: ${uncovered.size()} uncovered essentials, ${descendants.size()} boostable descendants")
    }

    /** v2.1.0 rescue: tier ranking, boost max(1.30 alphaEss, 0.14 - L0) on the top candidate. */
    private void legacyRescue(IntegrationState state, Map<String, Set<String>> descByUncovered, Map<String, Double> boosts) {
        for (Map.Entry<String, Set<String>> ue : descByUncovered.entrySet()) {
            List<Map.Entry<String, Double>> cands = candidates(state, ue.value)
            state.sortCandidatesForRescue(cands)
            int limit = Math.min(1, cands.size())
            for (int i = 0; i < limit; i++) {
                String fk = cands[i].key
                double bVal = Math.max(Math.min(maxBoostMultiplier * alphaEss, 1.30d * alphaEss), 0.14d - state.getLikelihoodLogOdds(fk))
                boosts[fk] = Math.max(boosts.getOrDefault(fk, 0.0d), bVal)
                record(state, ue.key, fk, cands[i].value, Double.NaN, Double.NaN)
            }
        }
    }

    /** Borderline candidate claims under a category closure: (functionKey, likelihood probability). */
    private List<Map.Entry<String, Double>> candidates(IntegrationState state, Set<String> uSet) {
        List<Map.Entry<String, Double>> cands = []
        for (String fkey : state.posteriorLogOdds.keySet()) {
            String[] parts = IntegrationState.splitFunctionKey(fkey)
            if (parts == null || parts[1] != 'GO') continue
            if (!uSet.contains(parts[2])) continue
            double likeProb = state.getLikelihoodProbability(fkey)
            if (likeProb >= minCandidate && likeProb <= IntegrationState.MAP_THRESHOLD) {
                cands << new AbstractMap.SimpleImmutableEntry<>(fkey, likeProb)
            }
        }
        cands
    }

    private void targetedRescue(IntegrationState state, Map<String, Set<String>> descByUncovered,
                                Set<String> allEssDomain, Map<String, Double> boosts) {
        // Confident Stage 0 GO terms per protein, and operon partners per protein, indexed once per iteration.
        Map<String, Set<String>> strongByProt = new HashMap<>()
        Map<String, Double> src = (state.likelihoodLogOdds != null && !state.likelihoodLogOdds.isEmpty()) ? state.likelihoodLogOdds : state.posteriorLogOdds
        for (Map.Entry<String, Double> e : src.entrySet()) {
            if (e.value < 0.0d) continue   // log-odds >= 0 <=> probability >= 0.50
            String[] parts = IntegrationState.splitFunctionKey(e.key)
            if (parts == null || parts[1] != 'GO') continue
            strongByProt.computeIfAbsent(parts[0], { new HashSet<String>() }) << parts[2]
        }
        Map<String, List<String>> mates = new HashMap<>()
        if (state.operons != null) {
            for (List<String> op : state.operons) {
                for (String p : op) mates[p] = op.findAll { it != p }
            }
        }
        double targetLogOdds = Math.log(rescueTarget / (1.0d - rescueTarget))

        for (Map.Entry<String, Set<String>> ue : descByUncovered.entrySet()) {
            String category = ue.key
            Set<String> uSet = ue.value
            Set<String> ctxTerms = new HashSet<>(allEssDomain)
            ctxTerms.addAll(ancestorsOf(state, category))
            ctxTerms.removeAll(ASPECT_ROOTS)
            String catAspect = state.goOntology != null ? state.goOntology.getAspect(category) : null

            List<Map.Entry<String, Double>> cands = candidates(state, uSet)
            if (cands.isEmpty()) continue
            Map<String, Double> score = new HashMap<>()
            for (Map.Entry<String, Double> c : cands) {
                String fk = c.key
                String[] parts = IntegrationState.splitFunctionKey(fk)
                String pid = parts[0]; String term = parts[2]
                double s = state.getLikelihoodLogOdds(fk)
                s += wSupport * Math.log(Math.max(1, state.getClaimCount(fk)))
                boolean free = state.getMaxLikelihoodProbForProtein(pid) < IntegrationState.MAP_THRESHOLD
                if (free) s += wFree
                if (!free && catAspect != null) {
                    Set<String> anc = ancestorsOf(state, term)
                    boolean competing = (strongByProt.get(pid) ?: Collections.<String>emptySet()).any { String g ->
                        !uSet.contains(g) && !anc.contains(g) && !ASPECT_ROOTS.contains(g) && state.goOntology.getAspect(g) == catAspect
                    }
                    if (competing) s -= wCompeting
                }
                List<String> ms = mates.get(pid)
                if (ms != null && !ms.isEmpty()) {
                    int hit = 0
                    for (String m : ms) {
                        Set<String> sg = strongByProt.get(m)
                        if (sg != null && sg.any { ctxTerms.contains(it) }) hit++
                    }
                    s += wContext * (hit / (double) ms.size())
                }
                score[fk] = s
            }
            cands.sort { a, b ->
                int c = Double.compare(score[b.key], score[a.key])
                c != 0 ? c : a.key <=> b.key
            }
            int limit = Math.min(maxRescuePerCategory, cands.size())
            for (int i = 0; i < limit; i++) {
                String fk = cands[i].key
                double l0 = state.getLikelihoodLogOdds(fk)
                double bVal = Math.max(Math.min(maxBoostMultiplier * alphaEss, 1.30d * alphaEss), targetLogOdds - l0)
                boosts[fk] = Math.max(boosts.getOrDefault(fk, 0.0d), bVal)
                Double prev = state.getPinnedFloor(ClaimKey.parse(fk))
                if (prev == null || prev < targetLogOdds) state.setPinnedFloor(ClaimKey.parse(fk), targetLogOdds)
                record(state, category, fk, cands[i].value, score[fk], rescueTarget)
            }
        }
    }

    private void record(IntegrationState state, String category, String fk, double likeProb, double score, double target) {
        String id = "${category}|${fk}".toString()
        if (!logged.add(id)) return
        state.rescueLog << [prior: name(), mode: rescueMode, category: category, functionKey: fk,
                            likelihood: likeProb, score: score, target: target]
    }

    @Override
    double logOddsBoost(String proteinId, String functionKey, IntegrationState state) {
        if (state.claimCounts != null && !state.claimCounts.isEmpty()) {
            return keyBoosts.getOrDefault(functionKey, 0.0d)
        }
        if (boostableDescendants.isEmpty()) return 0.0d
        String[] parts = IntegrationState.splitFunctionKey(functionKey)
        if (parts == null || parts[1] != 'GO') return 0.0d
        if (!boostableDescendants.contains(parts[2])) return 0.0d

        return Math.min(maxBoostMultiplier * alphaEss, alphaEss)
    }
}
