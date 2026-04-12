package gspa.integration.suggester

import gspa.integration.IntegratedAnnotationSet
import gspa.integration.IntegrationState
import gspa.integration.MetabolicGap
import gspa.model.AnnotationType
import gspa.ontology.PathwayGraph
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Phase 8 "Dark Matter / Contextual Gap" Suggester.
 *
 * <p>Takes a metabolic gap {@code (P, R)} and the genomic context
 * (operons + integrated posteriors) and assigns the missing function to
 * the most likely protein in the most likely operon. Implements the
 * four-layer Bayesian algorithm from plan §A.1:</p>
 *
 * <ol>
 *   <li><b>Layer 1</b> — Bayes factor {@code BF(O, P)} that operon O
 *       participates in pathway P, computed as a soft-weighted product
 *       over operon members' annotations (no "present / absent"
 *       threshold).</li>
 *   <li><b>Layer 2</b> — per protein
 *       {@code L_R(p) = L_lik + L_op + L_lm}, where {@code L_op =
 *       log(BF/(1+BF))} is shared across the operon and {@code L_lik}
 *       is the integrator's current posterior for {@code (p, f_R)} or a
 *       low "absent" default.</li>
 *   <li><b>Layer 3</b> — softmax over operon members:
 *       {@code q(p) = π_R(p) / Σ π_R(p')}.</li>
 *   <li><b>Layer 4</b> — singleton if the top q &gt; 0.5, otherwise a
 *       disjunctive suggestion over the smallest top-k whose cumulative
 *       q exceeds {@code coverageThreshold}.</li>
 * </ol>
 *
 * <p>Phase 8 is boost/suggest-only — suggestions never re-enter the
 * refinement loop. They are written to
 * {@link IntegratedAnnotationSet#suggestions} as a separate channel.</p>
 */
class DarkMatterSuggester {

    private static final Logger log = LoggerFactory.getLogger(DarkMatterSuggester)

    // ------------------------------------------------------------------
    // Hyperparameters (all learnable by the benchmark)
    // ------------------------------------------------------------------

    /** Minimum Bayes factor for candidate operons. Jeffreys "strong" = 10. */
    double bfMin = 10.0d

    /**
     * Likelihood-elevation factor {@code gamma_in_P}: how much more
     * likely an annotated function is under "operon is in pathway P"
     * than under the base rate. Default 50 (fairly strong).
     */
    double gammaInP = 50.0d

    /** Cumulative q coverage for disjunctive output. */
    double coverageThreshold = 0.9d

    /** Reference Bayes factor for the tanh squash in scoring. */
    double bfRef = 100.0d

    /** Default pathway confidence (curated pathways). Lower for predicted. */
    double defaultPathwayConfidence = 1.0d

    /** Cap on suggestionScore. */
    double maxSuggestionScore = 0.85d

    /**
     * Default likelihood log-odds used when the integrator has no
     * posterior for {@code (p, f_R)}. Negative because absence of a
     * specific GO term is the prior expectation. Corresponds to
     * probability ~0.0025. Not learned.
     */
    double absentLikelihoodLogOdds = -6.0d

    /**
     * Floor on the "free slot" probability used in the commitment
     * penalty. Prevents {@code log(0)} when a protein is annotated for
     * another pathway function with probability essentially 1.
     */
    double minFreeProbability = 1.0e-4d

    // ------------------------------------------------------------------
    // Main entry point
    // ------------------------------------------------------------------

    /**
     * Run the suggester against an already-converged integration state.
     *
     * @param state converged state (posteriorLogOdds populated by refiner)
     * @param integrated result of refinement; suggestions are attached to it
     * @return the same {@link IntegratedAnnotationSet} with {@code suggestions} populated
     */
    IntegratedAnnotationSet suggest(IntegrationState state, IntegratedAnnotationSet integrated) {
        integrated.suggestions = []

        if (state.metabolicGaps == null || state.metabolicGaps.isEmpty()) {
            log.info("DarkMatterSuggester: no metabolic gaps; skipping")
            return integrated
        }
        if (state.operons == null || state.operons.isEmpty()) {
            log.info("DarkMatterSuggester: no operons wired; skipping")
            return integrated
        }
        if (state.pathwayDatabase == null) {
            log.info("DarkMatterSuggester: no pathway database wired; skipping")
            return integrated
        }

        // Precompute required-GO-term sets per pathway.
        Map<String, Set<String>> pathwayTerms = new LinkedHashMap<>()
        for (Map.Entry<String, PathwayGraph> e : state.pathwayDatabase.pathways.entrySet()) {
            Set<String> req = e.value.getRequiredGoTerms() ?: new LinkedHashSet<String>()
            if (!req.isEmpty()) pathwayTerms[e.key] = req
        }

        // Assign stable operon IDs for provenance.
        Map<Integer, String> operonIds = new LinkedHashMap<>()
        state.operons.eachWithIndex { op, i -> operonIds[i] = "operon_${i}".toString() }

        // Reverse index: GO term → set of pathway IDs whose required set contains it.
        Map<String, Set<String>> goToPathways = new LinkedHashMap<>()
        pathwayTerms.each { pwId, terms ->
            terms.each { goId ->
                goToPathways.computeIfAbsent(goId, { new LinkedHashSet<String>() }) << pwId
            }
        }

        int emitted = 0
        for (MetabolicGap gap : state.metabolicGaps) {
            if (!gap.goTerm) continue                       // no target function → nothing to assign

            // Find pathways whose required-term set contains the gap's target GO.
            // The gap's own pathwayId is MetaCyc (from gapseq), but the pathway DB
            // may use KEGG IDs. Bridge via the GO term itself.
            Set<String> candidatePathways = goToPathways[gap.goTerm]
            if (candidatePathways == null || candidatePathways.isEmpty()) {
                // Also try via the gap's pathwayId directly (works if DB and gap use same namespace).
                Set<String> direct = pathwayTerms[gap.pathwayId]
                if (direct == null || direct.isEmpty()) continue
                candidatePathways = [gap.pathwayId] as Set
            }

            // Use the union of all candidate pathways' required terms.
            Set<String> pfTerms = new LinkedHashSet<>()
            candidatePathways.each { pwId -> pfTerms.addAll(pathwayTerms[pwId] ?: []) }

            for (int i = 0; i < state.operons.size(); i++) {
                List<String> operon = state.operons[i]
                if (operon == null || operon.size() < 1) continue

                // Layer 1: Bayes factor for this operon / pathway.
                double bf = computeBayesFactor(operon, pfTerms, state)
                if (bf < bfMin) continue

                // Layer 2: per-protein L_R decomposition.
                Map<String, PerProteinDecomposition> decomposed = buildDecompositions(
                    operon, gap.goTerm, bf, pfTerms, state,
                )
                // Filter out proteins whose π_R is effectively 0 (avoid divide-by-0).
                double z = decomposed.values().sum { it.piR } as double
                if (z <= 1e-12) continue

                // Layer 3: softmax within the operon.
                decomposed.values().each { d -> d.q = d.piR / z }

                // Layer 4: decide singleton vs disjunctive.
                List<PerProteinDecomposition> sorted = decomposed.values()
                    .toSorted { a, b -> b.q <=> a.q }

                double qMinEmit = 1.0d / operon.size()
                if (sorted[0].q < qMinEmit) continue        // operon is uninformative

                Suggestion s = decideAndBuild(gap, operon, operonIds[i], bf, sorted, decomposed, pfTerms)
                if (s != null) {
                    integrated.suggestions << s
                    emitted++
                }
            }
        }

        log.info("DarkMatterSuggester: emitted ${emitted} suggestions over ${state.metabolicGaps.size()} gaps, ${state.operons.size()} operons")
        integrated
    }

    // ------------------------------------------------------------------
    // Layer 1: Bayes factor
    // ------------------------------------------------------------------

    /**
     * Compute BF(O, P). Each annotated (p, f) pair contributes a soft
     * weighted factor of {@code gammaInP^sigmoid(posterior)} if
     * {@code f ∈ functions(P)}, and no factor otherwise.
     *
     * <p>In log space:</p>
     * <pre>
     *   log BF = log(gammaInP) * Σ_{p ∈ O, f ∈ functions(P)} sigmoid(posterior(p, f))
     * </pre>
     *
     * <p>This yields BF = 1 when no annotations lie in functions(P),
     * and grows exponentially with the total soft weight of pathway-
     * matching annotations.</p>
     */
    private double computeBayesFactor(List<String> operon, Set<String> pfTerms, IntegrationState state) {
        double totalWeight = 0.0d
        for (String pid : operon) {
            for (String key : state.functionKeysForProtein(pid)) {
                String[] parts = IntegrationState.splitFunctionKey(key)
                if (parts == null || parts[1] != 'GO') continue
                if (!pfTerms.contains(parts[2])) continue
                totalWeight += sigmoid(state.posteriorLogOdds.getOrDefault(key, 0.0d))
            }
        }
        Math.exp(totalWeight * Math.log(gammaInP))
    }

    // ------------------------------------------------------------------
    // Layer 2: per-protein posterior composition
    // ------------------------------------------------------------------

    private Map<String, PerProteinDecomposition> buildDecompositions(
            List<String> operon,
            String targetFunctionId,
            double bf,
            Set<String> pathwayFunctions,
            IntegrationState state) {
        double lOp = Math.log(bf / (1.0d + bf))
        Map<String, PerProteinDecomposition> out = new LinkedHashMap<>()
        for (String pid : operon) {
            String key = "${pid}|GO|${targetFunctionId}".toString()
            double lLik = state.posteriorLogOdds.containsKey(key)
                ? state.posteriorLogOdds[key]
                : absentLikelihoodLogOdds
            double lLm = 0.0d   // Phase 9 will populate this from a genomic-LM source
            // Commitment penalty: a protein already strongly annotated for
            // OTHER pathway functions is less likely to also be the one
            // performing the missing reaction. We measure this as the
            // log-probability that p's "pathway slot" is still free,
            //   log(1 - max_{f' ∈ functions(P) \ {f_R}} sigmoid(posterior(p, f'))).
            // For a fully-dark protein this is log(1) = 0; for a protein
            // committed to a different pathway function it is very negative.
            double maxOther = maxOtherPathwayPosterior(pid, targetFunctionId, pathwayFunctions, state)
            double freeProb = Math.max(minFreeProbability, 1.0d - maxOther)
            double lCommit = Math.log(freeProb)
            double lR = lLik + lOp + lLm + lCommit
            double pR = sigmoid(lR)
            out[pid] = new PerProteinDecomposition(
                proteinId: pid,
                likelihoodLogOdds: lLik + lCommit,   // absorb commitment into the likelihood term
                operonLogOdds: lOp,
                lmLogOdds: lLm,
                totalLogOdds: lR,
                piR: pR,
            )
        }
        out
    }

    /**
     * Max posterior probability of a pathway function OTHER than the
     * target on this protein. Returns 0 if the protein has no
     * annotations for pathway functions. Used by the commitment penalty.
     */
    private double maxOtherPathwayPosterior(
            String pid, String targetFunctionId, Set<String> pathwayFunctions, IntegrationState state) {
        double best = 0.0d
        for (String key : state.functionKeysForProtein(pid)) {
            String[] parts = IntegrationState.splitFunctionKey(key)
            if (parts == null || parts[1] != 'GO') continue
            String go = parts[2]
            if (go == targetFunctionId) continue
            if (!pathwayFunctions.contains(go)) continue
            double p = sigmoid(state.posteriorLogOdds.getOrDefault(key, 0.0d))
            if (p > best) best = p
        }
        best
    }

    // ------------------------------------------------------------------
    // Layer 4: decide singleton vs disjunctive + score
    // ------------------------------------------------------------------

    private Suggestion decideAndBuild(
            MetabolicGap gap,
            List<String> operon,
            String operonId,
            double bf,
            List<PerProteinDecomposition> sorted,
            Map<String, PerProteinDecomposition> allDecomposed,
            Set<String> pathwayFunctions) {
        double qTop = sorted[0].q
        double bfSquash = Math.tanh(Math.log(bf) / Math.log(bfRef))
        double pathwayConfidence = defaultPathwayConfidence

        if (qTop > 0.5d) {
            double score = Math.min(maxSuggestionScore, qTop * bfSquash * pathwayConfidence)
            return new SingletonSuggestion(
                proteinId: sorted[0].proteinId,
                q: qTop,
                functionId: gap.goTerm,
                functionType: AnnotationType.GO,
                pathwayId: gap.pathwayId,
                reactionId: gap.reactionId,
                operonId: operonId,
                bayesFactor: bf,
                suggestionScore: score,
                proteinScores: allDecomposed,
                provenance: "gap=${gap.pathwayId}/${gap.reactionId}, operon=${operonId}, BF=${String.format(Locale.ROOT, '%.2f', bf)}, q=${String.format(Locale.ROOT, '%.2f', qTop)}".toString(),
            )
        }

        // Build the credible set up to coverageThreshold.
        List<String> pids = []
        List<Double> qs = []
        double cumulative = 0.0d
        for (PerProteinDecomposition d : sorted) {
            pids << d.proteinId
            qs << d.q
            cumulative += d.q
            if (cumulative >= coverageThreshold) break
        }
        if (pids.size() < 2) return null     // cumulative already covers but top-1 was below 0.5 — degenerate; skip

        // suggestionScore via normalized entropy: (1 - H/log(k)) * tanh(logBF/logBFref) * pathwayConf.
        double H = 0.0d
        for (double q : qs) {
            if (q > 0.0d) H -= q * Math.log(q)
        }
        double k = pids.size()
        double uniform = Math.log(k)
        double normalizedEntropy = uniform > 0 ? (H / uniform) : 0.0d
        double concentration = 1.0d - normalizedEntropy
        double score = Math.min(maxSuggestionScore, concentration * bfSquash * pathwayConfidence)

        return new DisjunctiveSuggestion(
            proteinIds: pids,
            qValues: qs,
            functionId: gap.goTerm,
            functionType: AnnotationType.GO,
            pathwayId: gap.pathwayId,
            reactionId: gap.reactionId,
            operonId: operonId,
            bayesFactor: bf,
            suggestionScore: score,
            proteinScores: allDecomposed,
            provenance: "gap=${gap.pathwayId}/${gap.reactionId}, operon=${operonId}, BF=${String.format(Locale.ROOT, '%.2f', bf)}, k=${(int) k}, coverage=${String.format(Locale.ROOT, '%.2f', cumulative)}".toString(),
        )
    }

    private static double sigmoid(double x) {
        if (x >= 500.0) return 1.0d
        if (x <= -500.0) return 0.0d
        1.0d / (1.0d + Math.exp(-x))
    }
}
