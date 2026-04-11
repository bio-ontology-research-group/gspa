package gspa.integration

/**
 * Combines a set of {@link EvidenceClaim}s for the same (protein, function)
 * into a likelihood log-odds via Noisy-OR with correlation-group collapse.
 *
 * <p>Algorithm (plan §A.3):</p>
 * <ol>
 *   <li>Group claims by {@code evidenceType.correlationGroup()}.</li>
 *   <li>Within each group, take the single strongest claim —
 *       {@code argmax(reliability[type] * calibratedProb)}. This addresses
 *       non-independence between homology-family tools like DIAMOND and
 *       eggNOG that draw from overlapping databases.</li>
 *   <li>Compute Noisy-OR across groups:
 *       {@code p = 1 − ∏ (1 − w_t * p_g)}.</li>
 *   <li>Convert to log-odds, clip to [lMin, lMax].</li>
 * </ol>
 *
 * <p>Per-{@link EvidenceType} reliability weights are either learned by the
 * Phase 7.4 benchmark or set to the defaults in {@link #defaultReliability()}.</p>
 */
class EvidenceCombiner {

    /** Lower clip on log-odds to keep the combiner numerically stable. */
    double lMin = -12.0

    /** Upper clip on log-odds. */
    double lMax = 12.0

    /**
     * Optional bonus (log-odds) added when two claims within the same
     * correlation group come from genuinely distinct database sources
     * (metadata-tracked). Phase 7.4 learns this. Default 0 — disabled.
     */
    double secondOpinionBonus = 0.0

    /** Per-evidence-type reliability weights in [0, 1]. */
    Map<EvidenceType, Double> reliability

    EvidenceCombiner(Map<EvidenceType, Double> reliability = null) {
        this.reliability = reliability ?: defaultReliability()
    }

    /**
     * Default reliability weights. Used when the benchmark has not yet
     * written learned values. Chosen so sequence homology and orthology
     * receive modest weight (they're the most commonly invoked sources),
     * structure-based evidence receives slightly higher weight (it's more
     * independent), and context signals are conservative.
     */
    static Map<EvidenceType, Double> defaultReliability() {
        Map<EvidenceType, Double> r = new EnumMap<>(EvidenceType)
        EvidenceType.values().each { r[it] = 0.6d }
        r[EvidenceType.SEQUENCE_SIMILARITY]      = 0.70d
        r[EvidenceType.SEQUENCE_DOMAIN]          = 0.75d
        r[EvidenceType.STRUCTURE_SIMILARITY]     = 0.80d
        r[EvidenceType.STRUCTURE_DEEPLEARNING]   = 0.75d
        r[EvidenceType.ORTHOLOGY]                = 0.70d
        r[EvidenceType.SEQUENCE_DEEPLEARNING]    = 0.65d
        r[EvidenceType.PROTEIN_LM_EMBEDDING]     = 0.65d
        r[EvidenceType.GENOMIC_CONTEXT]          = 0.45d
        r[EvidenceType.METABOLIC_CONTEXT]        = 0.55d
        r[EvidenceType.GENOMIC_LANGUAGE_MODEL]   = 0.55d
        r[EvidenceType.LOCALIZATION]             = 0.50d
        r[EvidenceType.DOMAIN_SPECIFIC_AMR]      = 0.85d
        r[EvidenceType.DOMAIN_SPECIFIC_CAZY]     = 0.80d
        r[EvidenceType.DOMAIN_SPECIFIC_BGC]      = 0.75d
        r[EvidenceType.DOMAIN_SPECIFIC_VF]       = 0.75d
        r[EvidenceType.SEQUENCE_MOTIF]           = 0.50d
        r
    }

    /**
     * Combine a list of claims about the same (protein, function) into a
     * posterior log-odds. Claims must share the same
     * {@link EvidenceClaim#functionKey()} — no check is made, the caller is
     * responsible.
     *
     * @return log-odds in [{@link #lMin}, {@link #lMax}]
     */
    double combineLikelihood(List<EvidenceClaim> claimsForSameFunction) {
        if (claimsForSameFunction == null || claimsForSameFunction.isEmpty()) {
            return lMin
        }

        // Step 1 + 2: group by correlation group, pick the strongest in each group.
        Map<String, List<EvidenceClaim>> byGroup = [:]
        for (EvidenceClaim c : claimsForSameFunction) {
            String g = c.evidenceType.correlationGroup()
            byGroup.computeIfAbsent(g, { k -> new ArrayList<EvidenceClaim>() }) << c
        }

        double oneMinusProd = 1.0d
        for (Map.Entry<String, List<EvidenceClaim>> entry : byGroup.entrySet()) {
            EvidenceClaim best = null
            double bestScore = -1.0d
            EvidenceClaim secondBest = null
            double secondBestScore = -1.0d
            for (EvidenceClaim c : entry.value) {
                double w = reliability.getOrDefault(c.evidenceType, 0.5d)
                double s = w * c.calibratedProb
                if (s > bestScore) {
                    secondBest = best
                    secondBestScore = bestScore
                    best = c
                    bestScore = s
                } else if (s > secondBestScore) {
                    secondBest = c
                    secondBestScore = s
                }
            }

            double p = bestScore
            // Second-opinion bonus: only if the two claims are from distinct databases.
            if (secondOpinionBonus > 0.0 && secondBest != null
                    && distinctDatabases(best, secondBest)) {
                p += secondOpinionBonus * secondBestScore
            }
            p = Math.min(0.999d, Math.max(0.001d, p))
            oneMinusProd *= (1.0d - p)
        }

        // Step 3: Noisy-OR combination across groups.
        double pPost = 1.0d - oneMinusProd
        pPost = Math.min(0.999999d, Math.max(1.0e-6d, pPost))

        // Step 4: convert to log-odds + clip.
        double logOdds = Math.log(pPost / (1.0d - pPost))
        Math.min(lMax, Math.max(lMin, logOdds))
    }

    /**
     * Heuristic: two within-group claims count as "distinct databases"
     * if they name different source predictors AND their metadata records
     * different target databases (or metadata is absent). The benchmark
     * can refine this later.
     */
    private static boolean distinctDatabases(EvidenceClaim a, EvidenceClaim b) {
        if (a.source && b.source && a.source != b.source) {
            Object dbA = a.metadata?.get('database')
            Object dbB = b.metadata?.get('database')
            if (dbA == null || dbB == null) return true
            return dbA != dbB
        }
        return false
    }
}
