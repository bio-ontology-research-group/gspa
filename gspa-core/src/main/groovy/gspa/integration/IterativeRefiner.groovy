package gspa.integration

import gspa.model.Annotation
import gspa.model.AnnotationType
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Runs the evidence-integration fixed-point loop.
 *
 * <p>Phase 7.1: single-pass likelihood combination only (no priors).
 * Phase 7.3 adds the {@code PriorEngine}, damping, and the convergence
 * check that turns this into the multi-iteration refiner described in
 * plan §A.5.</p>
 *
 * <p>Result: {@link IntegratedAnnotationSet} with per-function provenance
 * and final posterior annotations whose score is the posterior
 * probability.</p>
 */
class IterativeRefiner {

    private static final Logger log = LoggerFactory.getLogger(IterativeRefiner)

    EvidenceCombiner combiner

    /** Not used in Phase 7.1; populated in Phase 7.3. */
    PriorEngine priorEngine = null

    /** Iteration cap (Phase 7.3). */
    int maxIter = 6

    /** Mean absolute delta in posterior probability to stop iterating. */
    double epsilon = 0.005

    /** Under-relaxation factor (Phase 7.3). */
    double damping = 0.5

    IterativeRefiner(EvidenceCombiner combiner) {
        this.combiner = combiner
    }

    /**
     * Refine a set of claims into posterior annotations.
     *
     * Phase 7.1 implementation: single pass. Combine each (protein,
     * function) group once, convert to probability, emit one Annotation
     * per function key.
     */
    IntegratedAnnotationSet refine(List<EvidenceClaim> claims, IntegrationState state) {
        Map<String, List<EvidenceClaim>> byKey = ClaimExtractor.groupByFunctionKey(claims)
        log.info("Refining ${claims.size()} claims across ${byKey.size()} (protein, function) groups")

        IntegratedAnnotationSet out = new IntegratedAnnotationSet()
        Map<String, Double> posteriorLogOdds = new LinkedHashMap<>()

        for (Map.Entry<String, List<EvidenceClaim>> entry : byKey.entrySet()) {
            String key = entry.key
            List<EvidenceClaim> group = entry.value

            double lLik = combiner.combineLikelihood(group)
            double lPost = lLik   // Phase 7.1: no prior contribution.
            double pPost = sigmoid(lPost)

            posteriorLogOdds[key] = lPost

            ClaimProvenance prov = new ClaimProvenance(
                functionKey: key,
                proteinId: group.first().proteinId,
                supportingClaims: new ArrayList<EvidenceClaim>(group),
                likelihoodLogOdds: lLik,
                finalLogOdds: lPost,
                finalProbability: pPost,
                convergenceIter: 0,
            )

            Annotation ann = buildAnnotation(group.first(), pPost)
            out.put(prov, ann)
        }

        state.updatePosteriors(posteriorLogOdds)
        out
    }

    /**
     * Build the final Annotation for a posterior. Uses the first claim's
     * metadata (type, functionId, aspect) — all claims in the group share
     * these by construction. The score is the final posterior probability.
     */
    private static Annotation buildAnnotation(EvidenceClaim claim, double posteriorProb) {
        new Annotation(
            type: claim.functionType,
            value: claim.functionId,
            score: posteriorProb,
            source: 'integrated',
            evidence: 'IEA',
            goAspect: claim.goAspect,
            evidenceType: claim.evidenceType,
        )
    }

    private static double sigmoid(double x) {
        if (x >= 500.0) return 1.0
        if (x <= -500.0) return 0.0
        1.0 / (1.0 + Math.exp(-x))
    }
}
