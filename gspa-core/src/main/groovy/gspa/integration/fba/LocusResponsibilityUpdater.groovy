package gspa.integration.fba

import groovy.transform.CompileStatic
import gspa.integration.ClaimProvenance
import gspa.integration.IntegratedAnnotationSet
import gspa.integration.IntegrationState
import gspa.model.AnnotationType

/**
 * Stage C: Backward Locus Responsibility Attribution r_j(h).
 *
 * Distributes the flux-certified repair signal pi_j^FBA * Delta_repair(j)
 * across candidate genes h in H(r) U {h_null} using:
 *   w_{j,r}(h) = (s_{h,r}/(1-s_{h,r}))^beta1 * BF_RLGC(h,r)^beta2 * exp(beta3*C_syn + beta4*S_PLM)
 *   r_{j,r}(h) = w_{j,r}(h) / (tau_null + sum_{h'} w_{j,r}(h'))
 * and updates log-odds posteriors in IntegratedAnnotationSet.
 */
@CompileStatic
class LocusResponsibilityUpdater {

    double beta1 = 1.0d
    double beta2 = 1.0d
    double beta3 = 1.5d
    double beta4 = 1.5d
    double tauNull = 0.05d
    double gammaFba = 2.2d

    static double logit(double p) {
        double pc = Math.min(Math.max(p, 1e-6d), 1.0d - 1e-6d)
        return Math.log(pc / (1.0d - pc))
    }

    static double sigmoid(double x) {
        if (x >= 0.0d) {
            double z = Math.exp(-x)
            return 1.0d / (1.0d + z)
        } else {
            double z = Math.exp(x)
            return z / (1.0d + z)
        }
    }

    /**
     * Apply updated gene-EC posterior probabilities from StoichiometricLocusRefiner
     * back into IntegratedAnnotationSet provenance and annotations.
     */
    int applyUpdatedScores(
        IntegratedAnnotationSet integrated,
        Map<String, Map<String, Double>> updatedGeneEcScores,
        Map<String, String> ecToGoMap = Collections.emptyMap()
    ) {
        if (integrated == null || updatedGeneEcScores == null || updatedGeneEcScores.isEmpty()) {
            return 0
        }
        int updatedCount = 0
        for (Map.Entry<String, Map<String, Double>> gEntry : updatedGeneEcScores.entrySet()) {
            String proteinId = gEntry.key
            Map<String, Double> ecMap = gEntry.value
            if (ecMap == null) continue
            for (Map.Entry<String, Double> eEntry : ecMap.entrySet()) {
                String rawEc = eEntry.key.replace('EC:', '').trim()
                double newProb = eEntry.value != null ? eEntry.value.doubleValue() : 0.0d
                for (String ecKey : [rawEc, "EC:${rawEc}".toString()]) {
                    ClaimProvenance prov = integrated.provenanceFor(proteinId, AnnotationType.EC, ecKey)
                    if (prov != null && newProb > prov.finalProbability + 1e-5d) {
                        double oldLogOdds = prov.finalLogOdds
                        double newLogOdds = logit(newProb)
                        prov.priorContributions['fba_locus_refiner'] = newLogOdds - oldLogOdds
                        prov.finalLogOdds = newLogOdds
                        prov.finalProbability = newProb
                        updatedCount++
                    }
                }
            }
        }
        return updatedCount
    }
}
