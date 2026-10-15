package gspa.integration.fba

import groovy.transform.CompileStatic

/**
 * Stage A2: Sparse Top-M Isozyme OR & Geometric-Mean Subunit AND Aggregation.
 *
 * Replaces METEOR's sum-of-log-odds over all paralogs (which inflates broad EC classes)
 * and product-of-probabilities over complex subunits (which crushes multi-subunit enzymes)
 * with:
 *   1. Top-M (default M=3) size-corrected Noisy-OR over isozymes for each catalytic role r
 *   2. Geometric-mean subunit bottleneck for heteromeric AND complexes k = {r_1, ..., r_L}
 *   3. Maximum over alternative mechanisms k in M(j): P_j = max_k P_k
 */
@CompileStatic
class SparseRoleAggregator {

    int topM = 3
    double tauSub = 0.25d
    double epsMin = 1e-4d

    /**
     * Aggregate gene-level role scores s_{g,r} into role presence probabilities p_r in (0, 1).
     *
     * @param geneRoleScores Map of geneId -> (roleId -> probability)
     * @param roleBackground Optional map of roleId -> background false-positive baseline b_r
     * @return Map of roleId -> aggregated role probability p_r
     */
    Map<String, Double> aggregateRoleProbabilities(
        Map<String, Map<String, Double>> geneRoleScores,
        Map<String, Double> roleBackground = Collections.emptyMap()
    ) {
        Map<String, List<Double>> roleTopScores = new LinkedHashMap<String, List<Double>>()
        for (Map.Entry<String, Map<String, Double>> geneEntry : geneRoleScores.entrySet()) {
            Map<String, Double> rMap = geneEntry.value
            if (rMap == null) continue
            for (Map.Entry<String, Double> rEntry : rMap.entrySet()) {
                String cleanRole = rEntry.key.replace('EC:', '').trim()
                double s = rEntry.value != null ? rEntry.value.doubleValue() : 0.0d
                if (s <= 0.0d) continue
                List<Double> list = roleTopScores.get(cleanRole)
                if (list == null) {
                    list = new ArrayList<Double>()
                    roleTopScores.put(cleanRole, list)
                }
                list.add(s)
            }
        }

        Map<String, Double> roleProbs = new LinkedHashMap<String, Double>()
        for (Map.Entry<String, List<Double>> entry : roleTopScores.entrySet()) {
            String role = entry.key
            List<Double> scores = entry.value
            scores.sort(Comparator.reverseOrder())
            int limit = Math.min(topM, scores.size())
            double bR = roleBackground.containsKey(role) ? roleBackground.get(role).doubleValue() : 0.0d
            double denom = Math.max(1.0d - bR, 1e-6d)
            double prodComp = 1.0d
            for (int i = 0; i < limit; i++) {
                double s = scores.get(i).doubleValue()
                double sTilde = Math.max(0.0d, (s - bR) / denom)
                sTilde = Math.min(Math.max(sTilde, 0.0d), 1.0d - epsMin)
                prodComp *= (1.0d - sTilde)
            }
            roleProbs.put(role, 1.0d - prodComp)
        }
        return roleProbs
    }

    /**
     * Evaluate mechanism probability P_k for a compiled DNF clause k = {r_1, ..., r_L}.
     */
    double evaluateMechanism(GprMechanismCompiler.MechanismClause clause, Map<String, Double> roleProbs) {
        if (clause == null || clause.requiredRoles.isEmpty()) return 0.0d
        if (!clause.isHeteromericComplex || clause.requiredRoles.size() == 1) {
            String r = clause.requiredRoles.get(0)
            return roleProbs.containsKey(r) ? roleProbs.get(r).doubleValue() : 0.0d
        }
        double minVal = 1.0d
        double logSum = 0.0d
        int L = clause.requiredRoles.size()
        for (String r : clause.requiredRoles) {
            double p = roleProbs.containsKey(r) ? roleProbs.get(r).doubleValue() : 0.0d
            if (p < minVal) minVal = p
            logSum += Math.log(Math.max(p, 1e-6d))
        }
        if (minVal >= tauSub) {
            return Math.exp(logSum / (double) L)
        }
        return minVal
    }
}
