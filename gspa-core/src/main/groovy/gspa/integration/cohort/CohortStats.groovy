package gspa.integration.cohort

import groovy.transform.CompileStatic

/** Small numeric helpers shared by the cohort components. */
@CompileStatic
class CohortStats {

    /** Jeffreys-smoothed proportion; defined (0.5) when there are no observations. */
    static double jeffreys(double successes, double n) {
        (successes + 0.5d) / (n + 1.0d)
    }

    /**
     * Approximate 90% Wilson score interval on effective counts. It is an
     * approximation (effective, not integer, counts), not an exact credible
     * interval. Returns [0, 1] when n == 0.
     */
    static double[] wilson(double successes, double n) {
        if (n <= 0d) return [0d, 1d] as double[]
        double z = 1.645d
        double p = successes / n
        double z2 = z * z
        double denom = 1d + z2 / n
        double centre = (p + z2 / (2d * n)) / denom
        double half = z * Math.sqrt(p * (1d - p) / n + z2 / (4d * n * n)) / denom
        [Math.max(0d, centre - half), Math.min(1d, centre + half)] as double[]
    }
}
