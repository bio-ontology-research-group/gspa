package gspa.integration.suggester

import gspa.model.GenomeLayout
import gspa.model.ProteinLocus
import groovy.transform.CompileStatic

/**
 * Kernel-smoothed genomic density over a set of weighted anchor
 * proteins. Used by the Phase 12 Reaction-Local Context Suggester to
 * replace the pathway × operon × BF-gate machinery with a continuous
 * "is this locus enriched for reaction-R neighbours?" signal.
 *
 * <p>The kernel is a boxcar-normalised Gaussian
 * {@code K(u) = exp(-u^2/2) / sqrt(2π)} with standard deviation equal
 * to the bandwidth {@code h} (in bp).</p>
 *
 * <p>This class is stateless wrt the genome; all lookups are
 * per-call. In practice a caller constructs one density field per
 * (gap, anchor-set) tuple and queries it at ~20-30 gene positions.</p>
 */
@CompileStatic
class GenomicDensityField {

    /** Protein-id → anchor weight. */
    final Map<String, Double> anchorWeights

    /** Bandwidth in bp. */
    final double bandwidth

    /** Search halfwidth (in bp) for kernel sums. */
    final int halfWidth

    GenomicDensityField(Map<String, Double> anchorWeights, double bandwidth) {
        this.anchorWeights = anchorWeights
        this.bandwidth = bandwidth
        // 4σ truncation covers >99.99% of kernel mass; cheaper than full sum.
        this.halfWidth = (int) Math.max(1.0d, 4.0d * bandwidth)
    }

    /**
     * Compute {@code D(x, R)} at genomic position {@code (contig, midpoint)}.
     * Anchors outside {@code ±halfWidth} contribute negligibly and are skipped.
     */
    double densityAt(GenomeLayout layout, String contig, int midpoint) {
        if (anchorWeights.isEmpty()) return 0.0d
        List<ProteinLocus> nearby = layout.windowAround(contig, midpoint, halfWidth)
        if (nearby.isEmpty()) return 0.0d
        double sum = 0.0d
        double inv2h2 = 1.0d / (2.0d * bandwidth * bandwidth)
        for (ProteinLocus loc : nearby) {
            Double w = anchorWeights[loc.proteinId]
            if (w == null || w == 0.0d) continue
            double dx = (double) (loc.midpoint() - midpoint)
            sum += w * Math.exp(-dx * dx * inv2h2)
        }
        sum
    }

    /** Convenience: density at a specific protein's locus. */
    double densityAt(ProteinLocus locus) {
        densityAt(null, locus.contig, locus.midpoint())
    }

    /**
     * Density at {@code locus} using an already-known {@link GenomeLayout}.
     */
    double densityAt(GenomeLayout layout, ProteinLocus locus) {
        densityAt(layout, locus.contig, locus.midpoint())
    }
}
