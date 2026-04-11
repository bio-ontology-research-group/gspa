package gspa.model

import groovy.transform.Canonical
import groovy.transform.builder.Builder
import groovy.transform.builder.SimpleStrategy

/**
 * Assembly quality metadata, particularly relevant for MAGs.
 * Integrates with CheckM/CheckM2 quality estimates.
 */
@Canonical
@Builder(builderStrategy = SimpleStrategy, prefix = '')
class AssemblyInfo {

    /** CheckM/CheckM2 estimated completeness (0-100) */
    Double completeness

    /** CheckM/CheckM2 estimated contamination (0-100) */
    Double contamination

    /** Assembly level: complete, chromosome, scaffold, contig */
    String assemblyLevel

    /** Source: isolate, MAG, SAG, etc. */
    String assemblySource

    /** MIMAG quality tier based on completeness and contamination */
    QualityTier getQualityTier() {
        if (completeness == null || contamination == null) return QualityTier.UNKNOWN
        if (completeness >= 90 && contamination < 5) return QualityTier.HIGH
        if (completeness >= 50 && contamination < 10) return QualityTier.MEDIUM
        return QualityTier.LOW
    }
}
