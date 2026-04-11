package gspa.metrics

import gspa.model.*
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Adjusts quality metrics for metagenome-assembled genomes (MAGs).
 *
 * MAGs are typically incomplete and may contain contamination from
 * other organisms. This adjuster normalizes quality scores using
 * CheckM/CheckM2 estimates and flags violations that may be caused
 * by contamination rather than annotation errors.
 */
class MagAdjuster {

    private static final Logger log = LoggerFactory.getLogger(MagAdjuster)

    /** CheckM2 completeness (0-100) */
    double completeness

    /** CheckM2 contamination (0-100) */
    double contamination

    /** GTDB-Tk lineage (if available) */
    String gtdbLineage

    /** Per-contig taxonomy from GTDB-Tk (contig ID -> lineage) */
    Map<String, String> contigLineage = [:]

    MagAdjuster(AssemblyInfo assemblyInfo) {
        this.completeness = assemblyInfo?.completeness ?: 100.0
        this.contamination = assemblyInfo?.contamination ?: 0.0
    }

    MagAdjuster(double completeness, double contamination) {
        this.completeness = completeness
        this.contamination = contamination
    }

    /** MIMAG quality tier */
    QualityTier getQualityTier() {
        if (completeness >= 90 && contamination < 5) return QualityTier.HIGH
        if (completeness >= 50 && contamination < 10) return QualityTier.MEDIUM
        return QualityTier.LOW
    }

    /**
     * Adjust completeness score: normalize by genome completeness.
     * If a MAG is 70% complete, we expect ~70% of essential functions.
     */
    double adjustCompleteness(double rawScore) {
        if (completeness <= 0) return rawScore
        double adjusted = rawScore / (completeness / 100.0)
        Math.min(adjusted, 1.0)
    }

    /**
     * Adjust coherence score: account for missing genes due to incompleteness.
     * Unsatisfied dependencies may be due to missing contigs, not bad annotations.
     */
    double adjustCoherence(double rawScore, int triggered, int unsatisfied) {
        if (completeness >= 95) return rawScore
        // Estimate how many unsatisfied dependencies are due to missing genes
        double missingFraction = 1.0 - (completeness / 100.0)
        int expectedMissing = (int)(triggered * missingFraction)
        int trueUnsatisfied = Math.max(0, unsatisfied - expectedMissing)
        triggered > 0 ? 1.0 - (trueUnsatisfied / (double) triggered) : 1.0
    }

    /**
     * Assess whether a consistency violation is likely caused by contamination
     * rather than annotation error.
     *
     * Heuristics:
     * - If contamination > 5% AND conflicting annotations come from different contigs,
     *   it's likely contamination
     * - If the organism's GTDB lineage is known, check if the violated taxon constraint
     *   involves a taxon outside the lineage
     */
    boolean isLikelyContamination(ConsistencyViolation violation, Genome genome) {
        if (contamination < 2.0) return false

        // Check if involved proteins are on different contigs
        if (violation.involvedProteins.size() >= 2) {
            Set<String> contigs = []
            violation.involvedProteins.each { proteinId ->
                def protein = genome.findProtein(proteinId)
                if (protein?.contig) {
                    contigs << protein.contig.id
                }
            }
            // Proteins on different contigs → likely contamination
            if (contigs.size() > 1 && contamination > 5.0) {
                return true
            }
        }

        // Check per-contig lineage if available
        if (gtdbLineage && !contigLineage.isEmpty()) {
            def involvedContigs = violation.involvedProteins.collect { pid ->
                genome.findProtein(pid)?.contig?.id
            }.findAll { it != null }.unique()

            return involvedContigs.any { contigId ->
                def lineage = contigLineage[contigId]
                lineage && lineage != gtdbLineage
            }
        }

        false
    }

    /**
     * Run per-contig quality assessment.
     * Identifies contigs that may be contamination based on:
     * - Anomalous GC content compared to genome average
     * - Different taxonomic assignments
     * - Unusual coding density
     */
    List<ContigQuality> assessContigs(Genome genome) {
        if (genome.contigs.isEmpty()) return []

        // Compute genome-wide GC
        double genomeGC = computeGC(genome.contigs.collectMany {
            it.sequence ? [it.sequence] : []
        }.join(''))

        genome.contigs.collect { contig ->
            double contigGC = contig.sequence ? computeGC(contig.sequence) : genomeGC
            double gcDeviation = Math.abs(contigGC - genomeGC)
            double codingDensity = contig.proteins.size() > 0 && contig.length > 0 ?
                contig.proteins.sum { Protein p -> p.end - p.start + 1 } / (double) contig.length : 0.0

            boolean suspect = gcDeviation > 0.10 // >10% GC deviation is suspicious
            String lineage = contigLineage[contig.id]
            if (lineage && gtdbLineage && lineage != gtdbLineage) {
                suspect = true
            }

            new ContigQuality(
                contigId: contig.id,
                length: contig.length,
                gcContent: contigGC,
                gcDeviation: gcDeviation,
                codingDensity: codingDensity,
                proteinCount: contig.proteins.size(),
                suspectContamination: suspect,
                lineage: lineage
            )
        }
    }

    /**
     * Build a full MAG quality report augmenting the standard QualityReport.
     */
    MagQualityReport buildMagReport(QualityReport baseReport, Genome genome) {
        def contigAssessment = assessContigs(genome)
        int suspectContigs = contigAssessment.count { it.suspectContamination }

        // Reclassify violations
        def reclassified = baseReport.violations.collect { v ->
            def vCopy = new ConsistencyViolation(
                type: v.type, severity: v.severity,
                description: v.description,
                involvedProteins: v.involvedProteins,
                involvedGoTerms: v.involvedGoTerms,
                suggestedAction: v.suggestedAction,
                justification: v.justification
            )
            if (isLikelyContamination(v, genome)) {
                vCopy.severity = ConsistencyViolation.Severity.WARNING
                vCopy.description = "[Likely contamination] ${v.description}"
                vCopy.suggestedAction = "Check if involved contigs are contamination; consider removing from bin"
            }
            vCopy
        }

        new MagQualityReport(
            baseReport: baseReport,
            qualityTier: qualityTier,
            checkMCompleteness: completeness,
            checkMContamination: contamination,
            adjustedCompleteness: adjustCompleteness(baseReport.completeness),
            contigAssessment: contigAssessment,
            suspectContigCount: suspectContigs,
            reclassifiedViolations: reclassified,
            gtdbLineage: gtdbLineage
        )
    }

    private static double computeGC(String seq) {
        if (!seq) return 0.0
        int gc = 0; int total = 0
        seq.toUpperCase().each { ch ->
            if (ch == 'G' || ch == 'C') gc++
            if (ch in ['A', 'T', 'G', 'C']) total++
        }
        total > 0 ? gc / (double) total : 0.0
    }
}

class ContigQuality {
    String contigId
    int length
    double gcContent
    double gcDeviation
    double codingDensity
    int proteinCount
    boolean suspectContamination
    String lineage
}

class MagQualityReport {
    QualityReport baseReport
    QualityTier qualityTier
    double checkMCompleteness
    double checkMContamination
    double adjustedCompleteness
    List<ContigQuality> contigAssessment = []
    int suspectContigCount
    List<ConsistencyViolation> reclassifiedViolations = []
    String gtdbLineage
}
