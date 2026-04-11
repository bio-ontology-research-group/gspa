package gspa.metrics

import com.fasterxml.jackson.databind.ObjectMapper
import com.fasterxml.jackson.databind.SerializationFeature
import gspa.model.QualityReport
import gspa.model.ConsistencyViolation
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Writes quality assessment reports in JSON and TSV formats.
 */
class QualityReportWriter {

    private static final Logger log = LoggerFactory.getLogger(QualityReportWriter)

    private static final ObjectMapper mapper = new ObjectMapper()
        .enable(SerializationFeature.INDENT_OUTPUT)

    /**
     * Write a quality report as JSON.
     */
    static void writeJson(QualityReport report, File output) {
        def reportMap = buildReportMap(report)
        mapper.writeValue(output, reportMap)
        log.info("Quality report written to: ${output}")
    }

    /**
     * Write a quality report as JSON string.
     */
    static String toJson(QualityReport report) {
        mapper.writeValueAsString(buildReportMap(report))
    }

    /**
     * Write a summary TSV (one line per genome, suitable for multi-genome comparison).
     */
    static void writeSummaryTsv(List<QualityReport> reports, File output) {
        output.withWriter { writer ->
            // Header
            writer.writeLine([
                'genome_id',
                'assessment_date',
                'total_proteins',
                'annotated_proteins',
                'annotation_coverage',
                'completeness',
                'process_coherence',
                'pathway_coherence',
                'complex_coherence',
                'consistent',
                'violation_count',
                'mean_ic',
                'composite_score',
                'present_essential',
                'missing_essential'
            ].join('\t'))

            reports.each { report ->
                writer.writeLine([
                    report.genomeId,
                    report.assessmentDate,
                    report.totalProteinCount,
                    report.annotatedProteinCount,
                    fmt(report.annotationCoverage),
                    fmt(report.completeness),
                    fmt(report.processCoherence),
                    fmt(report.pathwayCoherence),
                    fmt(report.complexCoherence),
                    report.consistent,
                    report.violations.size(),
                    fmt(report.meanIC),
                    fmt(report.compositeScore),
                    report.presentEssentialFunctions.size(),
                    report.missingEssentialFunctions.size()
                ].join('\t'))
            }
        }
        log.info("Summary TSV written to: ${output}")
    }

    /**
     * Write detailed violation report as TSV.
     */
    static void writeViolationsTsv(QualityReport report, File output) {
        output.withWriter { writer ->
            writer.writeLine([
                'genome_id',
                'violation_type',
                'severity',
                'description',
                'involved_proteins',
                'involved_go_terms',
                'suggested_action',
                'justification'
            ].join('\t'))

            report.violations.each { v ->
                writer.writeLine([
                    report.genomeId,
                    v.type,
                    v.severity,
                    v.description,
                    v.involvedProteins.join(','),
                    v.involvedGoTerms.join(','),
                    v.suggestedAction ?: '',
                    v.justification ?: ''
                ].join('\t'))
            }
        }
    }

    /**
     * Write missing essential functions as TSV (actionable for gap analysis).
     */
    static void writeMissingEssentialsTsv(QualityReport report, File output) {
        output.withWriter { writer ->
            writer.writeLine(['genome_id', 'missing_go_term', 'status'].join('\t'))

            report.missingEssentialFunctions.each { term ->
                writer.writeLine([report.genomeId, term, 'MISSING'].join('\t'))
            }
            report.presentEssentialFunctions.each { term ->
                writer.writeLine([report.genomeId, term, 'PRESENT'].join('\t'))
            }
        }
    }

    private static Map buildReportMap(QualityReport report) {
        [
            genome_id              : report.genomeId,
            assessment_date        : report.assessmentDate,
            summary                : [
                total_proteins      : report.totalProteinCount,
                annotated_proteins  : report.annotatedProteinCount,
                annotation_coverage : report.annotationCoverage,
                composite_score     : report.compositeScore,
            ],
            completeness           : [
                score               : report.completeness,
                present_count       : report.presentEssentialFunctions.size(),
                missing_count       : report.missingEssentialFunctions.size(),
                present_functions   : report.presentEssentialFunctions.sort(),
                missing_functions   : report.missingEssentialFunctions.sort(),
            ],
            coherence              : [
                process_coherence   : report.processCoherence,
                pathway_coherence   : report.pathwayCoherence,
                complex_coherence   : report.complexCoherence,
            ],
            consistency            : [
                consistent          : report.consistent,
                violation_count     : report.violations.size(),
                violations          : report.violations.collect { v ->
                    [
                        type             : v.type.toString(),
                        severity         : v.severity.toString(),
                        description      : v.description,
                        involved_proteins: v.involvedProteins,
                        involved_go_terms: v.involvedGoTerms,
                        suggested_action : v.suggestedAction,
                        justification    : v.justification,
                    ]
                },
            ],
            information_content    : [
                mean_ic             : report.meanIC,
            ],
            annotation_sources     : report.annotationCountBySource,
        ]
    }

    private static String fmt(double v) {
        if (v < 0) return 'N/A'
        String.format('%.4f', v)
    }
}
