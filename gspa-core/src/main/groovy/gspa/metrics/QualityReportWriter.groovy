package gspa.metrics

import com.fasterxml.jackson.databind.ObjectMapper
import com.fasterxml.jackson.databind.SerializationFeature
import gspa.model.QualityReport
import gspa.model.ConsistencyViolation
import gspa.ontology.GoOntology
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
     *
     * Pass {@code goOntology} to enrich the JSON with human-readable GO term
     * names alongside every ID. Without it, only IDs are emitted.
     */
    static void writeJson(QualityReport report, File output, GoOntology goOntology = null) {
        def reportMap = buildReportMap(report, goOntology)
        mapper.writeValue(output, reportMap)
        log.info("Quality report written to: ${output}")
    }

    /**
     * Write a quality report as JSON string.
     */
    static String toJson(QualityReport report, GoOntology goOntology = null) {
        mapper.writeValueAsString(buildReportMap(report, goOntology))
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

    private static Map buildReportMap(QualityReport report, GoOntology goOntology = null) {
        // Helper: turn a GO ID into {id, name} (or just {id} if no name).
        // Looks at the explicit goOntology first (for callers that want to
        // override), then falls back to the labels QualityScorer attached
        // to the report.
        def named = { String id ->
            String name = goOntology?.getLabel(id) ?: report.goLabels?.get(id)
            name ? [id: id, name: name] : [id: id]
        }
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
                score                       : report.completeness,
                present_count               : report.presentEssentialFunctions.size(),
                missing_count               : report.missingEssentialFunctions.size(),
                // Backwards-compat: bare-ID lists kept under their original keys.
                present_functions           : report.presentEssentialFunctions.sort(),
                missing_functions           : report.missingEssentialFunctions.sort(),
                // New: ID + human-readable name.
                present_functions_named     : report.presentEssentialFunctions.sort().collect(named),
                missing_functions_named     : report.missingEssentialFunctions.sort().collect(named),
            ],
            coherence              : [
                process_coherence            : report.processCoherence,
                pathway_coherence            : report.pathwayCoherence,
                complex_coherence            : report.complexCoherence,
                // Per-pair detail: required vs missing GO term, with names.
                process_unsatisfied_pairs    : report.incoherentProcessPairs.collect { p ->
                    [required: named(p.key), missing: named(p.value)]
                },
                // Per-pathway detail (sorted least-coherent first by the metric).
                pathway_detail               : report.incoherentPathways.collect { pw ->
                    [
                        id            : pw.id,
                        name          : pw.name,
                        pathway_term  : pw.pathway_term ? named(pw.pathway_term as String) : null,
                        completeness  : pw.completeness,
                        required      : pw.required,
                        n_present     : (pw.present as List).size(),
                        n_missing     : (pw.missing as List).size(),
                        present_terms : (pw.present as List).collect(named),
                        missing_terms : (pw.missing as List).collect(named),
                    ]
                },
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
                        involved_terms_named: (v.involvedGoTerms ?: []).collect(named),
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
