package gspa.metrics

import gspa.model.ConsistencyViolation
import gspa.model.ConsistencyViolation.Severity
import gspa.model.ConsistencyViolation.ViolationType
import gspa.model.QualityReport
import spock.lang.Specification

class QualityReportWriterSpec extends Specification {

    def "should write JSON report"() {
        given:
        def report = buildTestReport()
        def tmpFile = File.createTempFile('report_', '.json')

        when:
        QualityReportWriter.writeJson(report, tmpFile)
        def json = tmpFile.text

        then:
        json.contains('"genome_id"')
        json.contains('"test_genome"')
        json.contains('"completeness"')
        json.contains('"coherence"')
        json.contains('"consistency"')

        cleanup:
        tmpFile.delete()
    }

    def "should write summary TSV"() {
        given:
        def reports = [buildTestReport(), buildTestReport().tap { genomeId = 'genome2' }]
        def tmpFile = File.createTempFile('summary_', '.tsv')

        when:
        QualityReportWriter.writeSummaryTsv(reports, tmpFile)
        def lines = tmpFile.readLines()

        then:
        lines.size() == 3  // header + 2 data lines
        lines[0].startsWith('genome_id\t')
        lines[1].startsWith('test_genome\t')
        lines[2].startsWith('genome2\t')

        cleanup:
        tmpFile.delete()
    }

    def "should write violations TSV"() {
        given:
        def report = buildTestReport()
        report.violations << new ConsistencyViolation(
            type: ViolationType.TAXON_CONFLICT,
            severity: Severity.ERROR,
            description: 'Test violation',
            involvedProteins: ['p1', 'p2'],
            involvedGoTerms: ['GO:0006412'],
        )
        def tmpFile = File.createTempFile('violations_', '.tsv')

        when:
        QualityReportWriter.writeViolationsTsv(report, tmpFile)
        def lines = tmpFile.readLines()

        then:
        lines.size() == 2  // header + 1 violation
        lines[1].contains('TAXON_CONFLICT')
        lines[1].contains('p1,p2')

        cleanup:
        tmpFile.delete()
    }

    def "should generate valid JSON string"() {
        given:
        def report = buildTestReport()

        when:
        def json = QualityReportWriter.toJson(report)

        then:
        json.contains('test_genome')
        // Should be valid JSON (no exception)
        new groovy.json.JsonSlurper().parseText(json) != null
    }

    private QualityReport buildTestReport() {
        new QualityReport(
            genomeId: 'test_genome',
            assessmentDate: '2026-04-10',
            completeness: 0.85,
            processCoherence: 0.92,
            pathwayCoherence: 0.78,
            complexCoherence: 1.0,
            consistent: true,
            meanIC: 6.5,
            annotatedProteinCount: 400,
            totalProteinCount: 500,
            presentEssentialFunctions: ['GO:0006412', 'GO:0006260'] as Set,
            missingEssentialFunctions: ['GO:0051301'] as Set,
            annotationCountBySource: ['diamond': 350, 'interproscan': 280]
        )
    }
}
