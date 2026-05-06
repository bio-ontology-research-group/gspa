package gspa.cli

import gspa.metrics.HtmlReportWriter
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

@Command(name = 'report', description = 'Generate HTML quality report from existing JSON reports')
class ReportCommand implements Runnable {

    @Option(names = ['-o', '--output'], description = 'Output HTML file', defaultValue = 'quality_report.html')
    File output

    @Parameters(description = 'Input JSON quality report files')
    List<File> jsonFiles

    @Override
    void run() {
        if (!jsonFiles || jsonFiles.isEmpty()) {
            println "Error: provide at least one JSON quality report file."
            return
        }

        def mapper = new com.fasterxml.jackson.databind.ObjectMapper()
        List<gspa.model.QualityReport> reports = []

        jsonFiles.each { jsonFile ->
            println "Loading: ${jsonFile}"
            def map = mapper.readValue(jsonFile, Map)

            def report = new gspa.model.QualityReport(
                genomeId: map.genome_id ?: jsonFile.name,
                assessmentDate: map.assessment_date,
                completeness: map.completeness?.score ?: 0.0,
                processCoherence: map.coherence?.process_coherence ?: -1.0,
                pathwayCoherence: map.coherence?.pathway_coherence ?: -1.0,
                complexCoherence: map.coherence?.complex_coherence ?: -1.0,
                consistent: map.consistency?.consistent ?: true,
                meanIC: map.information_content?.mean_ic ?: 0.0,
                totalProteinCount: map.summary?.total_proteins ?: 0,
                annotatedProteinCount: map.summary?.annotated_proteins ?: 0,
                presentEssentialFunctions: (map.completeness?.present_functions ?: []) as Set,
                missingEssentialFunctions: (map.completeness?.missing_functions ?: []) as Set,
            )
            reports << report
        }

        if (reports.size() == 1) {
            HtmlReportWriter.writeHtml(reports[0], output)
        } else {
            HtmlReportWriter.writeMultiGenomeHtml(reports, output)
        }

        println "HTML report written to: ${output}"
    }
}
