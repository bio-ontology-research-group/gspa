package gspa.cli

import gspa.io.GafReader
import gspa.io.GffReader
import gspa.metrics.HtmlReportWriter
import gspa.metrics.QualityPipeline
import gspa.metrics.QualityReportWriter
import gspa.model.OrganismDomain
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

@Command(name = 'compare', description = 'Compare quality of multiple annotation sets for the same genome')
class CompareCommand implements Runnable {

    @Option(names = ['-i', '--input'], required = true, description = 'Input GFF3 file (genome structure)')
    File input

    @Option(names = ['--go-owl'], required = true, description = 'Path to GO OWL file')
    String goOwl

    @Option(names = ['-k', '--kingdom'], defaultValue = 'bacteria')
    String kingdom

    @Option(names = ['--lite'], description = 'Skip ELK (no process coherence)')
    boolean lite

    @Option(names = ['-o', '--output'], description = 'Output report', defaultValue = 'comparison.html')
    File output

    @Parameters(description = 'GAF annotation files to compare (one per method)')
    List<File> gafFiles

    @Override
    void run() {
        if (!gafFiles || gafFiles.isEmpty()) {
            println "Error: provide at least one GAF annotation file as positional argument."
            return
        }

        println "GSPA compare: ${gafFiles.size()} annotation sets"
        println ""

        // Build quality pipeline once
        def pipeline = new QualityPipeline()
            .goOwlFile(goOwl)
            .essentialFunctionsForDomain(OrganismDomain.fromName(kingdom))
        if (lite) {
            pipeline.initializeLite()
        } else {
            pipeline.initialize()
        }

        // Evaluate each annotation set
        List<gspa.model.QualityReport> reports = []
        gafFiles.each { gafFile ->
            String methodName = gafFile.name.replaceAll(/\.(gaf|tsv)(\.gz)?$/, '')
            println "Evaluating: ${methodName}..."

            // Load genome fresh for each method
            def genome = GffReader.readGff3(input)
            def goAnnotations = GafReader.readGaf(gafFile, methodName)

            genome.contigs.each { contig ->
                contig.featuresOfType(gspa.model.FeatureType.CDS).each { feature ->
                    def protein = new gspa.model.Protein(id: feature.id, sequence: '', sourceFeature: feature)
                    goAnnotations[feature.id]?.each { protein.annotations.add(it) }
                    contig.addProtein(protein)
                }
            }

            def report = pipeline.evaluate(genome)
            report.genomeId = methodName
            reports << report

            println "  Completeness: ${String.format('%.1f%%', report.completeness * 100)}, " +
                "Composite: ${String.format('%.3f', report.compositeScore)}"
        }

        // Write comparison report
        if (output.name.endsWith('.html')) {
            HtmlReportWriter.writeMultiGenomeHtml(reports, output)
        } else if (output.name.endsWith('.json')) {
            reports.each { r ->
                def f = new File(output.parentFile ?: new File('.'), "${r.genomeId}_quality.json")
                QualityReportWriter.writeJson(r, f)
            }
        } else {
            QualityReportWriter.writeSummaryTsv(reports, output)
        }

        println ""
        println "Comparison written to: ${output}"

        pipeline.dispose()
    }
}
