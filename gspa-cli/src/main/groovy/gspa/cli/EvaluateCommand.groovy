package gspa.cli

import gspa.io.GafReader
import gspa.io.GffReader
import gspa.metrics.QualityPipeline
import gspa.metrics.QualityReportWriter
import gspa.model.OrganismDomain
import picocli.CommandLine.Command
import picocli.CommandLine.Option

@Command(name = 'evaluate', description = 'Evaluate quality of existing annotations')
class EvaluateCommand implements Runnable {

    @Option(names = ['-i', '--input'], required = true, description = 'Input GFF3 file')
    File input

    @Option(names = ['-a', '--annotations'], description = 'GO annotations in GAF format')
    File annotations

    @Option(names = ['--go-owl'], required = true, description = 'Path to GO OWL file')
    String goOwl

    @Option(names = ['-k', '--kingdom'], description = 'Organism kingdom for essential functions profile',
            defaultValue = 'bacteria')
    String kingdom

    @Option(names = ['--complex-terms'], description = 'Complex terms classification TSV')
    File complexTerms

    @Option(names = ['--ec2go'], description = 'EC-to-GO mapping file')
    File ec2go

    @Option(names = ['--pathways'], description = 'Pathway definitions TSV')
    File pathways

    @Option(names = ['--taxonomy'], description = 'Taxonomy hierarchy TSV (child\\tparent)')
    File taxonomy

    @Option(names = ['--mag'], description = 'Input is a MAG')
    boolean mag

    @Option(names = ['--lite'], description = 'Skip ELK initialization (faster, no process coherence)')
    boolean lite

    @Option(names = ['--reasoner-cache'],
            description = 'Directory for cached has_part pairs (speeds up full ELK init)')
    File reasonerCache

    @Option(names = ['-o', '--output'], description = 'Output report file', defaultValue = 'quality_report.json')
    File output

    @Override
    void run() {
        println "GSPA evaluate"
        println "  Input: ${input}"
        println "  GO ontology: ${goOwl}"
        println "  Kingdom: ${kingdom}"
        println ""

        // Load genome from GFF3
        def genome = GffReader.readGff3(input)
        genome.mag = mag
        println "Loaded genome: ${genome.contigs.size()} contigs, ${genome.features.size()} features"

        // Load GO annotations from GAF
        if (annotations?.exists()) {
            println "Loading GO annotations from: ${annotations}"
            def goAnnotations = GafReader.readGaf(annotations)

            // Create proteins from CDS features and assign annotations
            genome.contigs.each { contig ->
                contig.featuresOfType(gspa.model.FeatureType.CDS).each { feature ->
                    def protein = new gspa.model.Protein(
                        id: feature.id,
                        sequence: '',
                        sourceFeature: feature
                    )
                    def proteinAnns = goAnnotations[feature.id]
                    if (proteinAnns) {
                        protein.annotations.addAll(proteinAnns)
                    }
                    contig.addProtein(protein)
                }
            }
            println "Assigned annotations to ${genome.proteinCount} proteins " +
                "(${genome.allGoTerms().size()} unique GO terms)"
        }

        // Build quality pipeline
        println ""
        println "Initializing quality pipeline..."
        def pipeline = new QualityPipeline()
            .goOwlFile(goOwl)
            .essentialFunctionsForDomain(OrganismDomain.fromName(kingdom))

        if (complexTerms?.exists()) {
            pipeline.complexTermsFile(complexTerms)
        }
        if (ec2go?.exists() && pathways?.exists()) {
            pipeline.pathwayFile(pathways, ec2go)
        } else if (ec2go?.exists()) {
            pipeline.ec2goFile(ec2go)
        }
        if (taxonomy?.exists()) {
            pipeline.taxonConstraintsFromOntology()
            pipeline.taxonomyFile(taxonomy)
        }

        if (reasonerCache != null) {
            pipeline.reasonerCacheDir(reasonerCache)
        }
        if (lite) {
            pipeline.initializeLite()
        } else {
            pipeline.initialize()
        }

        // Evaluate
        println "Running evaluation..."
        def report = pipeline.evaluate(genome)

        // Write outputs
        def outputDir = output.parentFile ?: new File('.')
        outputDir.mkdirs()

        if (output.name.endsWith('.json')) {
            QualityReportWriter.writeJson(report, output)
        } else {
            QualityReportWriter.writeSummaryTsv([report], output)
        }

        // Also write violations if any
        if (!report.violations.isEmpty()) {
            def violationsFile = new File(outputDir,
                output.name.replaceAll(/\.\w+$/, '_violations.tsv'))
            QualityReportWriter.writeViolationsTsv(report, violationsFile)
            println "Violations written to: ${violationsFile}"
        }

        println ""
        println "=== Quality Report ==="
        println "  Completeness:       ${String.format('%.1f%%', report.completeness * 100)}"
        if (report.processCoherence >= 0) {
            println "  Process coherence:  ${String.format('%.1f%%', report.processCoherence * 100)}"
        }
        if (report.pathwayCoherence >= 0) {
            println "  Pathway coherence:  ${String.format('%.1f%%', report.pathwayCoherence * 100)}"
        }
        if (report.complexCoherence >= 0) {
            println "  Complex coherence:  ${String.format('%.1f%%', report.complexCoherence * 100)}"
        }
        println "  Consistent:         ${report.consistent}"
        println "  Mean IC:            ${String.format('%.2f', report.meanIC)}"
        println "  Composite score:    ${String.format('%.3f', report.compositeScore)}"
        println "  Coverage:           ${report.annotatedProteinCount}/${report.totalProteinCount} proteins"
        if (!report.missingEssentialFunctions.isEmpty()) {
            println "  Missing essential:  ${report.missingEssentialFunctions.size()} functions"
        }
        println ""
        println "Report written to: ${output}"

        pipeline.dispose()
    }
}
