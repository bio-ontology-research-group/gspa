package gspa.cli

import gspa.config.ConfigLoader
import gspa.config.EssentialFunctions
import gspa.config.GspaConfig
import gspa.io.FastaReader
import gspa.io.GafReader
import gspa.io.GffReader
import gspa.io.GffWriter
import gspa.metrics.HtmlReportWriter
import gspa.metrics.QualityPipeline
import gspa.metrics.QualityReportWriter
import gspa.model.OrganismDomain
import gspa.predictor.AnnotationPipeline
import picocli.CommandLine
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

/**
 * GSPA - Genome-Scale Protein Annotation
 *
 * A comprehensive genome annotation pipeline combining multiple prediction
 * strategies with built-in quality assessment (completeness, coherence, consistency).
 */
@Command(
    name = 'gspa',
    description = 'Genome-Scale Protein Annotation pipeline',
    version = 'gspa 0.1.0-SNAPSHOT',
    mixinStandardHelpOptions = true,
    subcommands = [
        AnnotateCommand,
        EvaluateCommand,
        CompareCommand,
        ReportCommand,
        IntegrateCommand,
    ]
)
class GspaMain implements Runnable {

    @Override
    void run() {
        CommandLine.usage(this, System.out)
    }

    static void main(String[] args) {
        int exitCode = new CommandLine(new GspaMain()).execute(args)
        System.exit(exitCode)
    }
}

@Command(name = 'annotate', description = 'Annotate a genome or set of genomes')
class AnnotateCommand implements Runnable {

    @Option(names = ['-i', '--input'], required = true, description = 'Input FASTA/GenBank file')
    File input

    @Option(names = ['-o', '--output'], description = 'Output directory', defaultValue = 'gspa_output')
    File outputDir

    @Option(names = ['-c', '--config'], description = 'Configuration YAML file')
    File configFile

    @Option(names = ['-k', '--kingdom'], description = 'Organism kingdom: bacteria, archaea, eukaryote, virus, auto')
    String kingdom

    @Option(names = ['--mag'], description = 'Input is a MAG (adjust quality thresholds)')
    boolean mag

    @Option(names = ['--go-owl'], description = 'Path to GO OWL file for quality assessment')
    String goOwl

    @Option(names = ['-t', '--threads'], description = 'Number of threads', defaultValue = '0')
    int threads

    @Option(names = ['--db', '--database'], description = 'Path to sequence similarity database')
    String database

    @Override
    void run() {
        println "GSPA annotate: ${input}"
        println "Output: ${outputDir}"
        println "Kingdom: ${kingdom ?: 'auto'}"
        println "MAG mode: ${mag}"
        println ""

        def overrides = [:]
        if (kingdom) overrides['kingdom'] = kingdom
        if (goOwl) overrides['go-owl'] = goOwl
        if (database) overrides['database'] = database
        def config = ConfigLoader.buildConfig(configFile, kingdom, overrides)
        if (mag) config.input.type = 'mag'

        // Build and configure annotation pipeline
        def pipeline = new AnnotationPipeline(config)
        pipeline.configure()

        // Run annotation
        def genome = pipeline.annotate(input)

        // Write output
        pipeline.writeOutput(genome, outputDir)

        // Run quality evaluation if GO ontology provided
        if (goOwl && genome.allGoTerms().size() > 0) {
            println ""
            println "Running quality evaluation..."
            def qualityPipeline = new QualityPipeline(config)
                .goOwlFile(goOwl)
                .essentialFunctionsForDomain(config.resolveOrganismDomain())
                .initializeLite()

            def report = qualityPipeline.evaluate(genome)
            def reportFile = new File(outputDir, "${genome.id}_quality.json")
            QualityReportWriter.writeJson(report, reportFile)
            println "Quality report: ${reportFile}"

            qualityPipeline.dispose()
        }

        println ""
        println "Done. Output in: ${outputDir}"
    }
}

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
