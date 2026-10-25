package gspa.cli

import gspa.io.GafReader
import gspa.io.GffReader
import gspa.metrics.QualityPipeline
import gspa.metrics.QualityReportWriter
import gspa.model.OrganismDomain
import picocli.CommandLine.Command
import picocli.CommandLine.Option

@Command(name = 'evaluate', description = 'Evaluate quality of existing annotations', mixinStandardHelpOptions = true)
class EvaluateCommand implements Runnable {

    @Option(names = ['-i', '--input'], description = 'Input GFF3 file')
    File input

    @Option(names = ['-a', '--annotations'], description = 'GO annotations in GAF format')
    File annotations

    @Option(names = ['--batch-manifest'], description = 'TSV manifest (gff_path, gaf_path, out_json_path, kingdom, taxon) for shared-memory batch evaluation')
    File batchManifest

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

    @Option(names = ['--module-defs'], description = 'KEGG module definitions JSON with step logic (benchmark/pathway_coherence/fetch_kegg_module_definitions.py); enables evidence-gated module completeness.')
    File moduleDefs

    @Option(names = ['--modules'],
            description = 'Additional pathway TSV(s) to stack on top of --pathways ' +
                          '(KEGG Modules, MetaCyc, BioCyc, …). Comma-separated list.',
            split = ',')
    List<File> moduleFiles = []

    @Option(names = ['--taxonomy'], description = 'Taxonomy hierarchy TSV (child\\tparent)')
    File taxonomy

    @Option(names = ['--taxon'], description = 'Organism NCBI taxon (e.g. NCBITaxon_562, bacteria, eukaryote)')
    String organismTaxon

    @Option(names = ['--mag'], description = 'Input is a MAG')
    boolean mag

    @Option(names = ['--lite'], description = 'Skip ELK initialization (faster, no process coherence)')
    boolean lite

    @Option(names = ['--reasoner-cache'],
            description = 'Directory for cached has_part pairs (speeds up full ELK init)')
    File reasonerCache

    @Option(names = ['-o', '--output'], description = 'Output report file', defaultValue = 'quality_report.json')
    File output

    private gspa.model.Genome loadGenomeWithGaf(File gffFile, File gafFile, String taxStr, String kingStr) {
        def genome = GffReader.readGff3(gffFile)
        genome.mag = mag
        String normTax = taxStr ? gspa.ontology.TaxonConstraints.normalizeTaxon(taxStr) : (kingStr ? gspa.ontology.TaxonConstraints.normalizeTaxon(kingStr) : null)
        if (normTax) {
            genome.taxonomy = new gspa.model.TaxonomyInfo(ncbiTaxId: normTax)
        }
        if (gafFile?.exists()) {
            def goAnnotations = GafReader.readGaf(gafFile)
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
        }
        genome
    }

    @Override
    void run() {
        println "GSPA evaluate"
        println "  GO ontology: ${goOwl}"

        println "Initializing quality pipeline..."
        def pipeline = new QualityPipeline()
            .goOwlFile(goOwl)
            .essentialFunctionsForDomain(OrganismDomain.fromName(kingdom))

        if (complexTerms?.exists()) {
            pipeline.complexTermsFile(complexTerms)
        }
        if (ec2go?.exists() && pathways?.exists()) {
            def ec2goMap = gspa.ontology.PathwayLoader.loadEc2Go(ec2go)
            def pwDb = gspa.ontology.PathwayLoader.loadPathways(pathways, ec2goMap)
            for (File mod : (moduleFiles ?: [])) {
                if (mod != null && mod.exists()) {
                    gspa.ontology.PathwayLoader.loadPathwaysInto(pwDb, mod)
                }
            }
            pipeline.pathwayDatabase(pwDb)
            if (moduleDefs?.exists()) pipeline.moduleDefinitionsFile(moduleDefs, ec2goMap)
        } else if (ec2go?.exists()) {
            pipeline.ec2goFile(ec2go)
            if (moduleDefs?.exists()) pipeline.moduleDefinitionsFile(moduleDefs, gspa.ontology.PathwayLoader.loadEc2Go(ec2go))
        }
        if (taxonomy?.exists()) {
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

        if (batchManifest != null && batchManifest.exists()) {
            List<String> lines = batchManifest.readLines().findAll { it.trim() && !it.startsWith('#') }
            int idx = 0
            for (String line : lines) {
                String[] cols = line.split('\t')
                if (cols.length < 5 || cols[0] == 'gff_path') continue
                idx++
                File gffF = new File(cols[0])
                File gafF = new File(cols[1])
                File outF = new File(cols[2])
                String rowKing = cols[3]
                String rowTax = cols[4]
                pipeline.updateDomain(OrganismDomain.fromName(rowKing))
                def genome = loadGenomeWithGaf(gffF, gafF, rowTax, rowKing)
                def report = pipeline.evaluate(genome)
                def outDir = outF.parentFile ?: new File('.')
                outDir.mkdirs()
                QualityReportWriter.writeJson(report, outF)
                println "[BATCH-EVAL ${idx}] ${outF.name} -> comp=${String.format('%.3f', report.completeness)}, pw=${String.format('%.3f', report.pathwayCoherence)}, cx=${String.format('%.3f', report.complexCoherence)}, viols=${report.violations?.size() ?: 0}, score=${String.format('%.3f', report.compositeScore)}"
            }
            pipeline.dispose()
            return
        }

        if (input == null) {
            throw new IllegalArgumentException("Either -i/--input or --batch-manifest must be specified")
        }

        def genome = loadGenomeWithGaf(input, annotations, organismTaxon, kingdom)
        println "Loaded genome: ${genome.contigs.size()} contigs, ${genome.features.size()} features (taxon=${genome.taxonomy?.ncbiTaxId})"

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
