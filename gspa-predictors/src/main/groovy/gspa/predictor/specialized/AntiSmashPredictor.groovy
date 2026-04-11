package gspa.predictor.specialized

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.model.Genome
import gspa.predictor.AbstractToolPredictor
import gspa.predictor.GenomePredictor

/**
 * antiSMASH predictor for biosynthetic gene cluster (BGC) detection.
 * Identifies secondary metabolite clusters: NRPS, PKS, terpenes, etc.
 * Operates at the genome level (needs full nucleotide context).
 */
class AntiSmashPredictor extends AbstractToolPredictor implements GenomePredictor {

    int threads = Runtime.runtime.availableProcessors()

    /** Strictness level: strict, relaxed, loose */
    String strictness = 'relaxed'

    @Override
    String getName() { 'antismash' }

    @Override
    Set<AnnotationType> getOutputTypes() { [AnnotationType.SECONDARY_METABOLITE] as Set }

    @Override
    String getExecutable() { 'antismash' }

    @Override
    List<String> buildCommand(File inputFasta, File outputDir) {
        [
            executablePath ?: executable,
            inputFasta.absolutePath,
            '--output-dir', outputDir.absolutePath,
            '--cpus', threads.toString(),
            '--hmmdetection-strictness', strictness,
            '--minimal',  // skip resource-heavy analyses for speed
            '--genefinding-tool', 'prodigal',
        ]
    }

    @Override
    Map<String, List<Annotation>> parseOutput(File outputDir) {
        // antiSMASH outputs a JSON file with region/cluster information
        def jsonFiles = outputDir.listFiles()?.findAll { it.name.endsWith('.json') && it.name != 'regions.js' }
        if (!jsonFiles) return [:]

        Map<String, List<Annotation>> results = [:].withDefault { [] }

        jsonFiles.each { jsonFile ->
            try {
                def data = new groovy.json.JsonSlurper().parse(jsonFile)
                def records = data.records ?: []

                records.each { record ->
                    def regions = record.areas ?: record.regions ?: []
                    regions.each { region ->
                        String bgcType = region.products?.join(', ') ?: region.product ?: 'unknown'
                        def cdses = region.cds_members ?: region.genes ?: []

                        cdses.each { cds ->
                            String geneId = cds.locus_tag ?: cds.gene_id ?: cds.toString()
                            results[geneId] << new Annotation(
                                type: AnnotationType.SECONDARY_METABOLITE,
                                value: bgcType,
                                source: name,
                                score: 0.8,
                                metadata: [regionType: bgcType]
                            )
                        }
                    }
                }
            } catch (Exception e) {
                log.warn("Failed to parse antiSMASH JSON ${jsonFile}: ${e.message}")
            }
        }
        results
    }

    @Override
    Map<String, List<Annotation>> predictGenome(Genome genome) {
        def tmpDir = File.createTempDir("gspa_antismash_", '')
        try {
            def inputFasta = new File(tmpDir, 'genome.fna')
            inputFasta.withWriter { writer ->
                genome.contigs.each { contig ->
                    if (contig.sequence) {
                        writer.writeLine(">${contig.id}")
                        writer.writeLine(contig.sequence)
                    }
                }
            }
            def outputDir = new File(tmpDir, 'output')
            def command = buildCommand(inputFasta, outputDir)
            def result = execute(command, tmpDir)

            if (result.exitCode != 0) {
                log.error("antiSMASH failed: ${result.stderr}")
                return [:]
            }
            return parseOutput(outputDir)
        } finally {
            tmpDir.deleteDir()
        }
    }
}
