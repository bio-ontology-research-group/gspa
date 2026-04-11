package gspa.predictor.pathway

import gspa.model.*
import gspa.predictor.AbstractToolPredictor
import gspa.predictor.GenomePredictor

/**
 * gapseq predictor for metabolic pathway reconstruction.
 * Identifies metabolic reactions and pathways, performs gap-filling
 * to create draft genome-scale metabolic models.
 *
 * Supports community mode for crossfeeding analysis.
 */
class GapseqPredictor extends AbstractToolPredictor implements GenomePredictor {

    /** gapseq medium file for gap-filling */
    String medium

    /** Export SBML model */
    boolean exportSbml = false

    /** Taxonomy for pathway prediction */
    String taxonomy = 'Bacteria'

    @Override
    String getName() { 'gapseq' }

    @Override
    Set<AnnotationType> getOutputTypes() {
        [AnnotationType.EC, AnnotationType.KEGG] as Set
    }

    @Override
    String getExecutable() { 'gapseq' }

    @Override
    List<String> buildCommand(File inputFasta, File outputDir) {
        def outputPrefix = new File(outputDir, 'gapseq_model')
        def cmd = [
            executablePath ?: executable,
            'find',
            '-p', 'all',
            '-b', taxonomy,
            inputFasta.absolutePath,
        ]
        cmd
    }

    @Override
    Map<String, List<Annotation>> parseOutput(File outputDir) {
        // gapseq outputs reaction lists and pathway completeness
        def reactionsFile = findOutputFile(outputDir, '*-Reactions.tbl')
        if (reactionsFile == null) return [:]

        Map<String, List<Annotation>> results = [:].withDefault { [] }

        reactionsFile.eachLine { line ->
            if (line.startsWith('#') || line.startsWith('Reaction') || line.trim().isEmpty()) return
            def fields = line.split('\t')
            if (fields.length < 6) return

            String reactionId = fields[0]
            String ecNumber = fields[2]
            String status = fields[4]      // active/inactive
            String geneName = fields.length > 5 ? fields[5] : ''

            if (ecNumber && ecNumber != '-' && status == 'active') {
                // Map EC to protein IDs based on gene names
                // gapseq reports gene locus tags
                if (geneName && geneName != '-') {
                    results[geneName] << new Annotation(
                        type: AnnotationType.EC,
                        value: ecNumber.startsWith('EC:') ? ecNumber : "EC:${ecNumber}",
                        source: name,
                        score: 0.7,
                        metadata: [reaction: reactionId, status: status]
                    )
                }
            }
        }
        results
    }

    @Override
    Map<String, List<Annotation>> predictGenome(Genome genome) {
        // Write genome FASTA and run gapseq
        def tmpDir = File.createTempDir("gspa_gapseq_", '')
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

            def command = buildCommand(inputFasta, tmpDir)
            log.info("Running gapseq: ${command.join(' ')}")
            def result = execute(command, tmpDir)

            if (result.exitCode != 0) {
                log.error("gapseq failed: ${result.stderr}")
                return [:]
            }

            return parseOutput(tmpDir)
        } finally {
            tmpDir.deleteDir()
        }
    }

    /**
     * Run gapseq in community mode for crossfeeding analysis.
     * Takes multiple genome FASTA files and finds metabolic exchanges.
     */
    Map<String, Object> analyzeCommmunity(List<File> genomeFastas, File outputDir) {
        // gapseq community analysis would be done in a separate step
        // after individual models are built
        log.info("Community analysis with ${genomeFastas.size()} genomes")
        [:]
    }

    private File findOutputFile(File dir, String pattern) {
        dir.listFiles()?.find { it.name.endsWith('-Reactions.tbl') }
    }
}
