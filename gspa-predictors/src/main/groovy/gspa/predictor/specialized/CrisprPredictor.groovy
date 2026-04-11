package gspa.predictor.specialized

import gspa.model.Annotation
import gspa.model.AnnotationType
import gspa.model.Genome
import gspa.predictor.AbstractToolPredictor
import gspa.predictor.GenomePredictor

/**
 * CRISPR-Cas system predictor using minced (CRISPR finder).
 * Identifies CRISPR arrays and associated Cas genes.
 */
class CrisprPredictor extends AbstractToolPredictor implements GenomePredictor {

    /** Minimum number of repeats to call a CRISPR array */
    int minRepeats = 3

    @Override
    String getName() { 'crispr' }

    @Override
    Set<AnnotationType> getOutputTypes() { [AnnotationType.CRISPR] as Set }

    @Override
    String getExecutable() { 'minced' }

    @Override
    List<String> buildCommand(File inputFasta, File outputDir) {
        def outputFile = new File(outputDir, 'minced_results.gff')
        [
            executablePath ?: executable,
            '-minNR', minRepeats.toString(),
            '-gffFull',
            inputFasta.absolutePath,
            outputFile.absolutePath,
        ]
    }

    @Override
    Map<String, List<Annotation>> parseOutput(File outputDir) {
        def outputFile = new File(outputDir, 'minced_results.gff')
        if (!outputFile.exists()) return [:]

        Map<String, List<Annotation>> results = [:].withDefault { [] }
        String currentContig = null
        int arrayCount = 0

        outputFile.eachLine { line ->
            if (line.startsWith('#') || line.trim().isEmpty()) return
            def fields = line.split('\t')
            if (fields.length < 9) return

            String contigId = fields[0]
            String featureType = fields[2]

            if (featureType == 'repeat_region') {
                arrayCount++
                String arrayId = "CRISPR_${contigId}_${arrayCount}"
                // Assign to the contig as a whole
                results[contigId] << new Annotation(
                    type: AnnotationType.CRISPR,
                    value: "CRISPR array",
                    source: name,
                    score: 0.9,
                    metadata: [
                        start: fields[3] as int,
                        end: fields[4] as int,
                        arrayId: arrayId,
                    ]
                )
            }
        }
        results
    }

    @Override
    Map<String, List<Annotation>> predictGenome(Genome genome) {
        def tmpDir = File.createTempDir("gspa_crispr_", '')
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
            def result = execute(command, tmpDir)

            if (result.exitCode != 0) {
                log.error("minced failed: ${result.stderr}")
                return [:]
            }
            return parseOutput(tmpDir)
        } finally {
            tmpDir.deleteDir()
        }
    }
}
